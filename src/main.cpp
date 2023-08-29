#include "AlignmentResult.hpp"
#include "LibraryInfo.hpp"
#include "PancakeAligner.hpp"
#include "io/BamZmwReader.hpp"
#include "io/BamZmwReaderConfig.hpp"
#include "io/ZmwRecords.hpp"

#include <htslib/hts.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>
#include <pbbam/DataSet.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastaWriter.h>
#include <pbbam/PbbamVersion.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/PbiRawData.h>
#include <pbcopper/cli2/CLI.h>
#include <pbcopper/cli2/internal/BuiltinOptions.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/parallel/WorkQueue.h>
#include <pbcopper/utility/MemoryConsumption.h>
#include <pbcopper/utility/PbcopperVersion.h>
#include <pbcopper/utility/Ssize.h>
#include <pbcopper/utility/Stopwatch.h>
#include <zlib.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/version.hpp>

#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>

namespace PacBio {
namespace OptionNames {
// clang-format off
const CLI_v2::Option CcsQuery {
R"({
    "names" : ["ccs-query"],
    "description" : "second file is ccs data",
    "type" : "bool",
    "hidden"  : true
})"
};
// clang-format on
}  // namespace OptionNames
struct ActcSettings
{
    std::string InputCLRFile;
    std::string InputCCSFile;
    std::string OutputAlignmentFile;
    int32_t NumThreads = 1;
    int32_t ChunkCur = -1;
    int32_t ChunkAll = -1;
    bool ccsQuery = false;
};

CLI_v2::Interface CreateCLI()
{
    static const std::string description{"Align clr to ccs reads."};
    CLI_v2::Interface i{"actc", description, Actc::LibraryInfo().Release};

    Logging::LogConfig logConfig;
    logConfig.Header = "| ";
    logConfig.Delimiter = " | ";
    logConfig.Fields = Logging::LogField::TIMESTAMP | Logging::LogField::LOG_LEVEL;
    logConfig.Level = Logging::LogLevel::INFO;
    logConfig.LeftBlockWidth = 14;
    logConfig.AlignLevel = true;
    i.LogConfig(logConfig);

    const CLI_v2::PositionalArgument InputCLRFile{
        R"({
        "name" : "IN.subreads.bam",
        "description" : "Subreads BAM.",
        "type" : "file",
        "required" : true
    })"};
    const CLI_v2::PositionalArgument InputCCSFile{
        R"({
        "name" : "IN.ccs.bam",
        "description" : "CCS BAM.",
        "type" : "file",
        "required" : true
    })"};
    const CLI_v2::PositionalArgument Output{
        R"({
        "name" : "OUT.bam",
        "description" : "Aligned subreads to CCS BAM.",
        "type" : "file",
        "required" : true
    })"};
    i.AddPositionalArguments({InputCLRFile, InputCCSFile, Output});
    i.AddOption(IO::OptionNames::Chunk);

    i.AddOption(OptionNames::CcsQuery);

    const auto printVersion = [](const CLI_v2::Interface& interface) {
        const std::string actcVersion = []() {
            return Actc::LibraryInfo().Release + " (commit " + Actc::LibraryInfo().GitSha1 + ')';
        }();
        const std::string pbbamVersion = []() { return BAM::LibraryFormattedVersion(); }();
        const std::string pbcopperVersion = []() {
            return Utility::LibraryVersionString() + " (commit " + Utility::LibraryGitSha1String() +
                   ')';
        }();
        const std::string boostVersion = []() {
            std::string v = BOOST_LIB_VERSION;
            boost::replace_all(v, "_", ".");
            return v;
        }();
        const std::string htslibVersion = []() { return std::string{hts_version()}; }();
        const std::string zlibVersion = []() { return std::string{ZLIB_VERSION}; }();

        std::cout << interface.ApplicationName() << " " << interface.ApplicationVersion() << '\n';
        std::cout << '\n';
        std::cout << "Using:\n";
        std::cout << "  actc     : " << actcVersion << '\n';
        std::cout << "  pbbam    : " << pbbamVersion << '\n';
        std::cout << "  pbcopper : " << pbcopperVersion << '\n';
        std::cout << "  boost    : " << boostVersion << '\n';
        std::cout << "  htslib   : " << htslibVersion << '\n';
        std::cout << "  zlib     : " << zlibVersion << '\n';
    };
    i.RegisterVersionPrinter(printVersion);

    return i;
}

void WorkerThread(Parallel::WorkQueue<std::vector<BAM::BamRecord>>& queue, BAM::BamWriter& writer,
                  const int32_t numReads)
{
    int32_t counter = 0;
    double perc = 0;

    auto LambdaWorker = [&](std::vector<BAM::BamRecord>&& ps) {
        ++counter;
        if (1.0 * counter / numReads > (perc + 0.001)) {
            perc = counter * 1.0 / numReads;
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(2) << 100 * perc << '%';
            PBLOG_BLOCK_INFO("Progress", ss.str());
        }
        for (const auto& record : ps) {
            writer.Write(record);
        }
    };

    while (queue.ConsumeWith(LambdaWorker)) {
    }
}

void SetBamReaderDecompThreads(const int32_t numThreads)
{
    static constexpr char BAMREADER_ENV[] = "PB_BAMREADER_THREADS";
    const std::string decompThreads = std::to_string(numThreads);
    setenv(BAMREADER_ENV, decompThreads.c_str(), true);
}

int RunnerSubroutine(const CLI_v2::Results& options)
{
    Utility::Stopwatch globalTimer;
    SetBamReaderDecompThreads(options.NumThreads());

    ActcSettings settings;
    settings.ccsQuery = options[OptionNames::CcsQuery];
    const std::vector<std::string> files = options.PositionalArguments();

    const auto ReadType = [&settings](const std::string& inputFile) {
        const auto bamFiles = BAM::DataSet(inputFile).BamFiles();
        if (bamFiles.empty()) {
            PBLOG_BLOCK_FATAL("Input checker", "No BAM files available for: " + inputFile);
            std::exit(EXIT_FAILURE);
        }

        std::string readType;
        for (int32_t i = 0; i < Utility::Ssize(bamFiles); ++i) {
            for (const auto& rg : bamFiles.at(i).Header().ReadGroups()) {
                if (readType.empty()) {
                    readType = rg.ReadType();
                } else if (readType != rg.ReadType()) {
                    PBLOG_BLOCK_FATAL(
                        "Input checker",
                        "Do not mix and match different read types for input file : " + inputFile);
                    std::exit(EXIT_FAILURE);
                }
            }
        }
        if (readType.empty()) {
            PBLOG_BLOCK_FATAL(
                "Input checker",
                "Could not determine read type, read groups are missing : " + inputFile);
            std::exit(EXIT_FAILURE);
        }
        if (readType == "CCS") {
            if (!settings.InputCCSFile.empty() && !settings.ccsQuery) {
                PBLOG_BLOCK_FATAL("Input checker", "Multiple CCS files detected!");
                PBLOG_BLOCK_FATAL("Input checker", "1) " + settings.InputCCSFile);
                PBLOG_BLOCK_FATAL("Input checker", "2) " + inputFile);
                std::exit(EXIT_FAILURE);
            }
            if (settings.InputCCSFile.empty()) {
                settings.InputCCSFile = inputFile;
            } else {
                settings.InputCLRFile = inputFile;
            }
        } else if (readType == "SUBREAD") {
            if (!settings.InputCLRFile.empty()) {
                PBLOG_BLOCK_FATAL("Input checker", "Multiple CLR files detected!");
                PBLOG_BLOCK_FATAL("Input checker", "1) " + settings.InputCLRFile);
                PBLOG_BLOCK_FATAL("Input checker", "2) " + inputFile);
                std::exit(EXIT_FAILURE);
            }
            settings.InputCLRFile = inputFile;
        } else {
            PBLOG_BLOCK_FATAL("Input checker", "Unknown read type in : " + inputFile);
            std::exit(EXIT_FAILURE);
        }
    };

    ReadType(files[0]);
    ReadType(files[1]);
    settings.OutputAlignmentFile = files[2];
    settings.NumThreads = options.NumThreads();

    IO::BamZmwReaderConfig zmwReaderConfig{options};

    const auto clrFiles = BAM::DataSet(settings.InputCLRFile).BamFiles();
    if (clrFiles.size() != 1) {
        PBLOG_BLOCK_FATAL("Input checker", "CLR input must be exactly one BAM file! Found " +
                                               std::to_string(clrFiles.size()));
        std::exit(EXIT_FAILURE);
    }

    if (!clrFiles[0].PacBioIndexExists()) {
        PBLOG_BLOCK_FATAL("Input checker", "Missing PBI file for " + clrFiles[0].Filename());
        PBLOG_BLOCK_FATAL("Input checker",
                          "Please generate one with : pbindex " + clrFiles[0].Filename());
        PBLOG_BLOCK_FATAL("Input checker",
                          "You can get pbindex from bioconda: conda install -c bioconda pbbam");
        std::exit(EXIT_FAILURE);
    }

    std::unordered_map<int32_t, int64_t> holenumberToOffset;
    {
        BAM::PbiRawData pbi{clrFiles[0].PacBioIndexFilename()};
        const int32_t numPbiRecord = pbi.BasicData().holeNumber_.size();
        int32_t curHoleNumber = pbi.BasicData().holeNumber_[0];
        holenumberToOffset.insert({curHoleNumber, pbi.BasicData().fileOffset_[0]});
        for (int32_t i = 1; i < numPbiRecord; ++i) {
            if (pbi.BasicData().holeNumber_[i] != curHoleNumber) {
                curHoleNumber = pbi.BasicData().holeNumber_[i];
                holenumberToOffset.insert({curHoleNumber, pbi.BasicData().fileOffset_[i]});
            }
        }
    }

    BAM::BamReader clrFile{clrFiles[0].Filename()};
    {
        auto ccsFiles = BAM::DataSet(settings.InputCCSFile).BamFiles();
        if (ccsFiles.size() != 1) {
            PBLOG_BLOCK_FATAL("Input checker", "Expecting exactly one CCS BAM file, found " +
                                                   std::to_string(ccsFiles.size()));
            std::exit(EXIT_FAILURE);
        }
    }

    BAM::BamRecord clrRecord;
    clrFile.GetNext(clrRecord);

    BAM::BamHeader header = clrFile.Header().DeepCopy();
    int32_t numCcsReads = 0;
    {
        std::string outputFastaName = settings.OutputAlignmentFile;
        boost::replace_all(outputFastaName, ".bam", ".fasta");
        BAM::FastaWriter fasta{outputFastaName};
        IO::BamZmwReader ccsReader{settings.InputCCSFile, zmwReaderConfig};

        IO::ZmwRecords zmwRecords;
        PBLOG_BLOCK_INFO("Fasta CCS", "Start writing CCS reads to " + outputFastaName);
        while (ccsReader.GetNext(zmwRecords)) {
            if ((numCcsReads % 10000) == 0) {
                PBLOG_BLOCK_INFO("Fasta CCS", std::to_string(numCcsReads));
            }
            if (zmwRecords.InputRecords.empty()) {
                PBLOG_BLOCK_FATAL("CCS reader", "CCS ZMW " + std::to_string(zmwRecords.HoleNumber) +
                                                    " has no records!");
                continue;
            }

            if (std::ssize(zmwRecords.InputRecords) != 1) {
                PBLOG_BLOCK_FATAL("CCS reader", "CCS ZMW " + std::to_string(zmwRecords.HoleNumber) +
                                                    " has multiple records. Ignoring ZMW!");
                continue;
            }

            const auto ccsRecord = zmwRecords.InputRecords[0];
            const std::string seq = ccsRecord.Sequence();
            const std::string name = ccsRecord.FullName();
            header.AddSequence({name, std::to_string(seq.size())});
            fasta.Write(name, seq);
            ++numCcsReads;
        }
        PBLOG_BLOCK_INFO("Fasta CCS", std::to_string(numCcsReads));
    }

    BAM::ProgramInfo program("actc");
    program.Name("actc")
        .CommandLine(options.InputCommandLine())
        .Version(Actc::LibraryInfo().Release);
    header.AddProgram(program);

    BAM::BamWriter writer(settings.OutputAlignmentFile, header, BAM::BamWriter::DefaultCompression,
                          settings.NumThreads);

    Parallel::WorkQueue<std::vector<BAM::BamRecord>> workQueue(settings.NumThreads, 10);
    std::future<void> workerThread = std::async(std::launch::async, WorkerThread,
                                                std::ref(workQueue), std::ref(writer), numCcsReads);

    const auto Submit = [&header](const std::vector<BAM::BamRecord>& clrRecords,
                                  const BAM::BamRecord& ccsRecord, const int32_t curCcsIdx,
                                  const bool ccs) {
        const std::vector<AlnResults> alns =
            PancakeAlignerSubread(clrRecords, ccsRecord.Sequence());
        std::vector<BAM::BamRecord> alnRecords;
        int32_t subreadIdx = 0;
        for (const auto& aln : alns) {
            for (const auto& a : aln) {
                if (a->isAligned) {
                    alnRecords.emplace_back(
                        AlnToBam(curCcsIdx, header, *a, clrRecords[subreadIdx], ccs));
                }
            }
            ++subreadIdx;
        }
        return alnRecords;
    };

    int32_t curCcsIdx = 0;
    IO::BamZmwReader ccsReader{settings.InputCCSFile, zmwReaderConfig};
    IO::ZmwRecords zmwRecords;
    while (ccsReader.GetNext(zmwRecords)) {
        if (zmwRecords.InputRecords.empty()) {
            PBLOG_BLOCK_FATAL("CCS reader", "CCS ZMW " + std::to_string(zmwRecords.HoleNumber) +
                                                " has no records!");
            continue;
        }

        if (std::ssize(zmwRecords.InputRecords) != 1) {
            PBLOG_BLOCK_FATAL("CCS reader", "CCS ZMW " + std::to_string(zmwRecords.HoleNumber) +
                                                " has multiple records. Ignoring ZMW!");
            continue;
        }

        const auto ccsRecord = zmwRecords.InputRecords[0];

        PBLOG_BLOCK_DEBUG("CCS reader", ccsRecord.FullName());
        const int32_t holeNumber = ccsRecord.HoleNumber();
        if (clrRecord.HoleNumber() != holeNumber) {
            if (holenumberToOffset.find(holeNumber) == holenumberToOffset.cend()) {
                if (!settings.ccsQuery) {
                    PBLOG_BLOCK_FATAL("CLR reader", "ZMW " + std::to_string(holeNumber) +
                                                        " missing in CLR file )" +
                                                        clrFile.Filename());
                    std::exit(EXIT_FAILURE);
                } else {
                    PBLOG_BLOCK_WARN("CLR reader", "ZMW " + std::to_string(holeNumber) +
                                                       " missing in second file )" +
                                                       clrFile.Filename());
                    ++curCcsIdx;
                    continue;
                }
            }
            PBLOG_BLOCK_DEBUG("CLR parser", "SEEKING");
            clrFile.VirtualSeek(holenumberToOffset[holeNumber]);
            clrFile.GetNext(clrRecord);
        }
        std::vector<BAM::BamRecord> clrRecords;
        do {
            PBLOG_BLOCK_DEBUG("CLR parser", clrRecord.FullName());
            clrRecords.emplace_back(clrRecord);
            if (!clrFile.GetNext(clrRecord)) {
                break;
            }
        } while (clrRecord.HoleNumber() == holeNumber);

        workQueue.ProduceWith(Submit, std::move(clrRecords), ccsRecord, curCcsIdx,
                              settings.ccsQuery);

        ++curCcsIdx;
    }

    workQueue.FinalizeWorkers();
    workerThread.wait();
    workQueue.Finalize();

    globalTimer.Freeze();
    PBLOG_BLOCK_INFO("Run Time", globalTimer.ElapsedTime());
    PBLOG_BLOCK_INFO("CPU Time", Utility::Stopwatch::PrettyPrintNanoseconds(static_cast<int64_t>(
                                     Utility::Stopwatch::CpuTime() * 1000 * 1000 * 1000)));

    int64_t const peakRss = PacBio::Utility::MemoryConsumption::PeakRss();
    double const peakRssGb = peakRss / 1024.0 / 1024.0 / 1024.0;
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(3) << peakRssGb << " GB";
    PBLOG_BLOCK_INFO("Peak RSS", ss.str());

    return EXIT_SUCCESS;
}
}  // namespace PacBio

int main(int argc, char* argv[])
{
    return PacBio::CLI_v2::Run(argc, argv, PacBio::CreateCLI(), &PacBio::RunnerSubroutine);
}
