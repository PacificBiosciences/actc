#include "AlignmentResult.hpp"
#include "PancakeAligner.hpp"

#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>
#include <pbbam/DataSet.h>
#include <pbbam/FastaWriter.h>
#include <pbbam/PbiRawData.h>
#include <pbcopper/cli2/CLI.h>
#include <pbcopper/cli2/internal/BuiltinOptions.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/parallel/WorkQueue.h>
#include <pbcopper/utility/MemoryConsumption.h>
#include <pbcopper/utility/Ssize.h>
#include <pbcopper/utility/Stopwatch.h>
#include <boost/algorithm/string.hpp>

#include <cstdint>
#include <cstdlib>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>

namespace PacBio {

struct DaiDaiSettings
{
    std::string InputCLRFile;
    std::string InputCCSFile;
    std::string OutputAlignmentFile;
    int32_t NumThreads = 1;
};

CLI_v2::Interface CreateCLI()
{
    static const std::string description{"Align clr to ccs reads."};
    const auto version = "0.0.1";
    CLI_v2::Interface i{"actc", description, version};

    Logging::LogConfig logConfig;
    logConfig.Header = "| ";
    logConfig.Delimiter = " | ";
    logConfig.Fields = Logging::LogField::TIMESTAMP | Logging::LogField::LOG_LEVEL;
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

    return i;
}

void WorkerThread(Parallel::WorkQueue<std::vector<BAM::BamRecord>>& queue, BAM::BamWriter& writer,
                  const int32_t numReads)
{
    int32_t counter = 0;
    double perc = 0.1;

    auto LambdaWorker = [&](std::vector<BAM::BamRecord>&& ps) {
        ++counter;
        if (1.0 * counter / numReads > perc) {
            PBLOG_INFO << "Progress " << 100 * perc << '%';
            perc = std::min(perc + 0.1, 1.0);
        }
        for (const auto& record : ps) {
            writer.Write(record);
        }
    };

    while (queue.ConsumeWith(LambdaWorker)) {
    }
}

int RunnerSubroutine(const CLI_v2::Results& options)
{
    Utility::Stopwatch globalTimer;

    DaiDaiSettings settings;
    const std::vector<std::string> files = options.PositionalArguments();

    const auto ReadType = [&settings](const std::string& inputFile) {
        const auto bamFiles = BAM::DataSet(inputFile).BamFiles();
        if (bamFiles.empty()) {
            PBLOG_FATAL << "No BAM files available for: " << inputFile;
            std::exit(EXIT_FAILURE);
        }

        std::string readType;
        for (int32_t i = 0; i < Utility::Ssize(bamFiles); ++i) {
            for (const auto& rg : bamFiles.at(i).Header().ReadGroups()) {
                if (readType.empty()) {
                    readType = rg.ReadType();
                } else if (readType != rg.ReadType()) {
                    PBLOG_FATAL << "Do not mix and match different read types for input file : "
                                << inputFile;
                    std::exit(EXIT_FAILURE);
                }
            }
        }
        if (readType.empty()) {
            PBLOG_FATAL << "Could not determine read type, read groups are missing : " << inputFile;
            std::exit(EXIT_FAILURE);
        }
        if (readType == "CCS") {
            if (!settings.InputCCSFile.empty()) {
                PBLOG_FATAL << "Multiple CCS files detected!";
                PBLOG_FATAL << "1) " << settings.InputCCSFile;
                PBLOG_FATAL << "2) " << inputFile;
                std::exit(EXIT_FAILURE);
            }
            settings.InputCCSFile = inputFile;
        } else if (readType == "SUBREAD") {
            if (!settings.InputCLRFile.empty()) {
                PBLOG_FATAL << "Multiple CLR files detected!";
                PBLOG_FATAL << "1) " << settings.InputCLRFile;
                PBLOG_FATAL << "2) " << inputFile;
                std::exit(EXIT_FAILURE);
            }
            settings.InputCLRFile = inputFile;
        } else {
            PBLOG_FATAL << "Unknown read type in : " << inputFile;
            std::exit(EXIT_FAILURE);
        }
    };

    ReadType(files[0]);
    ReadType(files[1]);
    settings.OutputAlignmentFile = files[2];
    settings.NumThreads = options.NumThreads();

    const auto clrFiles = BAM::DataSet(settings.InputCLRFile).BamFiles();
    if (clrFiles.size() != 1) {
        PBLOG_FATAL << "CLR input must be exactly one BAM file! Found " << clrFiles.size();
        std::exit(EXIT_FAILURE);
    }

    if (!clrFiles[0].PacBioIndexExists()) {
        PBLOG_FATAL << "Missing PBI file for " << clrFiles[0].Filename();
        PBLOG_FATAL << "Please generate one with : pbindex " << clrFiles[0].Filename();
        PBLOG_FATAL << "You can get pbindex from bioconda: conda install -c bioconda pbbam";
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
    auto ccsFiles = BAM::DataSet(settings.InputCCSFile).BamFiles();
    if (ccsFiles.size() != 1) {
        PBLOG_FATAL << "Expecting exactly one CCS BAM file, found " << ccsFiles.size();
        std::exit(EXIT_FAILURE);
    }

    BAM::BamRecord clrRecord;
    clrFile.GetNext(clrRecord);
    const int32_t firstHoleNumber = clrRecord.HoleNumber();
    const auto GetNextRecord = [&]() {
        if (!clrFile.GetNext(clrRecord)) {
            PBLOG_VERBOSE << "End of CLR file, rewind";
            clrFile.VirtualSeek(holenumberToOffset[firstHoleNumber]);
            if (!clrFile.GetNext(clrRecord)) {
                PBLOG_FATAL << "I'm confused, no data in CLR file";
                std::exit(EXIT_FAILURE);
            }
        }
    };

    const std::string ccsFileName = ccsFiles[0].Filename();

    BAM::BamReader reader{ccsFileName};
    BAM::BamHeader header = reader.Header().DeepCopy();
    int32_t numCcsReads = 0;
    {
        std::string outputFastaName = settings.OutputAlignmentFile;
        boost::replace_all(outputFastaName, ".bam", ".fasta");
        BAM::FastaWriter fasta{outputFastaName};
        for (const auto& ccsRecord : BAM::BamReader{ccsFileName}) {
            const std::string seq = ccsRecord.Sequence();
            const std::string name = ccsRecord.FullName();
            header.AddSequence({name, std::to_string(seq.size())});
            fasta.Write(name, seq);
            ++numCcsReads;
        }
    }

    BAM::BamWriter writer(settings.OutputAlignmentFile, header);

    Parallel::WorkQueue<std::vector<BAM::BamRecord>> workQueue(settings.NumThreads, 3);
    std::future<void> workerThread = std::async(std::launch::async, WorkerThread,
                                                std::ref(workQueue), std::ref(writer), numCcsReads);

    const auto Submit = [&header](const std::vector<BAM::BamRecord>& clrRecords,
                                  const BAM::BamRecord& ccsRecord, const int32_t curCcsIdx) {
        const std::vector<AlnResults> alns =
            PancakeAlignerSubread(clrRecords, ccsRecord.Sequence());
        std::vector<BAM::BamRecord> alnRecords;
        int32_t subreadIdx = 0;
        for (const auto& aln : alns) {
            for (const auto& a : aln) {
                if (a->isAligned) {
                    alnRecords.emplace_back(
                        AlnToBam(curCcsIdx, header, *a, clrRecords[subreadIdx]));
                }
            }
            ++subreadIdx;
        }
        return alnRecords;
    };

    int32_t curCcsIdx = 0;
    for (const auto& ccsRecord : reader) {
        PBLOG_DEBUG << "CCS : " << ccsRecord.FullName();
        const int32_t holeNumber = ccsRecord.HoleNumber();
        if (clrRecord.HoleNumber() != holeNumber) {
            if (holenumberToOffset.find(holeNumber) == holenumberToOffset.cend()) {
                PBLOG_FATAL << "ZMW " << holeNumber << " missing in CLR file "
                            << clrFile.Filename();
                std::exit(EXIT_FAILURE);
            }
            PBLOG_DEBUG << "SEEKING";
            clrFile.VirtualSeek(holenumberToOffset[holeNumber]);
            GetNextRecord();
        }
        std::vector<BAM::BamRecord> clrRecords;
        do {
            PBLOG_DEBUG << "CLR : " << clrRecord.FullName();
            clrRecords.emplace_back(clrRecord);
            GetNextRecord();
        } while (clrRecord.HoleNumber() == holeNumber);

        workQueue.ProduceWith(Submit, std::move(clrRecords), ccsRecord, curCcsIdx);

        ++curCcsIdx;
    }

    workQueue.FinalizeWorkers();
    workerThread.wait();
    workQueue.Finalize();

    globalTimer.Freeze();
    PBLOG_INFO << "Run Time : " << globalTimer.ElapsedTime();
    PBLOG_INFO << "CPU Time : "
               << Utility::Stopwatch::PrettyPrintNanoseconds(
                      static_cast<int64_t>(Utility::Stopwatch::CpuTime() * 1000 * 1000 * 1000));

    int64_t const peakRss = PacBio::Utility::MemoryConsumption::PeakRss();
    double const peakRssGb = peakRss / 1024.0 / 1024.0 / 1024.0;
    PBLOG_INFO << "Peak RSS : " << std::fixed << std::setprecision(3) << peakRssGb << " GB";

    return EXIT_SUCCESS;
}
}  // namespace PacBio

int main(int argc, char* argv[])
{
    return PacBio::CLI_v2::Run(argc, argv, PacBio::CreateCLI(), &PacBio::RunnerSubroutine);
}
