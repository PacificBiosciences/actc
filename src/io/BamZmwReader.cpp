#include "BamZmwReader.hpp"

#include "ZmwRecords.hpp"

#include <pbbam/BamReader.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/Alarm.h>

#include <boost/algorithm/string/predicate.hpp>

#include <string>
#include <utility>

namespace PacBio {
namespace IO {
namespace {

struct UniqueZmw
{
    std::int32_t PbiIdx;
    std::int32_t HoleNumber;
    std::int64_t FileOffset;
};

// Get the unique ZMWs from the input dataset, applying filters
std::vector<UniqueZmw> UniqueZmws(const BAM::DataSet& ds, const bool dieOnError)
{
    const std::vector<BAM::BamFile> bamFiles = ds.BamFiles();
    // Only works with ONE input BAM file
    if (std::ssize(bamFiles) != 1) {
        if (dieOnError) {
            throw PB_CLI_ALARM("Chunking only works with one input BAM file!");
        } else {
            return {};
        }
    }

    // Input BAM file must have a PBI
    if (!bamFiles[0].PacBioIndexExists()) {
        if (dieOnError) {
            throw PB_CLI_ALARM(
                "PBI file is missing for input BAM file! Please create one using pbindex!");
        } else {
            return {};
        }
    }

    // Various filters
    bool fromZMWSet = false;
    bool toZMWSet = false;
    std::int32_t fromZMW{};
    std::int32_t toZMW{};
    bool moduloSet = false;
    std::uint32_t modulus = 0;
    BAM::PbiFilter zmwFilter{};

    using Operator = BAM::Compare::Type;

    // Read all available filters
    for (const BAM::Filter& filter : ds.Filters()) {
        for (const BAM::Property& property : filter.Properties()) {
            if (property.Name() != "zm") {
                continue;
            }

            const Operator op = BAM::Compare::TypeFromOperator(property.Operator());
            const std::int32_t value = [&]() {
                const std::string& valueString = property.Value();
                for (const char c : valueString) {
                    if (!std::isdigit(c)) {
                        return -1;
                    }
                }
                return std::stoi(valueString);
            }();

            switch (op) {
                case Operator::LESS_THAN_EQUAL:
                    toZMW = value + 1;
                    toZMWSet = true;
                    zmwFilter.Add(BAM::PbiZmwFilter{value, Operator::LESS_THAN_EQUAL});
                    break;
                case Operator::LESS_THAN:
                    toZMW = value;
                    toZMWSet = true;
                    zmwFilter.Add(BAM::PbiZmwFilter{value, Operator::LESS_THAN});
                    break;
                case Operator::GREATER_THAN:
                    fromZMW = value;
                    fromZMWSet = true;
                    zmwFilter.Add(BAM::PbiZmwFilter{value, Operator::GREATER_THAN});
                    break;
                case Operator::GREATER_THAN_EQUAL:
                    fromZMW = value - 1;
                    fromZMWSet = true;
                    zmwFilter.Add(BAM::PbiZmwFilter{value, Operator::GREATER_THAN_EQUAL});
                    break;
                case Operator::EQUAL: {
                    const auto& attributes = property.Attributes();
                    if (attributes.find("Modulo") != attributes.cend()) {
                        const std::string& hash = property.Attribute("Hash");
                        const BAM::FilterHash hashType = [&]() {
                            if (boost::iequals(hash, "uint32cast")) {
                                return BAM::FilterHash::UNSIGNED_LONG_CAST;
                            }
                            if (boost::iequals(hash, "boosthashcombine")) {
                                return BAM::FilterHash::BOOST_HASH_COMBINE;
                            }
                            throw PB_CLI_ALARM("Unsupported hash type: " + hash);
                        }();

                        modulus = std::stoul(property.Attribute("Modulo"));
                        const std::uint32_t modulusValue = std::stoul(property.Value());
                        zmwFilter.Add(BAM::PbiZmwModuloFilter{modulus, modulusValue, hashType});
                        moduloSet = true;
                    }
                    break;
                }
                default:
                    throw PB_CLI_ALARM(
                        "Unsupported operator type for ZMW range filter. "
                        "Supported are: <=, <, >, >=");
            }
        }
    }

    // ZMW filter range
    if (fromZMWSet && toZMWSet) {
        PBLOG_BLOCK_INFO("ZMW filter range",
                         '(' + std::to_string(fromZMW) + ',' + std::to_string(toZMW) + ')');
    }

    // ZMW downsample filter
    if (moduloSet) {
        const double factor = (static_cast<double>(100) / modulus);
        PBLOG_BLOCK_INFO("ZMW downsample", std::to_string(factor) + '%');
    }

    const BAM::PbiRawData index{bamFiles[0].PacBioIndexFilename()};
    const std::vector<std::int32_t>& zmws = index.BasicData().holeNumber_;
    const std::vector<std::int64_t>& fileOffset = index.BasicData().fileOffset_;
    const std::int32_t numRecords = std::ssize(zmws);
    if (numRecords == 0) {
        throw PB_ALARM("InputDataError", "No input records in PBI file!");
    }

    // Get the first record per ZMW
    std::vector<UniqueZmw> result;
    result.reserve(numRecords);
    if (zmwFilter.Accepts(index, 0)) {
        result.emplace_back(UniqueZmw{0, zmws[0], fileOffset[0]});
    }
    for (std::int32_t i = 1; i < numRecords; ++i) {
        if (zmws[i] != zmws[i - 1] && zmwFilter.Accepts(index, i)) {
            result.emplace_back(UniqueZmw{i, zmws[i], fileOffset[i]});
        }
    }

#ifndef NDEBUG
    PBLOG_DEBUG << "I INDEX ZMW OFFSET START END";
    for (std::int32_t i = 0; i < std::ssize(result); ++i) {
        PBLOG_DEBUG << i << ' ' << result[i].PbiIdx << ' ' << result[i].HoleNumber << ' '
                    << result[i].FileOffset << ' ' << index.BasicData().qStart_[result[i].PbiIdx]
                    << ' ' << index.BasicData().qEnd_[result[i].PbiIdx];
    }
#endif

    return result;
}

// Create BAM reader for a given file path (might have filters) and configuration (might chunk)
void CreateBamReader(const std::filesystem::path& filePath, const BamZmwReaderConfig& config,
                     std::int32_t& numZmws, std::unique_ptr<BAM::internal::IQuery>& query,
                     std::int32_t& endZmwHoleNumber)
{
    const BAM::DataSet dataset{filePath};
    const BAM::PbiFilter filter = BAM::PbiFilter::FromDataSet(dataset);

    // Indicate that there should be no end ZMW stopping by default / last chunk
    endZmwHoleNumber = -1;

    // Local variable to ease access to the chunking parameters
    const std::int32_t chunkNumerator = config.ChunkNumerator;
    const std::int32_t chunkDenominator = config.ChunkDenominator;

    // Is the first chunk requested?
    const bool firstChunk{chunkNumerator == 1};
    // Is the last chunk requested?
    const bool lastChunk{(chunkNumerator == chunkDenominator) && (chunkDenominator > 0)};
    // Access to the first ZMW of the chunk
    UniqueZmw startZmw{-1, -1, -1};
    // Is a specific chunk requested?
    if ((chunkNumerator > 0) && (chunkDenominator > 0)) {
        // Can we use ZMW chunking? Not with filters!
        if (!filter.IsEmpty()) {
            throw PB_CLI_ALARM(
                "Cannot combine dataset filters with --chunk. Please use ZMW chunking via "
                "dataset filters or remove filters from dataset.");
        }

        // Get the unique ZMWs
        const std::vector<UniqueZmw> zmwsUniq = UniqueZmws(dataset, true);
        const std::int32_t numZmwsAll = std::ssize(zmwsUniq);
        // Check if we have enough ZMWs
        if (numZmwsAll < chunkDenominator) {
            throw PB_CLI_ALARM(
                "Fewer ZMWs available than specified chunks: " + std::to_string(numZmwsAll) +
                " vs. " + std::to_string(chunkDenominator));
        }
        // Calculate the number of ZMWs per chunk
        const double chunkSize = 1.0 * numZmwsAll / chunkDenominator;
        // Calculate the start ZMW of the chunk
        const std::int32_t firstChunkIdx =
            firstChunk ? 0 : std::round(chunkSize * (chunkNumerator - 1));
        // Store the start ZMW
        startZmw = zmwsUniq[firstChunkIdx];
        // Calculate the end ZMW of the chunk
        const std::int32_t lastChunkIdx =
            lastChunk ? std::ssize(zmwsUniq) - 1 : std::round(chunkSize * (chunkNumerator));
        // Update the number of ZMWs in the chunk
        numZmws = lastChunkIdx - firstChunkIdx + static_cast<std::int32_t>(lastChunk);
        PBLOG_BLOCK_INFO("Chunk index",
                         std::to_string(chunkNumerator) + '/' + std::to_string(chunkDenominator));
        PBLOG_BLOCK_INFO("ZMW range", '[' + std::to_string(startZmw.HoleNumber) + ',' +
                                          std::to_string(zmwsUniq[lastChunkIdx].HoleNumber) +
                                          (lastChunk ? ']' : ')'));
        // Store the end ZMW, which is exclusive in the range
        if (!lastChunk) {
            endZmwHoleNumber = zmwsUniq[lastChunkIdx].HoleNumber;
        }
    }

    if (filter.IsEmpty()) {  // No filter used, we allow to chunk
        // Do not allow more than ONE BAM file
        if (std::ssize(dataset.BamFiles()) != 1) {
            throw PB_CLI_ALARM("Input must have exactly one BAM file.");
        }
        // Create a BAM reader for the first BAM file
        query = std::make_unique<BAM::BamReader>(dataset.BamFiles().front());
        auto& reader = dynamic_cast<BAM::BamReader&>(*query);
        if (startZmw.FileOffset != -1) {
            PBLOG_BLOCK_DEBUG("Chunking", "File offset " + std::to_string(startZmw.FileOffset));
            reader.VirtualSeek(startZmw.FileOffset);
        }
    } else {  // Filter used, we do not allow to chunk
        query = std::make_unique<BAM::PbiFilterQuery>(filter, dataset);
    }
}

}  // namespace

BamZmwReader::BamZmwReader(std::filesystem::path path, BamZmwReaderConfig config)
    : path_{std::move(path)}, config_{std::move(config)}
{
    CreateBamReader(path_, config_, numZmws_, reader_, endZmwHoleNumber_);
}

bool BamZmwReader::GetNext(ZmwRecords& zmw)
{
    if (endOfFile_) {
        return false;
    }

    // This is the first record of the BAM file
    if (!lastRecord_) {
        BAM::BamRecord record;
        if (!reader_->GetNext(record)) {
            PBLOG_BLOCK_WARN("BamZmwReader", "Input BAM is empty");
            return false;
        }
        PBLOG_BLOCK_TRACE("BamZmwReader", "BamRecord first " + record.FullName());
        lastRecord_ = std::move(record);
    }

    assert(lastRecord_);

    // Check if this is past the last ZMW of the chunk
    std::int32_t holeNumber = lastRecord_->HoleNumber();
    if ((endZmwHoleNumber_ != -1) && (endZmwHoleNumber_ == holeNumber)) {
        return false;
    }

    // Start a new vector of records for this ZMW
    std::vector<BAM::BamRecord> inputBamRecords;
    inputBamRecords.emplace_back(std::move(*lastRecord_));

    const auto changeOwnership = [&zmw, &inputBamRecords]() {
        zmw = ZmwRecords{};
        if (!inputBamRecords.empty()) {
            zmw.HoleNumber = inputBamRecords.front().HoleNumber();
            zmw.MovieName = inputBamRecords.front().MovieName();
            zmw.InputRecords = std::move(inputBamRecords);
        } else {
            PBLOG_BLOCK_TRACE("BamZmwReader", "Empty ZMW");
        }
    };

    // Read next records until we have a full ZMW
    for (const auto& read : *reader_) {
        PBLOG_BLOCK_TRACE("BamZmwReader", "BamRecord reading " + read.FullName());
        // Check if we've started a new ZMW
        if (holeNumber != read.HoleNumber()) {
            // Stop adding records to this ZMW and change ownership
            changeOwnership();
            // Start a new ZMW with this record
            lastRecord_ = read;
            // Return true to indicate that we have a new ZMW
            return true;
        }
        inputBamRecords.emplace_back(read);
    }

    // If we have reached this, we have reached the end of the file
    changeOwnership();
    PBLOG_BLOCK_TRACE("BamZmwReader",
                      "Last ZMW " + std::to_string(zmw.HoleNumber) + ' ' + zmw.MovieName +
                          " with " + std::to_string(std::ssize(zmw.InputRecords)) + " records");
    // For sanity
    lastRecord_ = std::nullopt;
    // Indicate that we have reached the end of the file
    endOfFile_ = true;
    return true;
}

}  // namespace IO
}  // namespace PacBio
