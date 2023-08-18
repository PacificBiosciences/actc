#include "BamZmwReaderConfig.hpp"

#include <pbcopper/cli2/internal/BuiltinOptions.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/Alarm.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace PacBio {
namespace IO {
namespace {

// Determine start and end hole numbers for a given chunk
std::pair<std::int32_t, std::int32_t> DetermineChunk(const std::string& chunk)
{
    const auto ChunkErrorAndDie = [&]() {
        throw PB_CLI_ALARM(
            "Wrong format for --chunk, please provide two integers separated by a slash like 2/10. "
            "First number must be less than the second number. Both must be positive and greater "
            "than 0.");
    };

    std::int32_t chunkNumerator{-1};
    std::int32_t chunkDenominator{-1};

    // Empty chunk means we want to process all holes and return early
    if (chunk.empty()) {
        return std::make_pair(chunkNumerator, chunkDenominator);
    }

    // Check if chunk is in the format i/N
    if (!boost::find_first(chunk, "/")) {
        ChunkErrorAndDie();
    }

    // Split chunk into numerator and denominator strings
    std::vector<std::string> splits;
    boost::split(splits, chunk, boost::is_any_of("/"));
    // Check if we have exactly two strings
    if (splits.size() != 2) {
        ChunkErrorAndDie();
    }
    try {
        // Try to convert numerator and denominator to integers
        chunkNumerator = boost::lexical_cast<float>(splits[0]);
        chunkDenominator = boost::lexical_cast<float>(splits[1]);

        // Check if numerator is less than denominator and both are greater than 0
        if (chunkNumerator > chunkDenominator || chunkNumerator <= 0 || chunkDenominator <= 0) {
            ChunkErrorAndDie();
        }
    } catch (const boost::bad_lexical_cast&) {
        // Conversion failed
        ChunkErrorAndDie();
    }

    return std::make_pair(chunkNumerator, chunkDenominator);
}

}  // namespace

BamZmwReaderConfig::BamZmwReaderConfig(const CLI_v2::Results& options)
{
    // Parse chunk string
    std::tie(ChunkNumerator, ChunkDenominator) = DetermineChunk(options[OptionNames::Chunk]);
}

}  // namespace IO
}  // namespace PacBio
