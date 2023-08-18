#ifndef Actc_IO_BAMZMWREADERCONFIG_HPP
#define Actc_IO_BAMZMWREADERCONFIG_HPP

#include <pbcopper/cli2/Option.h>
#include <pbcopper/cli2/Results.h>

#include <vector>

#include <cstdint>

namespace PacBio {
namespace IO {
namespace OptionNames {

// clang-format off
const CLI_v2::Option Chunk{
R"({
    "names" : ["chunk"],
    "description" : "Operate on a single chunk. Format i/N, where i in [1,N]. Examples: 3/24 or 9/9",
    "type" : "string",
    "default" : ""
})"
};
// clang-format on
}  // namespace OptionNames

struct BamZmwReaderConfig
{
    BamZmwReaderConfig(const CLI_v2::Results& options);

    static std::vector<CLI_v2::Option> GenerateOptionGroup();

    std::int32_t ChunkNumerator{-1};
    std::int32_t ChunkDenominator{-1};
};

}  // namespace IO
}  // namespace PacBio

#endif  // Actc_IO_BAMZMWREADERCONFIG_HPP
