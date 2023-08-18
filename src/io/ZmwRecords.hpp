#ifndef Actc_ZMW_HPP
#define Actc_ZMW_HPP

#include <pbbam/BamRecord.h>
#include <pbcopper/data/Strand.h>

#include <optional>
#include <string>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace IO {

// A ZMW
struct ZmwRecords
{
    // Immutable data
    std::string MovieName;
    std::size_t HoleNumber;
    std::vector<BAM::BamRecord> InputRecords;
};

}  // namespace IO
}  // namespace PacBio

#endif  // Actc_ZMW_HPP
