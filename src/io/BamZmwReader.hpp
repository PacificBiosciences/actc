#ifndef Actc_IO_BAMZMWREADER_HPP
#define Actc_IO_BAMZMWREADER_HPP

#include "BamZmwReaderConfig.hpp"
#include "ZmwRecords.hpp"

#include <pbbam/BamRecord.h>
#include <pbbam/internal/QueryBase.h>

#include <cstdint>

#include <filesystem>
#include <memory>
#include <optional>

namespace PacBio {
namespace IO {

class BamZmwReader : public BAM::internal::QueryBase<ZmwRecords>
{
public:
    explicit BamZmwReader(std::filesystem::path path, BamZmwReaderConfig config);
    bool GetNext(ZmwRecords& zmw) final;

private:
    std::filesystem::path path_;
    BamZmwReaderConfig config_;
    std::unique_ptr<BAM::internal::IQuery> reader_;
    std::int32_t numZmws_;
    std::int32_t endZmwHoleNumber_;
    std::optional<BAM::BamRecord> lastRecord_;
    bool endOfFile_{false};
};

}  // namespace IO
}  // namespace PacBio

#endif  // Actc_IO_BAMZMWREADER_HPP
