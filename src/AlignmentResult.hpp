#pragma once

#include <pbbam/BamHeader.h>
#include <pbbam/BamRecord.h>
#include <pbcopper/data/Cigar.h>

#include <cstdint>
#include <iosfwd>
#include <memory>

namespace PacBio {

class AlignmentResult
{
public:
    AlignmentResult() = default;

    AlignmentResult(int32_t rId, bool rReversed, int64_t rStart, int64_t rEnd, int64_t qStart,
                    int64_t qEnd, int64_t qLen, Data::Cigar cigar, uint8_t mapq, int32_t as,
                    bool isAligned, bool isSupplementary, bool isSecondary);

    AlignmentResult(bool rReversed, int64_t rStart, int64_t rEnd, int64_t qStart, int64_t qEnd,
                    int64_t qLen, Data::Cigar cigar);

    std::unique_ptr<AlignmentResult> Clip(int64_t clipQStart, int64_t clipQEnd, int64_t clipRStart,
                                          int64_t clipREnd) const;

    bool operator==(const AlignmentResult& op) const noexcept;

    int32_t rId = 0;
    bool rReversed = false;
    int64_t rStart = 0;
    int64_t rEnd = 0;
    int64_t qStart = 0;
    int64_t qEnd = 0;
    int64_t qLen = 0;
    Data::Cigar cigar;
    uint8_t mapq = 0;
    int32_t as = 0;
    bool isAligned = false;
    bool isSupplementary = false;
    bool isSecondary = false;
};

std::ostream& operator<<(std::ostream& out, const AlignmentResult& a);

using AlnResults = std::vector<std::unique_ptr<AlignmentResult>>;

BAM::BamRecord AlnToBam(const int32_t refId, const BAM::BamHeader& header,
                        const AlignmentResult& aln, const BAM::BamRecord& read, bool ccs);

}  // namespace PacBio
