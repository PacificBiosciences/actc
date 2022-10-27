#include "AlignmentResult.hpp"

#include <pbcopper/utility/SequenceUtils.h>
#include "AlignerUtils.hpp"

#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <tuple>

namespace PacBio {

AlignmentResult::AlignmentResult(int32_t rIdArg, bool rReversedArg, int64_t rStartArg,
                                 int64_t rEndArg, int64_t qStartArg, int64_t qEndArg,
                                 int64_t qLenArg, Data::Cigar cigarArg, uint8_t mapqArg,
                                 int32_t asArg, bool isAlignedArg, bool isSupplementaryArg,
                                 bool isSecondaryArg)
    : rId{rIdArg}
    , rReversed{rReversedArg}
    , rStart{rStartArg}
    , rEnd{rEndArg}
    , qStart{qStartArg}
    , qEnd{qEndArg}
    , qLen{qLenArg}
    , cigar{std::move(cigarArg)}
    , mapq{mapqArg}
    , as{asArg}
    , isAligned{isAlignedArg}
    , isSupplementary{isSupplementaryArg}
    , isSecondary{isSecondaryArg}
{}

AlignmentResult::AlignmentResult(bool rReversedArg, int64_t rStartArg, int64_t rEndArg,
                                 int64_t qStartArg, int64_t qEndArg, int64_t qLenArg,
                                 Data::Cigar cigarArg)
    : AlignmentResult{0,       rReversedArg,        rStartArg, rEndArg, qStartArg, qEndArg,
                      qLenArg, std::move(cigarArg), 60,        0,       true,      false,
                      false}
{}

std::unique_ptr<AlignmentResult> AlignmentResult::Clip(int64_t frontClipQuery,
                                                       int64_t backClipQuery, int64_t frontClipRef,
                                                       int64_t backClipRef) const
{
    if (!isAligned) {
        return nullptr;
    }

    int64_t newQStart = 0;
    int64_t newQEnd = 0;
    int64_t newRStart = 0;
    int64_t newREnd = 0;

    const std::vector<unsigned char> alnVec = ConvertCigarToEdlibAln(cigar);
    size_t vecStart = 0;
    size_t vecEnd = alnVec.size();

    int64_t qPos = 0;
    int64_t rPos = 0;

    // Find the front clipping.
    qPos = !rReversed ? qStart : (qLen - qEnd);
    rPos = rStart;
    for (size_t vecId = 0; vecId < alnVec.size(); ++vecId) {
        vecStart = vecId;
        const char op = LOOKUP_EDLIB_TO_CHAR[alnVec[vecId]];
        if ((qPos >= frontClipQuery) && (rPos >= frontClipRef) && op != 'D') {
            break;
        }
        switch (op) {
            case '=':
            case 'X':
                ++qPos;
                ++rPos;
                break;
            case 'I':
                ++qPos;
                break;
            case 'D':
                ++rPos;
                break;
            default:
                break;
        }
    }
    newQStart = qPos;
    newRStart = rPos - frontClipRef;

    // Find the back clipping.
    qPos = (!rReversed ? qEnd : (qLen - qStart)) - 1;
    rPos = rEnd - 1;
    for (int64_t vecId = static_cast<int64_t>(alnVec.size()) - 1; vecId >= 0; --vecId) {
        vecEnd = vecId;
        const char op = LOOKUP_EDLIB_TO_CHAR[alnVec[vecId]];
        if (qPos < backClipQuery && rPos < backClipRef && op != 'D') {
            break;
        }
        switch (op) {
            case '=':
            case 'X':
                --qPos;
                --rPos;
                break;
            case 'I':
                --qPos;
                break;
            case 'D':
                --rPos;
                break;
            default:
                break;
        }
    }
    newQEnd = qPos + 1;
    newREnd = rPos + 1 - frontClipRef;
    ++vecEnd;

    if ((newQEnd - newQStart) <= 0 || (newREnd - newRStart) <= 0 || vecEnd <= vecStart) {
        return nullptr;
    }

    auto newCigar = ConvertEdlibToCigar(
        std::vector<unsigned char>(alnVec.begin() + vecStart, alnVec.begin() + vecEnd));

    if (rReversed) {
        std::swap(newQStart, newQEnd);
        newQStart = qLen - newQStart;
        newQEnd = qLen - newQEnd;
    }

    newQStart -= frontClipQuery;
    newQEnd -= frontClipQuery;

    return std::make_unique<AlignmentResult>(rId, rReversed, newRStart, newREnd, newQStart, newQEnd,
                                             qLen, newCigar, mapq, as, isAligned, isSupplementary,
                                             isSecondary);
}

bool AlignmentResult::operator==(const AlignmentResult& op) const noexcept
{
    return std::tie(rId, rReversed, rStart, rEnd, qStart, qEnd, qLen, cigar, mapq, as, isAligned,
                    isSupplementary, isSecondary) ==
           std::tie(op.rId, op.rReversed, op.rStart, op.rEnd, op.qStart, op.qEnd, op.qLen, op.cigar,
                    op.mapq, op.as, op.isAligned, op.isSupplementary, op.isSecondary);
}

std::ostream& operator<<(std::ostream& out, const AlignmentResult& a)
{
    out << "query" << '\t' << a.qLen << '\t' << a.qStart << '\t' << a.qEnd << '\t'
        << (a.rReversed ? '-' : '+') << '\t' << a.rId << '\t' << a.rStart << '\t' << a.rEnd << '\t'
        << static_cast<int>(a.mapq) << '\t' << a.as << '\t' << a.isAligned << '\t'
        << a.isSupplementary << '\t' << a.isSecondary << '\t' << a.cigar.ToStdString();
    return out;
}

BAM::BamRecord AlnToBam(const int32_t refId, const BAM::BamHeader& header,
                        const AlignmentResult& aln, const BAM::BamRecord& read, const bool ccs)
{
    BAM::BamRecord record{header};
    record.Impl().SetSequenceAndQualities(read.Sequence());
    record.Impl().Name(read.FullName());
    record.Impl().Tags(read.Impl().Tags());
    std::string cigarStr = aln.cigar.ToStdString();
    int clipStart = aln.rReversed ? aln.qLen - aln.qEnd : aln.qStart;
    int clipEnd = aln.rReversed ? aln.qStart : aln.qLen - aln.qEnd;
    if (clipStart != 0) cigarStr = std::to_string(clipStart) + 'S' + cigarStr;
    if (clipEnd != 0) cigarStr += std::to_string(clipEnd) + 'S';
    record.Map(refId, aln.rStart, aln.rReversed ? BAM::Strand::REVERSE : BAM::Strand::FORWARD,
               BAM::Cigar::FromStdString(cigarStr), aln.mapq);

    const int32_t readstart = ccs ? 0 : read.QueryStart();

    if (!aln.rReversed) {
        record.Clip(BAM::ClipType::CLIP_TO_QUERY, readstart + clipStart,
                    readstart + read.Impl().SequenceLength() - clipEnd, true);
    } else {
        record.Clip(BAM::ClipType::CLIP_TO_QUERY, readstart + clipEnd,
                    readstart + read.Impl().SequenceLength() - clipStart, true);
    }
    return record;
}

}  // namespace PacBio
