#include "AlignerUtils.hpp"

#include <pbcopper/utility/SequenceUtils.h>

#include <cassert>
#include <cstddef>
#include <sstream>
#include <stdexcept>

namespace PacBio {

Data::Cigar ConvertEdlibToCigar(const std::vector<unsigned char>& aln)
{
    Data::Cigar cigar;
    if (aln.empty()) return cigar;

    constexpr uint32_t defaultCount = 1;
    uint32_t streakCount = 1;

    // First op. This removes need for checking if cigar is empty.
    char op = LOOKUP_EDLIB_TO_CHAR[aln[0]];
    cigar.emplace_back(Data::CigarOperation{op, defaultCount});

    for (size_t i = 1; i < aln.size(); ++i) {
        op = LOOKUP_EDLIB_TO_CHAR[aln[i]];
        Data::CigarOperation cig_op(op, defaultCount);

        if (cigar.back().Type() == cig_op.Type()) {
            ++streakCount;
            cigar.back().Length(streakCount);
        } else {
            streakCount = 1;
            cigar.emplace_back(std::move(cig_op));
        }
    }

    return cigar;
}

std::vector<unsigned char> ConvertCigarToEdlibAln(const Data::Cigar& cigar)
{
    std::vector<unsigned char> ret;
    for (const auto& cigar_op : cigar) {
        const int32_t op = CigarTypeToOp(cigar_op.Type());
        const int32_t count = cigar_op.Length();
        const unsigned char edlib_op = LOOKUP_MM2_TO_EDLIB[op];
        ret.insert(ret.end(), count, edlib_op);
    }
    return ret;
}

double CalcAlignmentIdentity(const Data::Cigar& cigar)
{
    int64_t num_eq = 0;
    int64_t num_x = 0;
    int64_t num_ins = 0;
    int64_t num_del = 0;

    for (const auto& cigar_op : cigar) {
        const auto op = cigar_op.Type();
        const int64_t count = static_cast<int64_t>(cigar_op.Length());

        switch (op) {
            case Data::CigarOperationType::SEQUENCE_MATCH:
                num_eq += count;
                break;
            case Data::CigarOperationType::SEQUENCE_MISMATCH:
                num_x += count;
                break;
            case Data::CigarOperationType::INSERTION:
                num_ins += count;
                break;
            case Data::CigarOperationType::DELETION:
                num_del += count;
                break;
            default:
                // PBLOG_DEBUG << "Unknown CIGAR operation encountered. Returning 0.0 instead of identity.\n";
                return 0.0;
        }
    }

    const int64_t qlen = num_eq + num_x + num_ins;
    double idt = (qlen > 0) ? (static_cast<double>(num_eq) / static_cast<double>(qlen)) : 0.0;

    return idt;
}

bool ConvertCigarToM5(const std::string& ref, const std::string& query, const int32_t rStart,
                      const int32_t rEnd, const int32_t qStart, const int32_t qEnd, const bool qRev,
                      const Data::Cigar& cigar, std::string& retRefAln, std::string& retQueryAln)
{
    retRefAln.clear();
    retQueryAln.clear();

    // Sanity check.
    if (cigar.empty()) {
        return false;
    }

#ifndef NDEBUG
    // Calculate the query and reference length from the CIGAR
    // to check for sanity.
    int32_t calcRefLen = 0;
    int32_t calcQueryLen = 0;
    for (const auto& cigarOp : cigar) {
        const auto op = cigarOp.Type();
        if (op == Data::CigarOperationType::ALIGNMENT_MATCH ||
            op == Data::CigarOperationType::SEQUENCE_MATCH ||
            op == Data::CigarOperationType::SEQUENCE_MISMATCH) {
            calcRefLen += cigarOp.Length();
            calcQueryLen += cigarOp.Length();
        } else if (op == Data::CigarOperationType::INSERTION ||
                   op == Data::CigarOperationType::SOFT_CLIP) {
            calcQueryLen += cigarOp.Length();
        } else if (op == Data::CigarOperationType::DELETION ||
                   op == Data::CigarOperationType::REFERENCE_SKIP) {
            calcRefLen += cigarOp.Length();
        }
    }

    assert(calcQueryLen == (qEnd - qStart));
    assert(calcRefLen == (rEnd - rStart));
#endif

    // Prepare the query for simpler usage.
    std::string querySub = query.substr(qStart, qEnd - qStart);
    if (qRev) {
        querySub = Utility::ReverseComplemented(querySub);
    }

    int32_t qPos = 0;
    int32_t rPos = rStart;

    int32_t localPos = 0;

    const int32_t maxReserved = (qEnd - qStart) + (rEnd - rStart);
    retRefAln.resize(maxReserved);
    retQueryAln.resize(maxReserved);

    for (const auto& cigarOp : cigar) {
        const auto op = cigarOp.Type();
        const int32_t count = cigarOp.Length();
        if (op == Data::CigarOperationType::ALIGNMENT_MATCH ||
            op == Data::CigarOperationType::SEQUENCE_MATCH ||
            op == Data::CigarOperationType::SEQUENCE_MISMATCH) {
            for (int32_t i = 0; i < count; ++i, ++qPos, ++rPos, ++localPos) {
                retQueryAln[localPos] = querySub[qPos];
                retRefAln[localPos] = ref[rPos];
            }
        } else if (op == Data::CigarOperationType::INSERTION ||
                   op == Data::CigarOperationType::SOFT_CLIP) {
            for (int32_t i = 0; i < count; ++i, ++qPos, ++localPos) {
                retQueryAln[localPos] = querySub[qPos];
                retRefAln[localPos] = '-';
            }
        } else if (op == Data::CigarOperationType::DELETION ||
                   op == Data::CigarOperationType::REFERENCE_SKIP) {
            for (int32_t i = 0; i < count; ++i, ++rPos, ++localPos) {
                retQueryAln[localPos] = '-';
                retRefAln[localPos] = ref[rPos];
            }
        } else {
            std::string msg{"ERROR: Unknown CIGAR op: "};
            msg += cigarOp.Char();
            throw std::runtime_error{msg};
        }
    }
    retRefAln.resize(localPos);
    retQueryAln.resize(localPos);

    return true;
}

}  // namespace PacBio
