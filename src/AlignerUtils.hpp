#pragma once

#include <pbcopper/data/Cigar.h>
#include <pbcopper/data/CigarOperation.h>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace PacBio {

/*
 * The following was copied from pbbam/CigarOperationType.h and modified.
*/
inline int32_t CigarCharToOp(const char c)
{
    switch (c) {
        case 'M':
            return 0;
        case 'I':
            return 1;
        case 'D':
            return 2;
        case 'N':
            return 3;
        case 'S':
            return 4;
        case 'H':
            return 5;
        case 'P':
            return 6;
        case '=':
            return 7;
        case 'X':
            return 8;
        default:
            return -1;
    }
}

/*
 * The following was copied from pbbam/CigarOperationType.h and modified.
*/
inline unsigned char CigarCharToEdlib(const char c)
{
    switch (c) {
        case 'M':
            return 0;
        case 'I':
            return 1;
        case 'D':
            return 2;
        case 'N':
            return 4;
        case 'S':
            return 4;
        case 'H':
            return 4;
        case 'P':
            return 4;
        case '=':
            return 0;
        case 'X':
            return 3;
        default:
            return 4;
    }
}

/*
 * The following was copied from pbbam/CigarOperationType.h and modified.
*/
inline int32_t CigarTypeToOp(Data::CigarOperationType c)
{
    switch (c) {
        case Data::CigarOperationType::ALIGNMENT_MATCH:
            return 0;
        case Data::CigarOperationType::INSERTION:
            return 1;
        case Data::CigarOperationType::DELETION:
            return 2;
        case Data::CigarOperationType::REFERENCE_SKIP:
            return 3;
        case Data::CigarOperationType::SOFT_CLIP:
            return 4;
        case Data::CigarOperationType::HARD_CLIP:
            return 5;
        case Data::CigarOperationType::PADDING:
            return 6;
        case Data::CigarOperationType::SEQUENCE_MATCH:
            return 7;
        case Data::CigarOperationType::SEQUENCE_MISMATCH:
            return 8;
        default:
            return -1;
    }
}

// Cigar ops                                            M, I, D, N, S, H, P, =, X, B
// Minimap2 (Same as Data::Cigar)                        0, 1, 2, 3, 4, 5, 6, 7, 8, 9
// Edlib                                                0, 1, 2, -, -, -, -, 0, 3, -
constexpr std::array<int32_t, 10> LOOKUP_MM2_TO_EDLIB{{0, 1, 2, 4, 4, 4, 4, 0, 3, 8}};
constexpr std::array<int32_t, 5> LOOKUP_EDLIB_TO_MM2{
    {7, 1, 2, 8, -1}};  // Everything > 3 is not defined by Edlib.
constexpr std::array<char, 5> LOOKUP_EDLIB_TO_CHAR{
    {'=', 'I', 'D', 'X', '?'}};  // '?' is a dummy value here, everything above 3 is undefined.

std::vector<unsigned char> ConvertCigarToEdlibAln(const Data::Cigar& cigar);

Data::Cigar ConvertEdlibToCigar(const std::vector<unsigned char>& aln);

double CalcAlignmentIdentity(const Data::Cigar& cigar);

bool ConvertCigarToM5(const std::string& ref, const std::string& query, int32_t rStart,
                      int32_t rEnd, int32_t qStart, int32_t qEnd, bool qRev,
                      const Data::Cigar& cigar, std::string& retRefAln, std::string& retQueryAln);

}  // namespace PacBio
