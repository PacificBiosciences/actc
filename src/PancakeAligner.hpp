#pragma once

#include "AlignmentResult.hpp"

#include <pacbio/pancake/MapperCLR.h>
#include <pbbam/BamRecord.h>

#include <cstdint>
#include <string>
#include <vector>

namespace PacBio {

std::vector<AlnResults> PancakeAligner(Pancake::MapperCLR& mapper,
                                       const std::vector<BAM::BamRecord>& reads,
                                       const std::string& reference);

Pancake::MapperCLRMapSettings InitPancakeMapSettingsSubread(const bool shortInsert);

Pancake::MapperCLRAlignSettings InitPancakeAlignSettingsSubread();

Pancake::MapperCLRSettings InitPancakeSettingsSubread(const bool shortInsert);

std::vector<AlnResults> PancakeAlignerSubread(const std::vector<BAM::BamRecord>& reads,
                                              const std::string& reference);
}  // namespace PacBio
