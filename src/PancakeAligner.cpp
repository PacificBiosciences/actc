#include "PancakeAligner.hpp"

#include <pbbam/BamRecord.h>

#include <cstdint>
#include <string>
#include <vector>

namespace PacBio {

std::vector<AlnResults> PancakeAligner(Pancake::MapperCLR& mapper,
                                       const std::vector<BAM::BamRecord>& reads,
                                       const std::string& reference)
{
    if (reads.empty()) {
        return {};
    }

    if (reference.empty()) {
        return {};
    }

    // Prepare the target for mapping.
    std::vector<std::string> refs = {reference};

    // Prepare the query sequences for mapping.
    std::vector<std::string> queries;
    for (const auto& r : reads) {
        queries.emplace_back(r.Sequence());
    }

    auto mappingResults = mapper.MapAndAlign(refs, queries);

    std::vector<AlnResults> ret(queries.size());

    // Convert the results.
    for (size_t i = 0; i < mappingResults.size(); ++i) {
        for (size_t j = 0; j < mappingResults[i].mappings.size(); ++j) {
            auto& m = mappingResults[i].mappings[j];
            auto& aln = m->mapping;
            if (aln->Brev) {
                std::reverse(aln->Cigar.begin(), aln->Cigar.end());
            }
            auto alnResult = std::make_unique<AlignmentResult>(
                aln->Bid, aln->Brev, aln->BstartFwd(), aln->BendFwd(), aln->Astart, aln->Aend,
                aln->Alen, std::move(aln->Cigar), 60, aln->Score, true, m->isSupplementary,
                m->priority > 0);
            ret[i].emplace_back(std::move(alnResult));
        }
    }

    return ret;
}

Pancake::MapperCLRMapSettings InitPancakeMapSettingsSubread(const bool shortInsert)
{
    Pancake::MapperCLRMapSettings settings;

    settings.seedParams.KmerSize = 15;
    settings.seedParams.MinimizerWindow = 5;
    settings.seedParams.Spacing = 0;
    settings.seedParams.UseHPCForSeedsOnly = true;

    settings.secondaryAllowedOverlapFractionQuery = 0.0;
    settings.secondaryAllowedOverlapFractionTarget = 0.5;

    settings.seedParamsFallback = settings.seedParams;
    settings.seedParamsFallback.KmerSize = 10;
    settings.seedParamsFallback.MinimizerWindow = 5;

    // Filter out the top percentile of most frequent minimizers.
    settings.freqPercentile = 0.000;
    // Determine the maximum occurrence cutoff automatically from the histogram so that all
    // seed hits for a query can fit into this much memory. If <= 0, it's turned off.
    settings.seedOccurrenceMaxMemory = 100'000'000;
    // Limit the maximum occurrence of a seed to this. (Upper bound.) If <= 0, it's turned off.
    settings.seedOccurrenceMax = 1000;
    // Do not filter seeds with occurrences lower than this. (Lower bound.)
    settings.seedOccurrenceMin = 10;

    settings.seedJoinDist = 10000;
    settings.maxFlankExtensionDist = settings.seedJoinDist;
    settings.minAlignmentSpan = 200;

    if (shortInsert) {
        settings.seedParams.KmerSize = 4;
        settings.seedParams.MinimizerWindow = 1;
        settings.seedParamsFallback = settings.seedParams;
        settings.minAlignmentSpan = 0;
        settings.minDPScore = 10;
        settings.minNumSeeds = 2;
        settings.minQueryLen = 0;
    }

    return settings;
}

Pancake::MapperCLRAlignSettings InitPancakeAlignSettingsSubread()
{
    Pancake::MapperCLRAlignSettings settings;

    settings.alnParamsGlobal.zdrop = 400;
    settings.alnParamsGlobal.zdrop2 = 200;
    settings.alnParamsGlobal.alignBandwidth = 500;
    settings.alnParamsGlobal.endBonus = 1000;
    settings.alnParamsGlobal.matchScore = 2;
    settings.alnParamsGlobal.mismatchPenalty = 4;
    settings.alnParamsGlobal.gapOpen1 = 4;
    settings.alnParamsGlobal.gapExtend1 = 2;
    settings.alnParamsGlobal.gapOpen2 = 24;
    settings.alnParamsGlobal.gapExtend2 = 1;
    settings.alignerTypeGlobal = Pancake::AlignerType::KSW2;
    settings.alignerTypeExt = Pancake::AlignerType::KSW2;
    settings.alnParamsExt = settings.alnParamsGlobal;

    return settings;
}

Pancake::MapperCLRSettings InitPancakeSettingsSubread(const bool shortInsert)
{
    Pancake::MapperCLRSettings settings;
    settings.map = InitPancakeMapSettingsSubread(shortInsert);
    settings.align = InitPancakeAlignSettingsSubread();

    return settings;
}

std::vector<AlnResults> PancakeAlignerSubread(const std::vector<BAM::BamRecord>& reads,
                                              const std::string& reference)
{
    Pancake::MapperCLRSettings settings =
        InitPancakeSettingsSubread(static_cast<int32_t>(reference.size()) < 200);
    Pancake::MapperCLR mapper(settings);
    return PancakeAligner(mapper, reads, reference);
}
}  // namespace PacBio
