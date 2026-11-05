#pragma once
#include "InferenceTypes.hpp"

namespace hic {
// All rewards use a common factor of two so the neutral panel baseline is
// exact. This is an integer representation, not an objective tuning weight.
constexpr Score phase_reward_scale=2;
inline int phase_orientation(unsigned char presence) {
    return int(presence&1)-int((presence>>1)&1);
}
// neutral=true: 2 unless both endpoints are heterozygous AND discordant with
// the linked relation, then 0. This removes the positive concordance bonus
// relative to an undefined orientation, but does NOT decouple phase from
// genotype: removing a discordant heterozygous endpoint can remove its penalty.
// Free haplotype exchange changes orientation without recombination charge.
// presence=true (linked presence): a linked molecule shows that some haplotype
// carries both alleles (relation +1) or that one haplotype carries each
// (relation -1). The factor is 2 when the diploid pair can have emitted it and
// 0 otherwise: homozygous and heterozygous-concordant pairs score alike (no
// heterozygosity bonus), a missing endpoint scores 0 (an observed allele is
// not dropped for free). Bounded in [0,2] like the other factors.
inline bool linked_presence(unsigned char a,unsigned char b,int relation) {
    return relation>0?(a&b&3)!=0:((a&1)&&(b&2))||((a&2)&&(b&1));
}
inline Score phase_factor(unsigned char a,unsigned char b,int relation,bool /*panel*/,bool neutral=false,bool presence=false) {
    if(presence)return linked_presence(a,b,relation)?2:0;
    const int x=phase_orientation(a),y=phase_orientation(b);
    if(neutral)return relation*x*y<0?0:2;
    // Panel evidence is conditional on heterozygosity. Undefined orientation
    // is neutral (1/2), rather than a failed phase observation (0). For an
    // unordered, uniformly oriented pair the expected panel score is 1/2,
    // independently of endpoint dosage. This is a bounded reward, not a
    // calibrated whole-molecule likelihood.
    // This is a phase score, not an observation likelihood: dosage is scored
    // separately. A homozygous/absent endpoint has no orientation and gets
    // neutral phase reward rather than being treated as discordant evidence.
    return 1+relation*x*y;
}
}
