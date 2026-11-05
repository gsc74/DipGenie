#pragma once
#include "SequenceEvidence.hpp"

namespace hic {
// Conditional phase evidence from two independently drawn panel haplotypes,
// given heterozygosity at both targets. A=n11*n00 and B=n10*n01 count cis
// and trans diplotypes, up to the same factor of two. No truth labels enter.
struct PanelPreference {int relation=0;uint64_t difference=0,total=0;};
inline PanelPreference panel_preference(const DonorMask&x,const DonorMask&y,int H) {
    (void)checked_pair_domain(size_t(H));
    if(!x.within(size_t(H))||!y.within(size_t(H)))throw std::invalid_argument("Panel carrier mask");
    const uint64_t nx=x.count(),ny=y.count(),n11=x.intersection_count(y);
    const uint64_t n00=uint64_t(H)-nx-(ny-n11);
    uint64_t a=n11*n00;
    uint64_t b=(nx-n11)*(ny-n11);
    return {a==b?0:(a>b?1:-1),a>b?a-b:b-a,a+b};
}
struct PanelSummary {Score units=0;size_t boundaries=0,links=0;};
inline PanelSummary add_panel_phase(PresenceModel&m,const SequenceIndex&index,
                                    const std::vector<size_t>&target_ids,
                                    const std::vector<uint64_t>*observed_counts=nullptr) {
    if(target_ids.size()!=m.targets.size())throw std::invalid_argument("Panel target map");
    if(observed_counts&&observed_counts->size()!=target_ids.size())throw std::invalid_argument("Panel support map");
    // One representative per carrier pattern and common-anchor segment.
    // Overlapping strings with the same panel evidence are not extra votes.
    std::map<int,std::map<DonorMask,size_t>> groups;
    const DonorMask all=DonorMask::full(size_t(m.H));
    for(size_t x=0;x<m.targets.size();++x) {
        // Keep the observed representative at each carrier pattern. Adding
        // zero-support dosage targets must not replace its phase evidence.
        if(observed_counts&&!observed_counts->at(x))continue;
        const auto&t=index.targets.at(target_ids[x]);
        if(!t.unique_locus||t.occurrences.size()!=1||!t.carriers||t.carriers==all)continue;
        groups[t.occurrences.front().first].emplace(t.carriers,x);
    }
    PanelSummary result;if(groups.size()<2)return result;
    auto previous=groups.begin();
    for(auto next=std::next(previous);next!=groups.end();previous=next++) {
        const auto&a=previous->second;const auto&b=next->second;
        const uint64_t combinations=uint64_t(a.size())*b.size();
        Score boundary_units=0;
        for(const auto&[mask_x,x]:a)for(const auto&[mask_y,y]:b) {
            const auto p=panel_preference(mask_x,mask_y,m.H);
            if(!p.relation||!p.total)continue;
            // Confidence |A-B|/(A+B), averaged over patterns at this boundary.
            // Conservative fixed-point rounding; total mass <= one molecule.
            const auto het_weight=[&](size_t z){const auto&t=m.targets[z];return t.hom?t.other_reward_units:t.reward_units;};
            const Score units=Score((__int128(p.difference)*het_weight(x)*het_weight(y))/
                                    (__int128(1024)*combinations*p.total));
            if(!units)continue;
            m.contacts.push_back({x,y,units,p.relation,true});
            boundary_units+=units;++result.links;
        }
        if(boundary_units>1024)throw std::logic_error("Panel boundary mass exceeds one");
        if(boundary_units){++result.boundaries;result.units+=boundary_units;}
    }
    return result;
}
}
