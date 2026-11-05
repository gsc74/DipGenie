#pragma once
#include "PresenceDP.hpp"

namespace hic {
struct ConstantPresenceReduction {
    PresenceModel model;
    Score offset=0;
    size_t fixed_targets=0,fixed_contacts=0;
};
// If every donor emits x in even one segment, every complete diploid pair
// contains x on both paths. Its dosage reward is constant. Every incident
// phase factor is neutral because its orientation is zero. Keep the original
// normalization weights: removing constants is not re-estimating the model.
inline ConstantPresenceReduction fold_constant_presence(const PresenceModel&m) {
    m.check();ConstantPresenceReduction out;
    out.model.H=m.H;out.model.blocks=m.blocks;out.model.budget=m.budget;out.model.weights=m.weights;out.model.free_swaps=m.free_swaps;out.model.neutral_phase=m.neutral_phase;out.model.presence_phase=m.presence_phase;out.model.ordered_repair=m.ordered_repair;
    std::vector<size_t> remap(m.targets.size(),SIZE_MAX);
    const auto full=DonorMask::full(size_t(m.H));
    // A constant endpoint has zero orientation, but its linked-presence factor
    // still depends on the other endpoint: keep such targets in the model.
    std::vector<char> incident(m.targets.size(),0);
    if(m.presence_phase)for(const auto&e:m.contacts)incident[e.x]=incident[e.y]=1;
    for(size_t x=0;x<m.targets.size();++x) {
        const auto&t=m.targets[x];bool constant=false;
        for(const auto&entry:t.occurrences)if(entry.second==full){constant=true;break;}
        if(constant&&!incident[x]) {
            ++out.fixed_targets;
            out.offset=score_sum(out.offset,score_product(score_product(phase_reward_scale,
                PresenceModel::genotype_reward(t,3)),m.weights.genotype));
        } else {
            remap[x]=out.model.targets.size();out.model.targets.push_back(t);
        }
    }
    for(const auto&e:m.contacts) {
        if(remap[e.x]==SIZE_MAX||remap[e.y]==SIZE_MAX) {
            ++out.fixed_contacts;
            // A constant endpoint has no orientation: factor 1 (v1) or 2 (neutral).
            out.offset=score_sum(out.offset,score_product(score_product(e.units,m.neutral_phase?2:1),m.weights.contact));
        } else out.model.contacts.push_back({remap[e.x],remap[e.y],e.units,e.relation,e.panel});
    }
    out.model.check();return out;
}
}
