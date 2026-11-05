#pragma once
#include "CompleteDosage.hpp"
#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <vector>

namespace hic_calibration {
// Copy-number evidence for a string confined to one common-anchor segment
// whose panel paths carry different numbers of copies (e.g. tandem-repeat
// units). Presence alone cannot distinguish such alleles; read counts can.
//
// Encoding in the unchanged three-state objective: level k (k = 1..K) is the
// event "this haplotype carries >= k copies". For a donor pair with copies
// (i,j), its level dosage is d_k=[i>=k]+[j>=k], and sum_k d_k = i+j. Each
// level is scored by its marginal posterior P(d_k | count), where the joint
// posterior over ordered pairs (i,j) has the panel-frequency prior n_i*n_j
// and the calibrated complete-count likelihood of total copy number i+j.
// Every level reward is a probability vector, as for single-copy targets, so
// rewards stay bounded, nonnegative and pair-local within one segment.
inline std::vector<std::array<long double,3>> copy_level_posterior(
        const CompleteDosage&model,uint64_t count,const std::map<uint32_t,uint64_t>&classes) {
    if(classes.empty())throw std::invalid_argument("Copy classes");
    uint32_t K=0;uint64_t paths=0;
    for(auto[c,n]:classes){if(!n)throw std::invalid_argument("Empty copy class");K=std::max(K,c);paths+=n;}
    if(!K||K>255)throw std::invalid_argument("Copy level range");
    std::vector<std::pair<std::pair<uint32_t,uint32_t>,long double>> joint;
    long double top=-INFINITY;
    for(auto[i,a]:classes)for(auto[j,b]:classes) {
        const long double v=std::log(static_cast<long double>(a))+std::log(static_cast<long double>(b))+model.copies_log(count,i+j);
        joint.push_back({{i,j},v});if(v>top)top=v;
    }
    if(!std::isfinite(top))throw std::overflow_error("Copy posterior");
    for(auto&e:joint)e.second=std::exp(e.second-top);
    std::vector<std::array<long double,3>> levels(K,{0,0,0});
    for(const auto&[pair,w]:joint)for(uint32_t k=1;k<=K;++k)
        levels[k-1][unsigned(pair.first>=k)+unsigned(pair.second>=k)]+=w;
    for(auto&p:levels){
        const long double sum=p[0]+p[1]+p[2];
        for(auto&v:p)v=std::min(1.L,std::max(0.L,v/sum));
    }
    (void)paths;
    return levels;
}
}
