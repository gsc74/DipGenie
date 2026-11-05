#pragma once
#include "ControlCalibration.hpp"
#include <algorithm>
#include <map>
#include <vector>

namespace hic_calibration {
// NB log likelihood up to constants depending only on count and calibration.
// One copy has half the mean of conserved diploid single-locus controls.
inline long double copy_log_likelihood(const Calibration&c,uint64_t count,uint64_t copies) {
    if(!copies||!(c.mean>0)||!std::isfinite(c.mean)||!(c.shape>0))
        throw std::invalid_argument("Copy-number likelihood inputs");
    const long double mean=c.mean*(static_cast<long double>(copies)/2);
    if(!std::isfinite(mean))throw std::overflow_error("Copy-number mean");
    const long double n=count;
    if(std::isinf(c.shape))return n*std::log(mean)-mean;
    // Stable even for very small/large mean/shape ratios.
    const long double u=std::log(mean)-std::log(c.shape);
    const auto softplus=[](long double x){return x>0?x+std::log1p(std::exp(-x)):std::log1p(std::exp(x));};
    return -c.shape*softplus(u)-n*softplus(-u);
}
inline long double log_add(long double a,long double b) {
    if(a==-INFINITY)return b;
    if(b==-INFINITY)return a;
    if(a<b)std::swap(a,b);
    return a+std::log1p(std::exp(b-a));
}
struct DosageConfidence {long double hom_probability,het_probability,log_odds;};
// P(D=d | counts, D in {1,2}), equal dosage-class priors. Within each
// dosage class, total copy number is drawn from two panel haplotypes.
// The SAME latent total copy number applies to every independent library;
// marginalizing it independently in each library would be inconsistent.
inline DosageConfidence dosage_confidence(
    const std::vector<std::pair<uint32_t,uint32_t>>&hist,uint32_t H,
    const std::vector<Calibration>&calibrations,const std::vector<uint64_t>&counts) {
    if(!H||hist.empty()||counts.empty()||counts.size()!=calibrations.size())
        throw std::invalid_argument("Dosage confidence dimensions");
    uint64_t positive=0;uint32_t previous=0;
    for(auto[k,n]:hist) {
        if(!k||!n||k<=previous||n>H||positive+n>H)
            throw std::invalid_argument("Panel copy histogram");
        previous=k;positive+=n;
    }
    if(hist.size()==1&&hist.front().first==1&&positive<H) {
        long double odds=0;
        for(size_t l=0;l<counts.size();++l)odds+=calibrations[l].hom_log_ratio(counts[l]);
        if(!std::isfinite(odds))throw std::overflow_error("Dosage log odds");
        const auto e=std::exp(-std::fabs(odds)),small=e/(1+e),large=1/(1+e);
        return odds>=0?DosageConfidence{large,small,odds}:DosageConfidence{small,large,odds};
    }
    std::map<uint64_t,long double> one,two;
    if(positive<H)for(auto[k,n]:hist)one[k]=n; // conditioning cancels zero-copy partner count
    for(auto[a,n]:hist)for(auto[b,m]:hist)two[uint64_t(a)+b]+=static_cast<long double>(n)*m;
    auto marginal=[&](const auto&mixture) {
        long double total=0,result=-INFINITY;
        for(auto[k,n]:mixture) {
            total+=n;long double value=std::log(n);
            for(size_t l=0;l<counts.size();++l)value+=copy_log_likelihood(calibrations[l],counts[l],k);
            result=log_add(result,value);
        }
        return total?result-std::log(total):-INFINITY;
    };
    const auto l1=marginal(one),l2=marginal(two);
    if(!std::isfinite(l2)||std::isnan(l1))throw std::overflow_error("Dosage likelihood");
    if(l1==-INFINITY)return{1,0,INFINITY};
    const auto odds=l2-l1;
    // Compute both probabilities stably; do not subtract a rounded 1.
    const auto e=std::exp(-std::fabs(odds)),small=e/(1+e),large=1/(1+e);
    return odds>=0?DosageConfidence{large,small,odds}:DosageConfidence{small,large,odds};
}
}
