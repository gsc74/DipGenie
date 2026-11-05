#pragma once
#include "ObservedDosage.hpp"

namespace hic_calibration {
// A normalized count distribution over ALL local panel strings, including
// zero observations. D=0 has a fitted zero mass and positive Poisson tail.
// D=1/2 use the calibrated NB body plus excess zero mass estimated exclusively
// from conserved controls. No inferred path or benchmark label enters the fit.
// Positive error counts reuse the tested truncated-Poisson calculation.
struct CompleteDosage : ObservedDosage {
    long double error_zero=.95L,dropout=0;
    uint64_t zero_observations=0;
    long double absence_log(uint64_t n)const {
        return n?std::log1p(-error_zero)+ObservedDosage::error_log(n):std::log(error_zero);
    }
    long double present_log(uint64_t n,unsigned d)const {
        const auto body=std::log1p(-dropout)+nb_log(n,d);
        return n?body:log_add(std::log(dropout),body);
    }
    std::array<long double,3> components(uint64_t n)const {
        return {absence_log(n),present_log(n,1),present_log(n,2)};
    }
    // Count likelihood for an arbitrary total diploid copy number.
    long double copies_log(uint64_t n,unsigned copies)const {
        return copies?present_log(n,copies):absence_log(n);
    }
    std::array<long double,3> posterior(uint64_t n)const {return posterior(n,prior);}
    // Posterior under a caller-supplied dosage prior (e.g. from the panel
    // frequency of the string); the count likelihoods are unchanged.
    std::array<long double,3> posterior(uint64_t n,const std::array<long double,3>&with)const {
        auto p=components(n);long double maximum=-INFINITY,sum=0;
        for(unsigned d=0;d<3;++d){
            if(!(with[d]>0))throw std::invalid_argument("Dosage prior must be positive");
            p[d]+=std::log(with[d]);maximum=std::max(maximum,p[d]);
        }
        for(auto&v:p){v=std::exp(v-maximum);sum+=v;}
        if(!(sum>0)||!std::isfinite(sum))throw std::overflow_error("Complete dosage posterior");
        for(auto&v:p)v/=sum;
        return p;
    }
    static CompleteDosage fit(const Calibration&c,const std::map<uint64_t,uint64_t>&histogram,
                             uint64_t control_zero,uint64_t control_total) {
        if(!control_total||control_zero>=control_total)
            throw std::invalid_argument("Complete dosage requires observed conserved controls");
        CompleteDosage model;model.calibration=c;
        const auto positive=-std::expm1(model.nb_log(0,2));
        const auto fraction=static_cast<long double>(control_zero)/control_total;
        model.dropout=std::max(0.L,1-(1-fraction)/positive);
        if(!(model.dropout>=0&&model.dropout<1))throw std::invalid_argument("Count dropout estimate");
        const long double cap=c.mean/2;
        model.error_rate=std::min(model.error_rate,cap);
        model.prior={.5L,.15L,.35L};
        struct Bin{uint64_t n,frequency;long double one,two;};std::vector<Bin>bins;
        for(auto[n,f]:histogram){
            if(!f||f>SIZE_MAX-model.observations)throw std::invalid_argument("Complete dosage histogram");
            model.observations+=f;if(!n)model.zero_observations=f;
            bins.push_back({n,f,model.present_log(n,1),model.present_log(n,2)});
        }
        if(bins.empty())return model;
        long double previous=-INFINITY;
        for(size_t iteration=0;iteration<500;++iteration){
            std::array<long double,3>mass{};long double count=0,positive_error=0,zero_error=0,ll=0;
            for(const auto&b:bins){
                std::array<long double,3>p{model.absence_log(b.n)+std::log(model.prior[0]),
                    b.one+std::log(model.prior[1]),b.two+std::log(model.prior[2])};
                const auto maximum=*std::max_element(p.begin(),p.end());long double total=0;
                for(auto&v:p){v=std::exp(v-maximum);total+=v;}
                ll+=b.frequency*(maximum+std::log(total));
                for(unsigned d=0;d<3;++d)mass[d]+=b.frequency*p[d]/total;
                const auto error=b.frequency*p[0]/total;
                if(b.n){positive_error+=error;count+=b.n*error;}else zero_error+=error;
            }
            if(!std::isfinite(ll)||(iteration&&ll<previous-1e-10L*std::max(1.L,std::fabs(previous))))
                throw std::logic_error("Complete dosage EM likelihood decreased/nonfinite");
            model.log_likelihood=ll;model.iterations=iteration+1;
            if((iteration&&ll-previous<1e-9L*std::max(1.L,std::fabs(ll)))||iteration+1==500)break;
            previous=ll;long double sum=0;
            for(unsigned d=0;d<3;++d){model.prior[d]=std::max(1e-15L,mass[d]/model.observations);sum+=model.prior[d];}
            for(auto&v:model.prior)v/=sum;
            if(mass[0]>0)model.error_zero=std::clamp(zero_error/mass[0],1e-15L,1-1e-15L);
            if(positive_error>0)model.error_rate=fit_error_rate(std::max(1.L,count/positive_error),cap);
        }
        return model;
    }
};
}
