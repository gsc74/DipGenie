#pragma once
#include "DosageConfidence.hpp"
#include <array>

namespace hic_calibration {
// Conditional model for observed (>0 count), unique single-locus strings.
// D=0 uses a zero-truncated Poisson error count; D=1/2 use calibrated
// negative-binomial counts conditioned on being positive. Class proportions
// and the error rate are fitted from input counts only, never truth/output.
// The error rate cannot exceed the one-copy mean. Unlike a geometric error
// tail against Poisson signal, this cannot explain high counts as absence.
struct ObservedDosage {
    Calibration calibration{};
    std::array<long double,3> prior{.1L,.45L,.45L};
    long double error_rate=.5L,log_likelihood=0;
    size_t observations=0,iterations=0;
    long double error_positive_mean()const {
        return error_rate==0?1.L:error_rate/(-std::expm1(-error_rate));
    }
    long double error_log(uint64_t n)const {
        if(!n)throw std::invalid_argument("Observed dosage requires positive count");
        if(error_rate<0||!std::isfinite(error_rate))throw std::invalid_argument("Error count rate");
        if(error_rate==0)return n==1?0.L:-INFINITY;
        const auto normalization=error_rate>1?
            error_rate+std::log1p(-std::exp(-error_rate)):std::log(std::expm1(error_rate));
        return static_cast<long double>(n)*std::log(error_rate)-std::lgamma(static_cast<long double>(n)+1)-normalization;
    }
    static long double fit_error_rate(long double mean,long double cap) {
        if(!(mean>=1)||!std::isfinite(mean)||!(cap>0)||!std::isfinite(cap))
            throw std::invalid_argument("Error count mean/rate cap");
        if(mean==1)return 0;
        if(mean>=cap/(-std::expm1(-cap)))return cap;
        long double lo=0,hi=std::min(mean,cap);
        for(unsigned i=0;i<100;++i){
            const auto mid=(lo+hi)/2;
            if(mid/(-std::expm1(-mid))<mean)lo=mid;else hi=mid;
        }
        return (lo+hi)/2;
    }
    long double nb_log(uint64_t n,unsigned copies)const {
        const auto value=copy_log_likelihood(calibration,n,copies);
        const long double r=calibration.shape;
        return value-std::lgamma(static_cast<long double>(n)+1)+
            (std::isinf(r)?0:std::lgamma(static_cast<long double>(n)+r)-std::lgamma(r));
    }
    std::array<long double,3> components(uint64_t n)const {
        if(!n)throw std::invalid_argument("Observed dosage requires positive count");
        std::array<long double,3> result;
        result[0]=error_log(n);
        for(unsigned d=1;d<=2;++d)result[d]=nb_log(n,d)-std::log(-std::expm1(nb_log(0,d)));
        return result;
    }
    std::array<long double,3> posterior(uint64_t n)const {
        auto p=components(n);long double maximum=-INFINITY,sum=0;
        for(unsigned d=0;d<3;++d){p[d]+=std::log(prior[d]);maximum=std::max(maximum,p[d]);}
        for(auto&v:p){v=std::exp(v-maximum);sum+=v;}
        if(!(sum>0)||!std::isfinite(sum))throw std::overflow_error("Observed dosage posterior");
        for(auto&v:p)v/=sum;
        return p;
    }
    static ObservedDosage fit(const Calibration&c,const std::map<uint64_t,uint64_t>&histogram) {
        ObservedDosage model;model.calibration=c;
        (void)model.nb_log(0,1); // validates calibration even for an empty histogram
        const long double cap=c.mean/2;
        model.error_rate=std::min(model.error_rate,cap);
        struct Bin{uint64_t n,frequency;long double one,two;};std::vector<Bin> bins;
        for(auto[n,f]:histogram){
            if(!n||!f||f>SIZE_MAX-model.observations)throw std::invalid_argument("Observed dosage histogram");
            model.observations+=f;
            auto p=model.components(n);bins.push_back({n,f,p[1],p[2]});
        }
        if(bins.empty())return model;
        long double previous=-INFINITY;
        for(size_t iteration=0;iteration<500;++iteration){
            std::array<long double,3> mass{};long double count=0,ll=0;
            for(const auto&b:bins){
                const long double error=model.error_log(b.n);
                std::array<long double,3> p{error+std::log(model.prior[0]),b.one+std::log(model.prior[1]),b.two+std::log(model.prior[2])};
                const auto maximum=*std::max_element(p.begin(),p.end());long double total=0;
                for(auto&v:p){v=std::exp(v-maximum);total+=v;}
                ll+=b.frequency*(maximum+std::log(total));
                for(unsigned d=0;d<3;++d)mass[d]+=b.frequency*p[d]/total;
                count+=static_cast<long double>(b.frequency)*b.n*p[0]/total;
            }
            if(!std::isfinite(ll)||ll<previous-1e-10L*std::max(1.L,std::fabs(previous)))
                throw std::logic_error("Observed dosage EM likelihood decreased/nonfinite");
            model.log_likelihood=ll;model.iterations=iteration+1;
            if((iteration&&ll-previous<1e-9L*std::max(1.L,std::fabs(ll)))||iteration+1==500)break;
            previous=ll;long double sum=0;
            for(unsigned d=0;d<3;++d){model.prior[d]=std::max(1e-15L,mass[d]/model.observations);sum+=model.prior[d];}
            for(auto&v:model.prior)v/=sum;
            if(mass[0]>0)model.error_rate=fit_error_rate(std::max(1.L,count/mass[0]),cap);
        }
        return model;
    }
};
// Largest-remainder quantization preserves exactly one fixed-point unit of
// probability over all three dosages, including ties and near-zero classes.
inline std::array<int64_t,3> dosage_units(const std::array<long double,3>&p,int64_t scale=1024) {
    if(scale<=0)throw std::invalid_argument("Dosage probability scale");
    std::array<int64_t,3> units{};std::array<long double,3> residual{};int64_t total=0;long double sum=0;
    for(unsigned d=0;d<3;++d){
        if(!std::isfinite(p[d])||p[d]<0||p[d]>1)throw std::invalid_argument("Dosage probability");
        sum+=p[d];const auto v=p[d]*scale;units[d]=int64_t(std::floor(v));residual[d]=v-units[d];total+=units[d];
    }
    if(std::fabs(sum-1)>1e-12L||total>scale||scale-total>3)throw std::invalid_argument("Dosage probability normalization");
    while(total<scale){auto it=std::max_element(residual.begin(),residual.end());++units[size_t(it-residual.begin())];*it=-1;++total;}
    return units;
}
// Clipped log-odds weights: w_d = max(0, cap - log(p_max/p_d)) for the
// posterior p (prior times count likelihood), scaled by scale/(3 cap) so the
// three weights sum to at most scale. The most probable dosage gets cap; a
// dosage that is cap nats less likely gets 0. Unlike posterior probabilities,
// weak evidence (an unobserved string at low depth) costs little, while an
// observed string weighs its full log-likelihood ratio.
inline std::array<int64_t,3> log_odds_units(const std::array<long double,3>&p,long double cap,int64_t scale=1024) {
    if(scale<=0||!(cap>0))throw std::invalid_argument("Log-odds dosage scale");
    long double top=-INFINITY;std::array<long double,3> l{};
    for(unsigned d=0;d<3;++d){
        if(!std::isfinite(p[d])||p[d]<0||p[d]>1)throw std::invalid_argument("Dosage probability");
        l[d]=p[d]>0?std::log(p[d]):-INFINITY;top=std::max(top,l[d]);
    }
    if(top==-INFINITY)throw std::invalid_argument("Dosage probability normalization");
    std::array<int64_t,3> units{};
    for(unsigned d=0;d<3;++d){
        const long double w=l[d]==-INFINITY?0.L:std::max(0.L,cap-(top-l[d]));
        units[d]=int64_t(std::llround(w*scale/(3*cap)));
    }
    return units;
}
}
