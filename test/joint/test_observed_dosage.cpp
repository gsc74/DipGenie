#include "ObservedDosage.hpp"
#include "PresenceDP.hpp"
#include "ConstantPresence.hpp"
#include <cassert>
#include <iostream>
#include <random>
using namespace hic;
using namespace hic_calibration;
int main(){
    Calibration c{40,200,10};ObservedDosage generating;generating.calibration=c;
    generating.prior={.06L,.24L,.70L};generating.error_rate=1.2;
    std::map<uint64_t,uint64_t> histogram;
    for(uint64_t n=1;n<300;++n){auto values=generating.components(n);long double probability=0;
        for(unsigned d=0;d<3;++d)probability+=generating.prior[d]*std::exp(values[d]);
        const auto count=uint64_t(std::llround(1000000*probability));if(count)histogram[n]=count;}
    const auto fit=ObservedDosage::fit(c,histogram);
    assert(fit.iterations>1&&fit.iterations<=500&&std::fabs(fit.error_rate-1.2)<.02L);
    for(unsigned d=0;d<3;++d)assert(std::fabs(fit.prior[d]-generating.prior[d])<.005L);
    assert(fit.posterior(1)[0]>.99L&&fit.posterior(40)[2]>.8L);
    long double likelihood=0;
    for(auto[n,count]:histogram){auto values=fit.components(n);long double probability=0;
        for(unsigned d=0;d<3;++d)probability+=fit.prior[d]*std::exp(values[d]);
        likelihood+=count*std::log(probability);
        const auto p=fit.posterior(n);const auto q=dosage_units(p);
        assert(q[0]+q[1]+q[2]==1024);for(unsigned d=0;d<3;++d)assert(std::fabs(q[d]/1024.L-p[d])<=1/1024.L);
    }
    assert(std::fabs(likelihood-fit.log_likelihood)<1e-6L);
    // Independent normalization and high-count regression. A heavier-tailed
    // error distribution must never steal arbitrarily large true-copy counts.
    for(long double shape:{10.L,static_cast<long double>(INFINITY)}) {
        ObservedDosage ordered;ordered.calibration={40,200,shape};ordered.error_rate=20;
        long double total=0;for(uint64_t n=1;n<500;++n)total+=std::exp(ordered.error_log(n));
        assert(std::fabs(total-1)<1e-12L);
        for(uint64_t n:{1000ULL,1000000ULL})assert(ordered.posterior(n)[0]<1e-10L);
    }
    for(long double rate:{0.L,1e-12L,.2L,1.L,20.L}) {
        ObservedDosage ordered;ordered.error_rate=rate;
        const auto recovered=ObservedDosage::fit_error_rate(ordered.error_positive_mean(),30);
        assert(std::fabs(recovered-rate)<1e-10L);
    }
    ObservedDosage point;point.error_rate=0;assert(point.error_log(1)==0&&point.error_log(2)==-INFINITY);
    assert(ObservedDosage::fit_error_rate(100,20)==20);
    const auto empty=ObservedDosage::fit(c,{});assert(!empty.observations&&!empty.iterations);
    for(auto p:{std::array<long double,3>{1,0,0},{0,1,0},{0,0,1},{1.L/3,1.L/3,1.L/3}}){auto q=dosage_units(p);assert(q[0]+q[1]+q[2]==1024);}
    for(const auto&hist:{std::map<uint64_t,uint64_t>{{1,100}},{{40,100}},{{1,1},{1000,1}}}){
        auto f=ObservedDosage::fit(c,hist);for(auto[n,count]:hist){(void)count;auto p=f.posterior(n);assert(std::isfinite(p[0]+p[1]+p[2]));}}
    bool bad=false;try{(void)fit.posterior(0);}catch(const std::invalid_argument&){bad=true;}assert(bad);
    std::mt19937_64 rng(17);
    for(int replicate=0;replicate<100;++replicate){
        PresenceModel m;m.H=2;m.blocks=3;m.budget=replicate%5;m.weights={3,2};
        for(int x=0;x<7;++x){PresenceTarget t;t.hom=rng()%2;t.reward_units=rng()%11;t.other_reward_units=rng()%11;t.absent_reward_units=rng()%11;
            for(int b=0;b<3;++b)if(rng()%2)t.occurrences.push_back({b,DonorMask(1+rng()%3)});
            if(t.occurrences.empty())t.occurrences={{0,DonorMask(1)}};
            m.targets.push_back(t);}
        m.contacts={{0,1,7,1,false},{1,2,4,-1,true},{3,5,9,1,false}};
        Score best=impossible;
        for(int a=0;a<4;++a)for(int b=0;b<4;++b)for(int d=0;d<4;++d){
            if(m.charge(a,b)+m.charge(b,d)>m.budget)continue;
            best=std::max(best,m.score({a,b,d}));}
        for(int threads:{1,2,4}){
            const auto exact=solve_presence(m,0,threads);assert(exact.score==best&&m.score(exact.path)==best);
            const auto pruned=solve_presence(m,3,threads,{},true);
            assert(pruned.score<=best&&pruned.upper>=best&&m.score(pruned.path)==pruned.score);
            const auto folded=fold_constant_presence(m);const auto r=solve_presence(folded.model,0,threads);
            assert(r.score+folded.offset==best&&best<=m.universal_bound());
        }
    }
    std::cout<<"PASS three-class fitted count model, independent likelihood replay, posterior quantization, degenerate inputs, exhaustive joint-DP optimum, beam bounds, constant folding and thread parity\n";
}
