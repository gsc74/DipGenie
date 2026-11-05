#include "CompleteDosage.hpp"
#include "BuildModel.hpp"
#include "ConstantPresence.hpp"
#include <cassert>
#include <iostream>
using namespace hic_calibration;
using namespace hic;
int main(){
    Calibration c{40,200,10};
    CompleteDosage generating;generating.calibration=c;generating.prior={.5L,.15L,.35L};
    generating.error_zero=.94L;generating.error_rate=1.1L;generating.dropout=.004L;
    std::map<uint64_t,uint64_t>hist;
    for(uint64_t n=0;n<1000;++n){auto logs=generating.components(n);long double p=0;
        for(unsigned d=0;d<3;++d)p+=generating.prior[d]*std::exp(logs[d]);
        auto f=uint64_t(std::llround(p*10000000));if(f)hist[n]=f;}
    const uint64_t controls=10000000,zero=uint64_t(std::llround(controls*std::exp(generating.present_log(0,2))));
    const auto fitted=CompleteDosage::fit(c,hist,zero,controls);
    assert(fitted.zero_observations==hist[0]);
    assert(std::fabs(fitted.error_rate-generating.error_rate)<.02L);
    assert(std::fabs(fitted.error_zero-generating.error_zero)<.002L);
    for(unsigned d=0;d<3;++d)assert(std::fabs(fitted.prior[d]-generating.prior[d])<.003L);
    long double ll=0;
    for(auto[n,f]:hist){auto logs=fitted.components(n);long double sum=0;
        for(unsigned d=0;d<3;++d)sum+=fitted.prior[d]*std::exp(logs[d]);ll+=f*std::log(sum);
        const auto p=fitted.posterior(n);const auto q=dosage_units(p);assert(q[0]+q[1]+q[2]==1024);
    }
    assert(std::fabs(ll-fitted.log_likelihood)<1e-5L);
    assert(fitted.posterior(0)[0]>.99L&&fitted.posterior(40)[2]>.8L);
    for(long double shape:{10.L,static_cast<long double>(INFINITY)})for(long double dropout:{0.L,.1L,.9L}){
        CompleteDosage m;m.calibration={40,200,shape};m.dropout=dropout;
        for(unsigned d=0;d<3;++d){long double mass=0;for(uint64_t n=0;n<1000;++n)mass+=std::exp(m.components(n)[d]);assert(std::fabs(mass-1)<1e-12L);}
        assert(m.posterior(1000000)[0]<1e-10L);
    }
    for(const auto&h:{std::map<uint64_t,uint64_t>{},{{0,1000}},{{1,1000}},{{20,1000}},{{0,100},{1,1},{10000,1}}}){
        const auto m=CompleteDosage::fit(c,h,0,1000);
        for(uint64_t n:{0ULL,1ULL,40ULL,10000ULL}){auto p=m.posterior(n);assert(std::isfinite(p[0]+p[1]+p[2]));}
    }
    bool bad=false;try{(void)CompleteDosage::fit(c,{{0,10}},10,10);}catch(const std::invalid_argument&){bad=true;}assert(bad);
    // A zero count must remain uncertain at shallow coverage, not become a
    // categorical exclusion. Excess control zeros weaken absence evidence.
    CompleteDosage shallow=generating;shallow.calibration={2,2,INFINITY};
    assert(shallow.posterior(0)[0]<.95L);
    CompleteDosage missing=fitted;missing.dropout=.5;
    assert(missing.posterior(0)[0]<fitted.posterior(0)[0]);
    // Raw model integration: zeros alter the actual optimum, with the same
    // recombination domain and valid exact/pruned a posteriori upper bounds.
    SequenceIndex index;
    WordTarget control;control.word=0;control.count=40;control.carriers=3;
    control.occurrences={{0,3}};control.copy_histogram={{1,2}};
    WordTarget absent;absent.word=1;absent.count=0;absent.carriers=1;
    absent.occurrences={{0,1}};absent.copy_histogram={{1,1}};
    index.targets={control,absent};
    for(int threads:{1,2,4}){
        const auto prepared=build_model(index,RawEvidence{},2,1,0,threads);
        assert(prepared.target_ids==std::vector<size_t>{1}&&prepared.zero_dosage_targets==1);
        const auto&m=prepared.model;assert(m.targets[0].absent_reward_units>m.targets[0].reward_units);
        const auto exact=solve_presence(m,0,threads),beam=solve_presence(m,1,threads,{},true);
        assert(exact.path==std::vector<int>{3}&&exact.score==m.score({3}));
        assert(beam.upper>=exact.score&&beam.score<=exact.score);
        auto folded=fold_constant_presence(m);const auto r=solve_presence(folded.model,0,threads);
        assert(r.score+folded.offset==exact.score);
    }
    std::cout<<"PASS complete count likelihood, zero-mass normalization, synthetic fit recovery, exact likelihood replay, shallow/dropout uncertainty, degenerate counts, raw model inclusion, joint optimum and pruning bounds\n";
}
