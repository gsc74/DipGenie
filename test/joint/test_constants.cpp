#include "ConstantPresence.hpp"
#include <cassert>
#include <iostream>
#include <random>
using namespace hic;
int main() {
    std::mt19937 rng(613);
    for(int trial=0;trial<100;++trial) {
        PresenceModel m;m.H=2+trial%2;m.blocks=4;m.budget=trial%7;m.weights={3,5};
        for(int x=0;x<6;++x) {
            PresenceTarget t;t.hom=rng()%2;t.reward_units=1+rng()%100;t.other_reward_units=rng()%100;
            for(int b=0;b<m.blocks;++b)if(rng()%2)t.occurrences.push_back({b,1+rng()%((1u<<m.H)-1)});
            if(t.occurrences.empty())t.occurrences.push_back({0,DonorMask::full(m.H)});
            m.targets.push_back(t);
        }
        for(int x=0;x<5;++x)m.contacts.push_back({size_t(x),size_t(x+1),1+rng()%30,rng()%2?1:-1,bool(rng()%2)});
        const auto folded=fold_constant_presence(m);
        Score best=impossible;
        std::vector<int> path(m.blocks);
        const int domain=m.D(),paths=domain*domain*domain*domain;
        for(int code=0;code<paths;++code) {
            int c=code;for(auto&z:path){z=c%domain;c/=domain;}
            const auto old=m.score(path),now=folded.model.score(path);
            assert((old==impossible&&now==impossible)||old==folded.offset+now);
            best=std::max(best,old);
        }
        for(int threads:{1,4})for(size_t beam:{size_t(0),size_t(2),size_t(8)}) {
            auto result=solve_presence(folded.model,beam,threads,{},true);
            assert(!result.path.empty());
            assert(m.score(result.path)==result.score+folded.offset);
            assert(best<=result.upper+folded.offset);
            if(!beam)assert(best==result.score+folded.offset);
            if(threads>1) {
                const auto serial=solve_presence(folded.model,beam,1,{},true);
                assert(result.path==serial.path&&result.score==serial.score&&result.upper==serial.upper);
            }
            if(beam)assert(result.peak<=beam);
        }
    }
    PresenceModel wide;wide.H=65;wide.blocks=2;wide.budget=2;wide.weights={Score(1)<<65,1};
    wide.targets={{true,{{0,DonorMask::full(65)}},Score(1)<<30,0}};
    const auto fold=fold_constant_presence(wide);
    assert(fold.fixed_targets==1&&fold.model.targets.empty());
    assert(wide.score({0,0})==fold.offset);
    auto r=solve_presence(fold.model,4);assert(wide.score(r.path)==r.score+fold.offset);
    std::cout<<"constant-folding: exhaustive score, exact optimum, pruned bounds, threads and wide masks PASS\n";
}
