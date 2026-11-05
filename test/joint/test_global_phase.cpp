#include <cmath>
#include "GlobalPhaseMove.hpp"
#include <cassert>
#include <random>
#include <iostream>
using namespace hic;
using namespace hic::global_phase;
int main() {
    std::mt19937 rng(113);
    size_t paths=0,improved=0;
    {DonorMask a,b;a.set(1);a.set(128);b.set(63);b.set(64);a|=b;
     assert(a.count()==4&&a.test(1)&&a.test(63)&&a.test(64)&&a.test(128));}
    for(int test=0;test<240;++test) {
        PresenceModel m;m.H=2+rng()%2;m.blocks=1+rng()%5;m.budget=rng()%7;
        std::vector<Block> b(m.blocks);std::vector<std::vector<Score>> costs(m.blocks);
        for(int i=0;i<m.blocks;++i) {
            b[i]={i,i+1,0,{0}};
            for(int j=0;j<3;++j){int z=int(rng()%m.D());if(std::find(b[i].choices.begin(),b[i].choices.end(),z)==b[i].choices.end())b[i].choices.push_back(z);}
            for(size_t j=0;j<b[i].choices.size();++j)costs[i].push_back(int(rng()%21)-10);
        }
        Score exact=impossible;std::vector<int> picked(m.blocks);
        std::function<void(int,int,Score)> enumerate=[&](int i,int r,Score s) {
            if(i==m.blocks){if(r<=m.budget)exact=std::max(exact,s);++paths;return;}
            for(size_t j=0;j<b[i].choices.size();++j){picked[i]=int(j);enumerate(i+1,r+(i?m.charge(b[i-1].choices[picked[i-1]],b[i].choices[j]):0),s+costs[i][j]);}
        };
        enumerate(0,0,0);
        const auto a=repair(m,b,costs,1),c=repair(m,b,costs,4);
        assert(a.score==exact&&c.score==exact&&a.choices==c.choices&&a.trace_bytes==c.trace_bytes);
    }
    // More than one resource tile: exercises actual concurrent repair tasks.
    {
        PresenceModel m;m.H=3;m.blocks=90;m.budget=178;
        std::vector<Block> blocks(m.blocks);std::vector<std::vector<Score>> cost(m.blocks);
        for(int b=0;b<m.blocks;++b){blocks[b]={b,b+1,0,{0,1,3,4,8}};
            for(int j=0;j<5;++j)cost[b].push_back(int(rng()%100)-50);}
        auto a=repair(m,blocks,cost,1),b=repair(m,blocks,cost,4);
        assert(a.score==b.score&&a.choices==b.choices&&a.trace_bytes==b.trace_bytes);
    }
    for(int test=0;test<100;++test) {
        int n=3+rng()%12;std::vector<Edge> es;
        for(int a=0;a<n;++a)for(int b=a+1;b<n;++b)if(rng()%3==0)es.push_back({a,b,int(rng()%21)-10});
        std::vector<int>s(n);for(auto&v:s)v=rng()%2?1:-1;
        auto output=intervals(es,s);assert(energy(es,output)>=energy(es,s));
        auto greedy=coordinate(es,s);assert(energy(es,greedy)>=energy(es,s));
        for(int a=0;a<n;++a)for(int b=a;b<n;++b) {
            auto trial=output;for(int i=a;i<=b;++i)trial[i]=-trial[i];
            assert(energy(es,trial)<=energy(es,output));
        }
    }
    // Repeated strings, neutral endpoints, cis/trans contacts and R=0.
    for(int test=0;test<160;++test) {
        PresenceModel m;m.H=3;m.blocks=6;m.budget=test%7;m.weights={2,3};
        for(int b=0;b<m.blocks;++b)for(int a=0;a<2;++a) {
            DonorMask mask;mask.set(size_t((b+a)%m.H));
            PresenceTarget t{bool(rng()%2),{{b,mask}},Score(1+rng()%4),Score(rng()%3)};
            if(test%3==0&&b<5)t.occurrences.push_back({b+1,mask});
            m.targets.push_back(t);
        }
        for(size_t a=0;a<m.targets.size();++a)for(size_t b=a+1;b<m.targets.size();++b)
            if(rng()%4==0)m.contacts.push_back({a,b,Score(1+rng()%8),rng()%2?1:-1,bool(rng()%2)});
        m.check();Result incumbent;incumbent.path.assign(m.blocks,1);incumbent.score=m.score(incumbent.path);incumbent.upper=m.universal_bound();
        Result parallel=incumbent;const Score initial=incumbent.score;
        const auto report=move(m,incumbent,2,1);const auto other=move(m,parallel,2,4);
        assert(incumbent.path==parallel.path&&incumbent.score==parallel.score&&report.repaired==other.repaired);
        assert(incumbent.score>=initial&&incumbent.score==m.score(incumbent.path)&&incumbent.score<=incumbent.upper);
        if(report.edges) {
            const size_t resources=1+std::min(size_t(m.budget),2*(report.blocks-1));
            // Checkpointed traceback: score checkpoints every s blocks plus one
            // interval of 16-bit predecessors (DensePairRepair.hpp).
            const size_t B=report.blocks,stride=resources*size_t(m.D());
            const size_t s=std::min(B,std::max<size_t>(1,size_t(std::ceil(std::sqrt(double(B)*double(sizeof(Score))/2.0)))));
            size_t checkpoints=1;for(size_t b=1;b<B;++b)if(b%s==0&&b+1<B)++checkpoints;
            assert(report.trace_bytes==checkpoints*stride*sizeof(Score)+s*stride*sizeof(uint16_t));
            assert(report.repaired_recombinations<=m.budget);
        }
        if(incumbent.score>initial)++improved;
        std::vector<int> exchanged=incumbent.path;for(auto&z:exchanged)z=z%m.H*m.H+z/m.H;
        assert(m.score(exchanged)==incumbent.score);
    }
    assert(improved>0);
    std::cout<<"global phase move tests PASS; exhaustive repair paths="<<paths<<" improving full-objective cases="<<improved<<"\n";
}
