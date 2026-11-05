#include "ConstantPresence.hpp"
#include "GlobalPhaseMove.hpp"
#include "FreeOrientation.hpp"
#include <cassert>
#include <iostream>
#include <random>
using namespace hic;
// Exhaustive checks of the unordered-pair recombination charge (free haplotype
// exchange) against the ordered v1 charge on random small models.
static Score brute(const PresenceModel&m) {
    Score best=impossible;std::vector<int> path(m.blocks);const int D=m.D();
    long total=1;for(int b=0;b<m.blocks;++b)total*=D;
    for(long code=0;code<total;++code){long c=code;for(auto&z:path){z=int(c%D);c/=D;}best=std::max(best,m.score(path));}
    return best;
}
int main() {
    std::mt19937 rng(20261027);
    // Charge identities.
    for(int H:{1,2,3,5}) {
        PresenceModel m;m.H=H;m.free_swaps=true;
        for(int a=0;a<m.D();++a)for(int b=0;b<m.D();++b) {
            const int s=(b%H)*H+b/H;
            assert(m.charge(a,b)==std::min(m.ordered_charge(a,b),m.ordered_charge(a,s)));
            assert(m.charge(a,b)==m.charge(b,a)&&m.charge(a,b)<=m.ordered_charge(a,b));
            assert(m.charge(a,s)==m.charge(a,b));
        }
    }
    // Linked presence: some haplotype must be able to emit the linked pair.
    for(unsigned a=0;a<4;++a)for(unsigned b=0;b<4;++b) {
        const bool cis=(a&b)!=0,trans=((a&1)&&(b&2))||((a&2)&&(b&1));
        assert(phase_factor(a,b,1,false,false,true)==(cis?2:0));
        assert(phase_factor(a,b,-1,false,false,true)==(trans?2:0));
        assert(phase_factor(a,b,1,false,true,true)==phase_factor(a,b,1,false,false,true));
    }
    assert(phase_factor(3,3,1,false,false,true)==2&&phase_factor(1,1,1,false,false,true)==2);  // no het bonus
    assert(phase_factor(0,1,1,false,false,true)==0&&phase_factor(1,2,1,false,false,true)==0);  // missing/trans
    size_t strictly_better=0;
    for(int trial=0;trial<270;++trial) {
        PresenceModel m;m.H=2+trial%3;m.blocks=m.H==4?4:5;m.budget=trial%5;m.weights={3,5};
        for(int x=0;x<7;++x) {
            PresenceTarget t;t.hom=rng()%2;t.reward_units=1+rng()%100;t.other_reward_units=rng()%100;t.absent_reward_units=rng()%60;
            for(int b=0;b<m.blocks;++b)if(rng()%2)t.occurrences.push_back({b,1+rng()%((1u<<m.H)-1)});
            if(t.occurrences.empty())t.occurrences.push_back({int(rng()%m.blocks),1+rng()%((1u<<m.H)-1)});
            m.targets.push_back(t);
        }
        for(int x=0;x<9;++x){size_t a=rng()%7,b=rng()%7;if(a==b)b=(a+1)%7;m.contacts.push_back({a,b,1+rng()%40,rng()%2?1:-1,bool(rng()%2)});}
        m.check();
        // symmetric, neutral and linked-presence factors, both charges
        m.neutral_phase=(trial/3)%3==1;m.presence_phase=(trial/3)%3==2;m.ordered_repair=(trial/9)%2;  // both repair charges
        PresenceModel free=m;free.free_swaps=true;
        const Score ordered_best=brute(m),free_best=brute(free);
        assert(free_best>=ordered_best);strictly_better+=free_best>ordered_best;
        for(const auto*model:{&m,&free}) {
            const Score best=model==&m?ordered_best:free_best;
            const auto folded=fold_constant_presence(*model);
            assert(folded.model.free_swaps==model->free_swaps&&folded.model.neutral_phase==model->neutral_phase&&folded.model.presence_phase==model->presence_phase);
            {   // folding is exact for every complete path
                std::vector<int> path(model->blocks);const int D=model->D();long total=1;for(int b=0;b<model->blocks;++b)total*=D;
                for(long code=0;code<total;code+=7){long c=code;for(auto&z:path){z=int(c%D);c/=D;}
                    const auto a=model->score(path),b=folded.model.score(path);
                    assert((a==impossible&&b==impossible)||a==b+folded.offset);}
            }
            for(int threads:{1,3})for(size_t beam:{size_t(0),size_t(1),size_t(3),size_t(16)}) {
                const auto r=solve_presence(folded.model,beam,threads,{},beam!=0);
                assert(!r.path.empty()&&model->score(r.path)==r.score+folded.offset);
                assert(best<=r.upper+folded.offset); // completion bound valid under this charge
                if(!beam)assert(r.score+folded.offset==best);
            }
            const auto both=solve_presence_bidirectional(folded.model,0,2);
            assert(both.best.score+folded.offset==best);
            // Global move keeps feasibility under the model's own charge.
            auto incumbent=solve_presence(*model,2,2,{},true);
            const auto before=incumbent.score;
            (void)global_phase::move(*model,incumbent,std::numeric_limits<size_t>::max(),2);
            assert(model->score(incumbent.path)==incumbent.score&&incumbent.score>=before&&incumbent.score<=best);
            if(model->free_swaps) {
                const auto prior=incumbent.score;
                const auto r=free_orientation::improve(*model,incumbent);
                assert(model->score(incumbent.path)==incumbent.score&&incumbent.score>=prior&&incumbent.score<=best);
                assert(r.after==incumbent.score&&r.before==prior);
            }
        }
    }
    assert(strictly_better>0);
    // Phase-only case: one heterozygous marker per segment, read links demand
    // an orientation change between segments 1 and 2. With R=0 the ordered
    // charge must keep one orientation; free exchange satisfies every link.
    {
        PresenceModel m;m.H=2;m.blocks=4;m.budget=0;m.weights={1,1};
        for(int b=0;b<4;++b){PresenceTarget t;t.hom=false;t.reward_units=10;t.occurrences={{b,DonorMask(1)}};m.targets.push_back(t);}
        m.contacts={{0,1,100,1,false},{1,2,100,-1,false},{2,3,100,1,false}};
        PresenceModel f=m;f.free_swaps=true;
        const auto a=solve_presence(m,0),b=solve_presence(f,0);
        assert(b.score>a.score&&f.score(b.path)==b.score);
        assert(f.charge(b.path[1],b.path[2])==0&&m.ordered_charge(b.path[1],b.path[2])==2);
    }
    // Neutral factor values.
    for(int r:{-1,1})for(unsigned char a=0;a<4;++a)for(unsigned char b=0;b<4;++b){
        auto sw=[](unsigned char x){return (unsigned char)(((x&1)<<1)|((x&2)>>1));};
        const auto v=phase_factor(a,b,r,false,true);
        const int x=phase_orientation(a),y=phase_orientation(b);
        assert(v==(r*x*y<0?0:2)&&v==phase_factor(sw(a),sw(b),r,false,true));
        if(!x||!y)assert(v==2);
    }
    // Orientation search repairs a planted wrong-orientation interval exactly.
    for(bool neutral:{false,true}) {
        PresenceModel m;m.H=3;m.blocks=40;m.budget=0;m.weights={1,1};m.free_swaps=true;m.neutral_phase=neutral;
        for(int b=0;b<40;++b){PresenceTarget t;t.hom=false;t.reward_units=10;t.other_reward_units=0;t.occurrences={{b,DonorMask(1)}};m.targets.push_back(t);}
        for(int b=0;b+1<40;++b)m.contacts.push_back({size_t(b),size_t(b+1),50,1,false});
        for(int b=0;b+5<40;b+=3)m.contacts.push_back({size_t(b),size_t(b+5),20,1,false});
        Result r;r.path.assign(40,0*3+1);for(int b=12;b<27;++b)r.path[b]=1*3+0;  // donors (0,1), flipped in [12,27)
        r.score=m.score(r.path);r.upper=m.universal_bound();
        const auto before=r.score;
        const auto rep=free_orientation::improve(m,r);
        const auto best=solve_presence(m,0);
        assert(rep.accepted>=1&&r.score>before&&r.score==best.score);
        int o=r.path[0];for(int z:r.path)assert(z==o);
    }
    // Continuity tie-break: with no phase evidence every exchange is free and
    // score-neutral; the pass must remove all of them without changing score.
    {
        PresenceModel m;m.H=4;m.blocks=30;m.budget=6;m.weights={1,1};m.free_swaps=true;
        for(int b=0;b<30;++b){PresenceTarget t;t.hom=false;t.reward_units=10;t.other_reward_units=3;t.occurrences={{b,DonorMask(1u<<(b%4))}};m.targets.push_back(t);}
        Result r;r.path.resize(30);
        for(int b=0;b<30;++b){int a=b<10?0:1,c=b<20?2:3;r.path[b]=(b%3==1)?c*4+a:a*4+c;}   // exchanges sprinkled
        r.score=m.score(r.path);r.upper=m.universal_bound();
        const auto before=r.score;const int e0=free_orientation::exchanges(m,r.path);
        const auto rep=free_orientation::improve(m,r);
        assert(e0>0&&rep.exchanges_before==e0&&rep.exchanges_after==0&&r.score>=before&&m.score(r.path)==r.score);
        for(size_t b=1;b<r.path.size();++b)assert(m.ordered_charge(r.path[b-1],r.path[b])==m.charge(r.path[b-1],r.path[b]));
    }
    std::cout<<"PASS unordered-pair charge: identities, exhaustive exact/beam optima, completion bounds, folding, bidirectional, global move, phase-only exchange, neutral factor, exact folding, orientation search, continuity tie-break; strictly better in "<<strictly_better<<" trials\n";
}
