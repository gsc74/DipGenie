#include "GlobalPhaseMove.hpp"
#include <cassert>
#include <random>
#include <iostream>
using namespace hic;
using namespace hic::global_phase;
static DensePairResult scalar(const std::vector<std::vector<Score>>&c,int H,int budget,bool free_swaps=false){
    const size_t B=c.size(),D=size_t(H*H);const int R=int(std::min(size_t(budget),2*(B-1)));
    auto charge=[&](int a,int b){const int direct=(a/H!=b/H)+(a%H!=b%H);
        return free_swaps?std::min(direct,(a/H!=b%H)+(a%H!=b/H)):direct;};
    const size_t stride=size_t(R+1)*D;std::vector<int>trace(B*stride,-1);
    std::vector<Score> p(stride,impossible),n(stride,impossible);std::copy(c[0].begin(),c[0].end(),p.begin());
    for(size_t b=1;b<B;++b){
        std::fill(n.begin(),n.end(),impossible);
        for(int r=0;r<=R;++r)for(size_t z=0;z<D;++z){
            if(c[b][z]==impossible)continue;
            Score best=impossible;int prior=-1;
            for(size_t a=0;a<D;++a){const int q=r-charge(int(a),int(z));if(q<0)continue;
                if(p[size_t(q)*D+a]>best){best=p[size_t(q)*D+a];prior=int(a);}}
            if(prior>=0){n[size_t(r)*D+z]=best+c[b][z];trace[b*stride+size_t(r)*D+z]=prior;}
        }
        p.swap(n);
    }
    DensePairResult out;int r=-1,z=-1;
    for(int q=0;q<=R;++q)for(size_t a=0;a<D;++a)if(p[size_t(q)*D+a]>out.score){out.score=p[size_t(q)*D+a];r=q;z=int(a);}
    if(z<0)return out;
    out.pairs.resize(B);
    for(size_t b=B;b--;){out.pairs[b]=z;if(!b)break;const int a=trace[b*stride+size_t(r)*D+size_t(z)];assert(a>=0);r-=charge(a,z);z=a;}
    assert(r==0);return out;
}
int main(){
    std::mt19937 rng(214051);size_t cases=0;
    {
        // Review regression: orientation swaps must cost zero in the free
        // exchange domain, including the complete-domain repair dispatch.
        PresenceModel m;m.H=2;m.blocks=2;m.budget=0;m.free_swaps=true;
        std::vector<Block> blocks{{0,1,0,{0,1,2,3}},{1,2,0,{0,1,2,3}}};
        const std::vector<std::vector<Score>> costs{{0,10,0,0},{0,0,10,0}};
        auto actual=repair(m,blocks,costs,1);
        assert(actual.score==20&&actual.choices==std::vector<int>({1,2}));
    }
    for(int i=0;i<360;++i){
        const int H=1+rng()%5,B=1+rng()%(i<180?7:23),R=rng()%15;
        std::vector<std::vector<Score>>c(B,std::vector<Score>(H*H));
        for(auto&row:c)for(auto&x:row)x=rng()%13==0?impossible:Score(int(rng()%15)-7);
        auto a=scalar(c,H,R),b=dense_pair_repair(c,H,R,1),p=dense_pair_repair(c,H,R,4);
        auto wide=dense_detail::solve<uint32_t>(c,H,R,2,0);
        assert(a.score==b.score&&a.pairs==b.pairs&&a.pairs==p.pairs&&a.score==p.score&&a.pairs==wide.pairs&&a.score==wide.score);
        // Checkpointed traceback must reproduce the full-table result exactly
        // for every interval, including one block and the whole sequence.
        for(size_t k:{size_t(1),size_t(2),size_t(3),size_t(5),size_t(B),size_t(B+4)}){
            auto q=dense_pair_repair(c,H,R,3,k);assert(q.score==a.score&&q.pairs==a.pairs);
        }
        ++cases;
    }
    {
        std::vector<std::vector<Score>>c(45,std::vector<Score>(9));
        for(auto&row:c)for(auto&x:row)x=Score(int(rng()%9)-4)*(Score(1)<<80);
        auto a=scalar(c,3,88),b=dense_pair_repair(c,3,88,4);assert(a.score==b.score&&a.pairs==b.pairs);++cases;
    }
    {
        // Force the 32-bit predecessor representation without a quadratic
        // scalar oracle. The unique high-reward pair changes BOTH donors.
        const int H=256,D=H*H;std::vector<std::vector<Score>>c(2,std::vector<Score>(D));
        c[0][0]=10;c[1][D-1]=20;
        auto a=dense_pair_repair(c,H,2,4);assert(a.score==30&&a.pairs==std::vector<int>({0,D-1}));
        // B=2, R=2: interval s=min(B,ceil(sqrt(2*16/4)))=2, one score checkpoint.
        assert(a.trace_bytes==3*size_t(D)*sizeof(Score)+2*3*size_t(D)*sizeof(uint32_t));
        auto b=dense_pair_repair(c,H,1,4);assert(b.score==20);++cases;
    }
    {
        std::vector<std::vector<Score>>c(3,std::vector<Score>(4,impossible));
        auto a=dense_pair_repair(c,2,4,2);assert(a.score==impossible&&a.pairs.empty());
        auto b=dense_pair_repair(c,2,4,2,0,true);assert(b.score==impossible&&b.pairs.empty());
    }
    size_t free_cases=0;
    for(int i=0;i<480;++i) {
        const int H=1+rng()%5,B=1+rng()%(i<240?7:23),R=i%19;
        std::vector<std::vector<Score>> c(B,std::vector<Score>(H*H));
        for(auto&row:c)for(auto&x:row)x=rng()%11==0?impossible:Score(int(rng()%11)-5);
        if(i%7==0)for(auto&row:c)for(auto&x:row)if(x!=impossible)x*=Score(1)<<80;
        const auto exact=scalar(c,H,R,true),first=dense_pair_repair(c,H,R,1,0,true);
        assert(first.score==exact.score);
        for(int threads:{1,2,4,8})for(size_t interval:{size_t(1),size_t(3),size_t(B+4)}) {
            const auto actual=dense_pair_repair(c,H,R,threads,interval,true);
            assert(actual.score==exact.score&&actual.pairs==first.pairs);
        }
        const auto wide=dense_detail::solve<uint32_t>(c,H,R,2,0,true);
        assert(wide.score==first.score&&wide.pairs==first.pairs);
        // Dispatch must use the model charge, not the shape of its domain.
        PresenceModel model;model.H=H;model.blocks=B;model.budget=R;model.free_swaps=true;
        std::vector<Block> blocks(B);
        for(int b=0;b<B;++b){blocks[b].begin=b;blocks[b].end=b+1;
            for(int z=0;z<H*H;++z)blocks[b].choices.push_back(z);}
        const auto dispatched=repair(model,blocks,c,2);
        assert(dispatched.score==first.score&&dispatched.choices==first.pairs);
        ++free_cases;
    }
    {
        const int H=256,D=H*H;
        std::vector<std::vector<Score>> c(2,std::vector<Score>(D));
        c[0][1]=10;c[1][H]=20;
        const auto x=dense_pair_repair(c,H,0,4,0,true);
        assert(x.score==30&&x.pairs==std::vector<int>({1,H}));
        ++free_cases;
    }
    for(int bad=0;bad<4;++bad){bool rejected=false;try{
        if(bad==0)dense_pair_repair({},2,2,1);
        if(bad==1)dense_pair_repair({{1,2,3}},2,2,1);
        if(bad==2)dense_pair_repair({{1,2,3,4}},2,-1,1);
        if(bad==3)dense_pair_repair({{1,2,3,4}},2,2,0);
    }catch(const std::exception&){rejected=true;}assert(rejected);}
    std::cout<<"PASS full donor-pair DP: "<<cases<<" ordered and "<<free_cases
             <<" min-swapped oracle cases; deterministic 1/2/4/8 threads, checkpoint intervals, negative/wide scores, unreachable states, 16/32-bit traces\n";
}
