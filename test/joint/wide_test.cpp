#include "PanelSelection.hpp"
#include <cassert>
#include <random>
int main(){
    using hic::Score;const Score big=Score(1)<<80;
    hic::Result certificate;certificate.path={0};certificate.score=3;certificate.upper=4;
    assert(certificate.ratio_text()=="1.333334");certificate.score=big*3;certificate.upper=big*4;assert(certificate.ratio_text()=="1.333334");
    hic::PackedScoreRows rows(4);std::vector<Score> input{0,1,big,big+3},decoded;
    rows.append(input);assert(rows.width(0)==16);rows.decode(0,decoded);assert(input==decoded);
    rows.append({big,big+4,big+8,big+12});assert(rows.width(1)==1);rows.decode(1,decoded);assert(decoded==std::vector<Score>({big,big+4,big+8,big+12}));
    // Cached resource rows must reproduce every scalar query, including
    // final blocks, clipping, tile reconstruction and wide integer scores.
    for(int H:{1,2,7})for(int blocks:{1,7})for(int budget:{0,1,20}) {
        hic::PackedScoreRows local(size_t(H*H));
        std::vector<Score> other(size_t(blocks+1),0);
        for(int b=0;b<blocks;++b) {
            std::vector<Score> row(size_t(H*H));
            for(int z=0;z<H*H;++z)row[z]=big*(1+(b+z*3)%11);
            local.append(row);other[b]=big*(blocks-b);
        }
        hic::PackedCompletion completion(H,blocks,budget,std::move(local),other,2,4);
        for(int b=0;b<blocks;++b) {
            completion.prepare(b);
            for(int remaining=0;remaining<=budget;++remaining)for(Score cap:{Score(0),big,big*1000}) {
                std::vector<Score> row(size_t(H*H));
                completion.remaining_row(b,remaining,cap,row.data());
                for(int z=0;z<H*H;++z)assert(row[z]==completion.remaining(b,z,remaining,cap));
            }
        }
    }
    hic::PresenceModel m;m.H=3;m.blocks=3;m.budget=2;m.weights={3,2};
    m.targets={{false,{{0,1}},5,2},{true,{{1,2}},6,1},{false,{{2,3}},7,0}};
    m.contacts={{0,2,11,1,false},{1,2,3,-1,true}};
    auto a=hic::solve_presence_bidirectional(m,3,1);m.weights.genotype*=big;m.weights.contact*=big;
    auto b=hic::solve_presence_bidirectional(m,3,4);
    assert(a.best.path==b.best.path&&b.best.score==a.best.score*big&&b.best.upper==a.best.upper*big);
    assert(a.best.ratio_text()==b.best.ratio_text());
    bool caught=false;try{(void)hic::score_product(hic::score_limit,2);}catch(const std::overflow_error&){caught=true;}assert(caught);
    hic::ReadSketch sketch(31);std::map<uint64_t,uint64_t> full;std::mt19937 rng(42);
    for(size_t i=0;i<10000;++i){auto w=uint64_t(rng()%1000);sketch.add(w);++full[w];}
    std::vector<std::pair<uint64_t,uint64_t>> order;for(auto[w,n]:full){(void)n;order.emplace_back(hic::sketch_hash(w),w);}std::sort(order.begin(),order.end());
    assert(sketch.counts.size()==31);for(size_t i=0;i<31;++i)assert(sketch.counts.at(order[i].second)==full.at(order[i].second));
    std::cout<<"PASS 128-bit packed rows, exact large-score DP/parallel/certificate invariance, overflow guard and streaming sketch counts\n";
}
