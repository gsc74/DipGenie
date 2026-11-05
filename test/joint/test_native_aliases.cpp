#include "NativeLocalizer.hpp"
#include <cassert>
#include <iostream>
#include <random>
using namespace hic;

// Frozen pre-optimization oracle. Tree keys own their complete footprints.
static std::vector<PlacementCandidate> reference(std::vector<PlacementCandidate> candidates) {
    using Key=std::tuple<uint32_t,std::vector<PlacementSlice>,uint32_t,uint32_t,int64_t>;
    std::map<Key,size_t> aliases;std::vector<PlacementCandidate> distinct;
    for(auto&c:candidates) {
        if(c.graph_context.empty()){distinct.push_back(std::move(c));continue;}
        Key key{c.component,c.graph_context,c.context_clip_left,c.context_clip_right,c.context_jitter};
        const auto [at,inserted]=aliases.emplace(std::move(key),distinct.size());
        if(inserted){distinct.push_back(std::move(c));continue;}
        auto&kept=distinct[at->second];
        if(c.support>kept.support)kept=std::move(c);
        else if(c.support==kept.support&&c.anchor.status==AnchorStatus::Known&&
                (kept.anchor.status!=AnchorStatus::Known||
                 std::tie(c.anchor.start,c.anchor.word)<std::tie(kept.anchor.start,kept.anchor.word)))
            kept.anchor=c.anchor;
    }
    return distinct;
}
static void same(const std::vector<PlacementCandidate>&a,const std::vector<PlacementCandidate>&b) {
    assert(a.size()==b.size());
    for(size_t i=0;i<a.size();++i) {
        const auto&x=a[i];const auto&y=b[i];
        assert(same_placement_alias(x,y));
        assert(std::tie(x.space,x.support,x.reverse,x.diagonal_low,x.diagonal_high)==
               std::tie(y.space,y.support,y.reverse,y.diagonal_low,y.diagonal_high));
        assert(std::tie(x.anchor.status,x.anchor.contig,x.anchor.start,x.anchor.word)==
               std::tie(y.anchor.status,y.anchor.contig,y.anchor.start,y.anchor.word));
        assert(x.intervals==y.intervals&&x.seeds.size()==y.seeds.size());
        for(size_t k=0;k<x.seeds.size();++k) {
            const auto&p=x.seeds[k];const auto&q=y.seeds[k];
            assert(std::tie(p.word,p.position,p.query,p.context,p.reverse)==
                   std::tie(q.word,q.position,q.query,q.context,q.reverse));
        }
    }
}
struct Collisions {size_t operator()(const PlacementCandidate&)const{return 0;}};
int main() {
    std::mt19937_64 rng(20261031);
    for(size_t trial=0;trial<1000;++trial) {
        std::vector<PlacementCandidate> input;
        for(size_t j=0,n=trial%257;j<n;++j) {
            PlacementCandidate c;
            if(j&&rng()%3)c=input[rng()%j];
            else {
                c.component=rng()%4;c.context_clip_left=rng()%3;c.context_clip_right=rng()%3;
                c.context_jitter=int64_t(rng()%9)-4;
                for(size_t k=0,n=rng()%12;k<n;++k)
                    c.graph_context.push_back({rng()%100,rng()%4,rng()%50,50+rng()%50,bool(rng()%2)});
            }
            c.support=rng()%7;c.space=rng()%101;c.reverse=rng()%2;
            c.diagonal_low=int64_t(rng()%1000)-500;c.diagonal_high=c.diagonal_low+3;
            c.anchor={rng()%2?AnchorStatus::Known:AnchorStatus::Unprojected,uint32_t(rng()%3),rng()%5,rng()%7};
            c.intervals={{uint32_t(j),uint32_t(j+10)}};
            c.seeds={{rng(),rng(),uint32_t(j),uint32_t(trial),bool(rng()%2)}};
            input.push_back(std::move(c));
        }
        const auto expected=reference(input);
        same(expected,merge_placement_aliases(input));
        // Exact equality, never hash equality, determines physical identity.
        same(expected,merge_placement_aliases(input,Collisions{}));
    }
    PlacementCandidate a;a.graph_context={{5,0,0,100,false}};a.support=250;
    a.anchor={AnchorStatus::Known,0,50,7};
    auto b=a;b.anchor.start=40;b.space=99;b.seeds={{1,2,3,4,false}};
    auto c=b;c.graph_context[0].copy=1;
    auto d=a;d.graph_context.clear();
    const auto merged=merge_placement_aliases({a,b,c,d,d});
    assert(merged.size()==4&&merged[0].anchor.start==40&&merged[0].space==a.space&&merged[0].seeds.empty());
    same(reference({a,b,c,d,d}),merged);
    std::cout<<"PASS native alias merging: 1000 differential cases, forced hash collisions, exact occurrence keys, stable support/anchor ties and unkeyed contexts\n";
}
