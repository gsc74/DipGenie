#pragma once
#include "NativeLocalizer.hpp"
#include <numeric>

namespace hic {
inline constexpr const char* native_graph_context_version="graph-footprint-local-repeat-v1";
// Continuous panel contexts with exact graph-coordinate maps. Coordinates on
// nonreference paths are NEVER treated as reference coordinates: an anchor is
// inherited only by an exact k-mer with the same physical graph footprint as
// an indexed occurrence on the designated reference.
//
// Unique-node copies share physical identity across donor paths. If a node is
// revisited, occurrence correspondence between different cyclic routes is not
// implied by its DNA or node ID. Such visits are conservatively distinguished
// by (canonical local route between nonrepeated flanks, visit rank). Identical
// local contexts share identities even when donors differ far away. Different
// repeat routes retain ambiguity rather than invent copy correspondence.
inline std::vector<PlacementSequence> native_graph_contexts(const NativeGraph&g,
        size_t reference,uint32_t component=0,uint32_t reference_contig=0,uint32_t space_base=0) {
    if(reference>=g.paths.size()||g.paths.size()>=UINT32_MAX||g.labels.empty())
        throw std::invalid_argument("Native graph reference/path dimensions");
    const size_t H=g.paths.size(),V=g.labels.size();
    std::vector<uint32_t> maximum(V),counts(V);
    std::vector<bool> reversed(H);
    for(size_t h=0;h<H;++h) {
        const auto&walk=g.paths[h].walk;
        if(walk.empty())throw std::invalid_argument("Empty native graph path");
        std::vector<uint32_t> touched;
        for(size_t i=0;i<walk.size();++i) {
            const uint32_t node=walk[i]/2;
            if(node>=V||g.labels[node].empty())throw std::invalid_argument("Native graph path label");
            if(i&&!g.edges.count(NativeGraph::edge(walk[i-1],walk[i])))
                throw std::invalid_argument("Native graph context crosses a missing edge");
            if(!counts[node])touched.push_back(node);
            if(counts[node]==UINT32_MAX)throw std::overflow_error("Native graph visit count");
            maximum[node]=std::max(maximum[node],++counts[node]);
        }
        for(auto node:touched)counts[node]=0;
        // Canonical route orientation is independent of donor name and order.
        for(size_t i=0;i<walk.size();++i) {
            const uint32_t a=walk[i],b=walk[walk.size()-i-1]^1;
            if(a!=b){reversed[h]=b<a;break;}
        }
    }
    auto at=[&](size_t h,size_t i) {
        const auto&w=g.paths[h].walk;return reversed[h]?(w[w.size()-i-1]^1):w[i];
    };
    auto less=[&](size_t a,size_t b) {
        const auto n=std::min(g.paths[a].walk.size(),g.paths[b].walk.size());
        for(size_t i=0;i<n;++i)if(at(a,i)!=at(b,i))return at(a,i)<at(b,i);
        return g.paths[a].walk.size()<g.paths[b].walk.size();
    };
    std::vector<size_t> order(H);std::iota(order.begin(),order.end(),0);
    std::sort(order.begin(),order.end(),[&](size_t a,size_t b){return less(a,b)||(!less(b,a)&&a<b);});
    std::vector<PlacementSequence> result;
    std::map<std::vector<uint32_t>,uint32_t> repeat_contexts;
    for(size_t first=0;first<H;) {
        size_t last=first+1,h=order[first];
        while(last<H&&!less(h,order[last])&&!less(order[last],h))++last;
        bool is_reference=false;
        for(size_t i=first;i<last;++i)if(order[i]==reference){h=reference;is_reference=true;break;}
        if(result.size()>UINT32_MAX-space_base)throw std::overflow_error("Native graph context space");
        PlacementSequence context;
        context.name=g.paths[h].name;context.sequence=g.spell(h);context.component=component;
        context.space=space_base+uint32_t(result.size());
        if(is_reference)context.reference_contig=reference_contig;
        const auto&walk=g.paths[h].walk;
        std::vector<uint64_t> copies(walk.size());
        for(size_t a=0;a<walk.size();) {
            if(maximum[walk[a]/2]<=1){++a;continue;}
            size_t b=a+1;while(b<walk.size()&&maximum[walk[b]/2]>1)++b;
            // Include the immediately surrounding nonrepeated physical nodes;
            // explicit sentinels distinguish an unbounded end from a flank.
            std::vector<uint32_t> route;
            route.push_back(a?walk[a-1]:UINT32_MAX);
            route.insert(route.end(),walk.begin()+a,walk.begin()+b);
            route.push_back(b<walk.size()?walk[b]:UINT32_MAX);
            auto opposite=route;std::reverse(opposite.begin(),opposite.end());
            for(auto&handle:opposite)if(handle!=UINT32_MAX)handle^=1;
            const bool reverse=opposite<route;
            if(reverse)route=std::move(opposite);
            if(repeat_contexts.size()>=UINT32_MAX-1)throw std::overflow_error("Native repeat contexts");
            const auto entry=repeat_contexts.emplace(std::move(route),uint32_t(repeat_contexts.size()+1)).first;
            for(size_t i=a;i<b;++i)++counts[walk[i]/2];
            std::unordered_map<uint32_t,uint32_t> seen;
            for(size_t i=a;i<b;++i) {
                const auto node=walk[i]/2,rank=seen[node]++;
                const auto ordinal=reverse?counts[node]-1-rank:rank;
                copies[i]=(uint64_t(entry->second)<<32)|(uint64_t(ordinal)+1);
            }
            for(size_t i=a;i<b;++i)counts[walk[i]/2]=0;
            a=b;
        }
        uint64_t position=0;
        for(size_t i=0;i<walk.size();++i) {
            const auto handle=walk[i];const uint32_t node=handle/2;
            if(g.labels[node].size()>UINT64_MAX-position)throw std::overflow_error("Native graph context length");
            position+=g.labels[node].size();
            context.runs.push_back({position,node,copies[i],0,bool(handle&1)});
        }
        result.push_back(std::move(context));first=last;
    }
    return result;
}
} // namespace hic
