#include "NativeGraphPlacement.hpp"
#include <cassert>
#include <iostream>
#include <random>

using namespace hic;
static std::string dna(size_t n,uint64_t seed) {
    std::mt19937_64 rng(seed);std::string s(n,'A');for(auto&c:s)c="ACGT"[rng()%4];return s;
}
static uint32_t node(NativeGraph&g,const std::string&s,bool reverse=false) {
    const auto id=uint32_t(g.labels.size());g.names.push_back(std::to_string(id));
    g.labels.push_back(reverse?reverse_complement(s):s);return 2*id+unsigned(reverse);
}
static void path(NativeGraph&g,std::vector<uint32_t>w,std::string name="panel") {
    for(size_t i=1;i<w.size();++i){g.edges.insert(NativeGraph::edge(w[i-1],w[i]));g.edges.insert(NativeGraph::edge(w[i]^1,w[i-1]^1));}
    g.paths.push_back({name,std::move(w)});
}
static std::vector<uint32_t> reverse_walk(std::vector<uint32_t>w) {
    std::reverse(w.begin(),w.end());for(auto&v:w)v^=1;return w;
}
static void same(const EndPlacement&a,const EndPlacement&b) {
    assert(a.status==b.status&&a.component==b.component&&a.best_support==b.best_support&&a.competing_support==b.competing_support);
    assert(a.anchor.status==b.anchor.status);
    if(a.anchor.status==AnchorStatus::Known)assert(std::tie(a.anchor.contig,a.anchor.start,a.anchor.word)==std::tie(b.anchor.contig,b.anchor.start,b.anchor.word));
    assert(a.retained.size()==b.retained.size());
}
static uint64_t word(const std::string&s,unsigned k) {
    uint64_t key=UINT64_MAX;words(s,k,[&](uint64_t x,size_t,bool){key=x;});return key;
}
int main() {
    NativeLocalizer::Policy policy;policy.k=15;policy.w=7;
    NativeGraph tiny;const auto ref=dna(600,319);
    std::vector<uint32_t> walk;
    for(size_t begin=0;begin<ref.size();) {
        const auto n=std::min(size_t(1+begin%7),ref.size()-begin);
        walk.push_back(node(tiny,ref.substr(begin,n),begin%2));begin+=n;
    }
    path(tiny,walk,"CHM13");path(tiny,reverse_walk(walk),"reverse-copy");
    for(int i=0;i<50;++i)path(tiny,walk,"donor"+std::to_string(i));
    const auto contexts=native_graph_contexts(tiny,0,7,19);
    assert(contexts.size()==1&&contexts[0].runs.size()==walk.size());
    assert(contexts[0].sequence==ref&&contexts[0].reference_contig==19);
    PlacementSequence flat;flat.name="flat";flat.sequence=ref;flat.component=7;flat.reference_contig=19;
    NativeLocalizer linear({flat},policy,1);
    for(int threads:{1,2,4,8}) {
        NativeLocalizer native(contexts,policy,threads);
        for(size_t begin:{size_t(0),size_t(43),size_t(201),size_t(450)}) {
            const auto read=ref.substr(begin,150);
            const auto placed=native.locate(read);
            same(placed,linear.locate(read));same(placed,native.locate(reverse_complement(read)));
            assert(placed.status==PlacementStatus::Unique&&placed.best_support==150);
            assert(placed.anchor.word==word(ref.substr(placed.anchor.start,policy.k),policy.k));
        }
        // Neither a graph node boundary nor a reverse-oriented node resets support.
        auto first=native.locate(ref.substr(43,150)),last=native.locate(ref.substr(350,150));
        auto distance=selected_anchor_distance(first,last);
        assert(distance.status==AnchorDistanceStatus::Known);
        assert(distance.separation==last.anchor.start-first.anchor.start);
        assert(selected_anchor_distance(last,first).separation==distance.separation);
    }
    // Same physical prefix in different donor routes is ONE locus/context.
    // A novel allele is localized continuously with an exact shared anchor;
    // an interior novel read is valid but has no reference anchor.
    NativeGraph bubble;
    const auto a=node(bubble,dna(200,101)),b=node(bubble,dna(250,102),true),
               alt=node(bubble,dna(300,103)),c=node(bubble,dna(200,104));
    path(bubble,{a,b,c},"CHM13");path(bubble,{a,alt,c},"alternate");
    path(bubble,reverse_walk({a,alt,c}),"alternate-reverse");
    auto panel=native_graph_contexts(bubble,0);
    assert(panel.size()==2);
    for(int threads:{1,2,4,8}) {
        NativeLocalizer index(panel,policy,threads);
        const auto shared=bubble.spell_vertex(a).substr(20,150);
        auto placed=index.locate(shared);
        assert(placed.status==PlacementStatus::Unique&&placed.best_support==150&&placed.retained.size()==1);
        same(placed,index.locate(reverse_complement(shared)));
        minimizers(shared,policy.k,policy.w,[&](auto x,size_t){assert(index.physical_loci(x)==1);});
        const auto across=bubble.spell(1).substr(150,150);
        placed=index.locate(across);
        assert(placed.status==PlacementStatus::Unique&&placed.best_support==150&&placed.anchor.status==AnchorStatus::Known);
        assert(placed.anchor.start>=150&&placed.anchor.start+policy.k<=200);
        same(placed,index.locate(reverse_complement(across)));
        const auto novel=bubble.spell_vertex(alt).substr(50,150);
        placed=index.locate(novel);
        assert(placed.status==PlacementStatus::Unique&&placed.anchor.status==AnchorStatus::Unprojected);
        same(placed,index.locate(reverse_complement(novel)));
        assert(selected_anchor_distance(placed,placed).status==AnchorDistanceStatus::Unknown);
        // Indels must not be treated as affine reference projection. The
        // shared suffix is at 450 in the reference, 500 in the alternate.
        const auto suffix=bubble.spell_vertex(c).substr(25,150);
        placed=index.locate(suffix);
        assert(placed.status==PlacementStatus::Unique&&placed.anchor.status==AnchorStatus::Known);
        assert(placed.anchor.start>=475&&placed.anchor.start<625);
        same(placed,index.locate(reverse_complement(suffix)));
    }
    // Equal DNA on a distinct physical vertex is a real repeat, not an alias.
    NativeGraph repeated=bubble;
    const auto duplicate=node(repeated,bubble.spell_vertex(a));
    path(repeated,{duplicate},"true-repeat");
    NativeLocalizer repeats(native_graph_contexts(repeated,0),policy,2);
    const auto shared=bubble.spell_vertex(a).substr(20,150);
    auto ambiguous=repeats.locate(shared);
    assert(ambiguous.status==PlacementStatus::Ambiguous&&ambiguous.component_unique());
    minimizers(shared,policy.k,policy.w,[&](auto x,size_t){assert(repeats.physical_loci(x)==2);});
    auto capped=policy;capped.maximum_loci=1;
    NativeLocalizer filtered(native_graph_contexts(repeated,0),capped,2);
    assert(filtered.locate(shared).status==PlacementStatus::NoSeeds);
    // Cyclic paths: identical donor copies are collapsed, but the two visits
    // to a repeated vertex are retained, including in reverse orientation.
    NativeGraph cycle;
    auto left=node(cycle,dna(100,901)),loop=node(cycle,dna(180,902),true),right=node(cycle,dna(100,903));
    path(cycle,{left,loop,loop,right},"CHM13");
    for(int i=0;i<51;++i)path(cycle,i%2?reverse_walk(cycle.paths[0].walk):cycle.paths[0].walk);
    auto cyclic=native_graph_contexts(cycle,0);
    assert(cyclic.size()==1&&cyclic[0].runs[1].copy!=cyclic[0].runs[2].copy);
    for(int threads:{1,2,4,8}) {
        NativeLocalizer index(cyclic,policy,threads);
        const auto middle=cycle.spell_vertex(loop).substr(20,140);
        auto r=index.locate(middle);
        assert(r.status==PlacementStatus::Ambiguous&&r.anchor.status==AnchorStatus::Ambiguous&&r.retained.size()==2);
        same(r,index.locate(reverse_complement(middle)));
        minimizers(middle,policy.k,policy.w,[&](auto x,size_t){assert(index.physical_loci(x)==2);});
        // The exact junction disambiguates the two copies without losing a
        // cross-boundary seed or resetting the full-read extension.
        const auto junction=cycle.spell(0).substr(240,140);
        r=index.locate(junction);
        assert(r.status==PlacementStatus::Unique&&r.best_support==140&&r.anchor.status==AnchorStatus::Known);
        same(r,index.locate(reverse_complement(junction)));
    }
    // A distant donor difference must not multiply repeat loci. The repeated
    // local route and its flanking unique nodes are identical here.
    auto distant=cycle;auto before=node(distant,dna(75,976));
    path(distant,{before,left,loop,loop,right},"distant-variant");
    NativeLocalizer distant_index(native_graph_contexts(distant,0),policy,4);
    const auto interior=cycle.spell_vertex(loop).substr(20,140);
    auto r=distant_index.locate(interior);
    assert(r.status==PlacementStatus::Ambiguous&&r.retained.size()==2);
    minimizers(interior,policy.k,policy.w,[&](auto x,size_t){assert(distant_index.physical_loci(x)==2);});
    const auto junction=cycle.spell(0).substr(240,140);
    r=distant_index.locate(junction);
    assert(r.status==PlacementStatus::Unique&&r.anchor.status==AnchorStatus::Known);
    // Slicing a map and putting it back together preserves physical identity.
    PlacementSequence sliced;sliced.sequence=dna(100,29);sliced.runs={{40,77,0,0,false},{100,77,0,40,false}};
    assert(placement_slices(sliced,20,80)==(std::vector<PlacementSlice>{{77,0,20,80,false}}));
    sliced.runs={{40,77,0,60,true},{100,77,0,0,true}};
    assert(placement_slices(sliced,20,80)==(std::vector<PlacementSlice>{{77,0,20,80,true}}));
    bool rejected=false;try{auto bad=cyclic;bad[0].runs.back().end--;NativeLocalizer index(bad,policy);}catch(const std::invalid_argument&){rejected=true;}assert(rejected);
    rejected=false;try{auto bad=bubble;bad.edges.clear();(void)native_graph_contexts(bad,0);}catch(const std::invalid_argument&){rejected=true;}assert(rejected);
    // Independent base-expanded oracle for physical-locus multiplicities.
    // This deliberately does not call placement_slices or its canonicalizer.
    for(unsigned k:{1u,2u,3u,5u}) {
        auto p=policy;p.k=k;p.w=1;p.maximum_loci=UINT32_MAX;
        NativeGraph small;
        const auto x=node(small,"AACCTGA"),y=node(small,"GCACT",true),z=node(small,"AACCTGA");
        path(small,{x,y,x,z},"reference");path(small,{x,z,y},"different-route");
        path(small,reverse_walk(small.paths[0].walk),"reverse-donor");
        auto contexts=native_graph_contexts(small,0);
        using Base=std::tuple<uint64_t,uint64_t,uint64_t>;
        std::map<uint64_t,std::set<std::vector<Base>>> loci;
        for(const auto&s:contexts) {
            std::vector<Base> bases;uint64_t begin=0;
            for(const auto&r:s.runs) {
                for(uint64_t i=0;i<r.end-begin;++i)
                    bases.emplace_back(r.node,r.copy,r.offset+(r.reverse?r.end-begin-i-1:i));
                begin=r.end;
            }
            words(s.sequence,k,[&](uint64_t w,size_t start,bool) {
                std::vector<Base> span(bases.begin()+start,bases.begin()+start+k);
                auto opposite=span;std::reverse(opposite.begin(),opposite.end());
                loci[w].insert(std::min(span,opposite));
            });
        }
        NativeLocalizer indexed(contexts,p,4);
        for(const auto&[w,locations]:loci)assert(indexed.physical_loci(w)==locations.size());
    }
    std::cout<<"PASS native graph contexts: tiny/reversed nodes, donor aliases, exact reference projections, indels, repeats, cyclic visits and 1/2/4/8 threads\n";
}
