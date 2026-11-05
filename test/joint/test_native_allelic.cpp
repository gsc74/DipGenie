#include "NativeReadEvidence.hpp"
#include <cassert>
#include <iostream>
#include <random>

using namespace hic;
static std::string dna(size_t n,uint64_t seed) {
    std::mt19937_64 rng(seed);std::string s(n,'A');
    for(auto&c:s)c="ACGT"[rng()%4];
    return s;
}
static uint32_t node(NativeGraph&g,const std::string&s) {
    const auto id=uint32_t(g.labels.size());g.names.push_back(std::to_string(id));
    g.labels.push_back(s);return 2*id;
}
static void path(NativeGraph&g,std::vector<uint32_t>w) {
    for(size_t i=1;i<w.size();++i) {
        g.edges.insert(NativeGraph::edge(w[i-1],w[i]));
        g.edges.insert(NativeGraph::edge(w[i]^1,w[i-1]^1));
    }
    g.paths.push_back({"path"+std::to_string(g.paths.size()),std::move(w)});
}
static void same(const EndPlacement&a,const EndPlacement&b) {
    assert(a.status==b.status&&a.component==b.component);
    assert(a.best_support==b.best_support&&a.competing_support==b.competing_support);
    assert(a.anchor.status==b.anchor.status&&a.retained.size()==b.retained.size());
    if(a.anchor.status==AnchorStatus::Known)
        assert(std::tie(a.anchor.contig,a.anchor.start,a.anchor.word)==
               std::tie(b.anchor.contig,b.anchor.start,b.anchor.word));
}
int main() {
    NativeGraph graph;
    const auto left=node(graph,dna(500,891)),a=node(graph,"A"),c=node(graph,"C"),
               right=node(graph,dna(500,892));
    path(graph,{left,a,right});path(graph,{left,c,right});
    NativeLocalizer::Policy policy; // production k31/w25, 250-base read
    for(int threads:{1,2,4,8}) {
        NativeLocalizer index(native_graph_contexts(graph,0),policy,threads);
        // A one-SNP alternative must not accumulate 249 bases from two
        // disjoint extensions and compete with the exact 250-base route.
        for(size_t h:{size_t(0),size_t(1)}) {
            const auto read=graph.spell(h).substr(400,250);
            const auto placed=index.locate(read);
            assert(placed.status==PlacementStatus::Unique&&placed.best_support==250);
            assert(placed.anchor.status==AnchorStatus::Known);
            same(placed,index.locate(reverse_complement(read)));
            auto error=read;error[40]=error[40]=='A'?'G':'A';
            const auto imperfect=index.locate(error);
            // The 40-base error-side flank need not contain a selected seed;
            // the 209-base SNP-bearing side must still localize correctly.
            assert(imperfect.status==PlacementStatus::Unique&&imperfect.best_support>=209&&imperfect.best_support<=249);
            for(const auto&span:imperfect.candidates.front().intervals)
                assert(span.end<=40||span.begin>40);
            same(imperfect,index.locate(reverse_complement(error)));
        }
        auto unknown=graph.spell(0).substr(400,250);unknown[100]='N';
        const auto unresolved=index.locate(unknown);
        assert(unresolved.status==PlacementStatus::Ambiguous);
        same(unresolved,index.locate(reverse_complement(unknown)));
        // Same-locus maximum ties are retained, regardless of donor count.
        auto copies=graph;for(int n=0;n<51;++n)path(copies,graph.paths[0].walk);
        NativeLocalizer copied(native_graph_contexts(copies,0),policy,threads);
        same(unresolved,copied.locate(unknown));
        same(index.locate(graph.spell(1).substr(400,250)),copied.locate(graph.spell(1).substr(400,250)));
        // Near-identical PARALOGS remain ambiguous. They have different
        // physical seed loci even though their DNA differs by only one base.
        NativeGraph paralogs;
        auto p0=node(paralogs,graph.spell(0)),p1=node(paralogs,graph.spell(1));
        path(paralogs,{p0});path(paralogs,{p1});
        NativeLocalizer repeated(native_graph_contexts(paralogs,0),policy,threads);
        const auto multi=repeated.locate(graph.spell(0).substr(400,250));
        assert(multi.status==PlacementStatus::Ambiguous&&multi.best_support==250&&multi.competing_support==249);
        same(multi,repeated.locate(reverse_complement(graph.spell(0).substr(400,250))));
    }
    // Paired allelic reads carry their exact marker identities into the
    // objective. This is not a whole-read phase label guessed by the mapper.
    const auto mid=node(graph,dna(900,893)),g=node(graph,"G"),t=node(graph,"T"),
               tail=node(graph,dna(500,894));
    graph.paths.clear();path(graph,{left,a,mid,g,tail});path(graph,{left,c,mid,t,tail});
    NativeLocalizer index(native_graph_contexts(graph,0),policy,4);
    SequenceIndex targets;targets.k=31;targets.w=25;
    for(size_t h=0;h<2;++h)minimizers(graph.spell(h),31,25,[&](uint64_t word,size_t) {
        if(targets.ids.count(word))return;
        targets.ids[word]=targets.targets.size();WordTarget target;target.word=word;
        target.unique_locus=true;targets.targets.push_back(target);
    });
    for(size_t h=0;h<2;++h) {
        FastqRead one,two;one.sequence=graph.spell(h).substr(400,250);two.sequence=graph.spell(h).substr(1300,250);
        const auto base=native_pair_evidence(targets,index,one,two,"pair",MoleculeKind::HiC,0);
        assert(!base.marker.ends[0].empty()&&!base.marker.ends[1].empty());
        assert(base.marker.geometry==MoleculeGeometry::Known);
        for(bool r1:{false,true})for(bool r2:{false,true})for(bool swap:{false,true}) {
            auto x=one,y=two;if(r1)x.sequence=reverse_complement(x.sequence);if(r2)y.sequence=reverse_complement(y.sequence);
            if(swap)std::swap(x,y);
            const auto placed=native_pair_evidence(targets,index,x,y,"pair",MoleculeKind::HiC,0);
            assert(placed.marker.reference_distance==base.marker.reference_distance);
            for(int e=0;e<2;++e)assert(placed.marker.ends[e]==base.marker.ends[swap?1-e:e]);
        }
    }
    std::cout<<"PASS native allelic extension: exact SNPs, errors, unresolved alleles, donor copies, true paralogs, paired marker retention and 1/2/4/8-thread/RC parity\n";
}
