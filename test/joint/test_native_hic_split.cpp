#include "NativeReadEvidence.hpp"
#include <cassert>
#include <iostream>
#include <random>

using namespace hic;
static std::string dna(size_t n,uint64_t seed) {
    std::mt19937_64 rng(seed);std::string s(n,'A');for(auto&c:s)c="ACGT"[rng()%4];return s;
}
static PlacementCandidate candidate(uint32_t component,uint32_t begin,uint32_t end,uint32_t support,uint64_t anchor) {
    PlacementCandidate c;c.component=component;c.support=support;c.intervals={{begin,end}};
    c.anchor.status=AnchorStatus::Known;c.anchor.contig=0;c.anchor.start=anchor;c.anchor.word=anchor;
    return c;
}
static EndPlacement end_with(PlacementStatus status,std::vector<PlacementCandidate> candidates) {
    EndPlacement e;e.status=status;e.candidates=std::move(candidates);
    e.anchor.status=status==PlacementStatus::Split?AnchorStatus::Split:AnchorStatus::InsufficientSupport;
    return e;
}
static SequenceIndex target_index(const std::vector<std::string>&sequences,unsigned k,unsigned w) {
    SequenceIndex index;index.k=k;index.w=w;std::set<uint64_t> selected;
    for(const auto&s:sequences)minimizers(s,k,w,[&](auto word,size_t){selected.insert(word);});
    for(auto word:selected){index.ids[word]=index.targets.size();WordTarget t;t.word=word;index.targets.push_back(t);}
    return index;
}
int main() {
    const NativePairPolicy on;
    {
        auto e=end_with(PlacementStatus::Split,{candidate(0,80,150,90,9000),candidate(0,0,80,80,100)});
        assert(resolve_hic_five_prime(e,on));
        assert(e.status==PlacementStatus::Unique&&e.retained==std::vector<size_t>{1});
        assert(e.anchor.status==AnchorStatus::Known&&e.anchor.start==100&&e.best_support==80&&e.component==0);
    }
    // Partial support (no second candidate) resolves likewise.
    {
        auto e=end_with(PlacementStatus::InsufficientSupport,{candidate(0,0,70,70,500)});
        assert(resolve_hic_five_prime(e,on)&&e.anchor.start==500);
    }
    // A competing placement of the same read span with >= 80% support: reject.
    {
        auto e=end_with(PlacementStatus::Split,{candidate(0,0,80,80,100),candidate(0,5,80,64,7000),candidate(0,80,150,70,9000)});
        const auto before=e.status;
        assert(!resolve_hic_five_prime(e,on)&&e.status==before&&e.anchor.status==AnchorStatus::Split);
    }
    // A weaker competitor (< 80%) or one covering under half of the span: accept.
    {
        auto e=end_with(PlacementStatus::Split,{candidate(0,0,80,80,100),candidate(0,5,80,63,7000)});
        assert(resolve_hic_five_prime(e,on)&&e.anchor.start==100);
        auto f=end_with(PlacementStatus::Split,{candidate(0,0,80,80,100),candidate(0,50,130,80,7000)});
        assert(resolve_hic_five_prime(f,on)&&f.anchor.start==100);
    }
    // Too little exact support at the 5' end: reject, even if the 3' part is large.
    {
        auto e=end_with(PlacementStatus::Split,{candidate(0,0,59,59,100),candidate(0,59,150,91,9000)});
        assert(!resolve_hic_five_prime(e,on));
    }
    // Equal 5' starts: larger support wins.
    {
        auto e=end_with(PlacementStatus::Split,{candidate(0,0,75,75,100),candidate(0,0,90,90,300)});
        // The two cover the same span with 75 >= 0.8*90 support: genuinely ambiguous.
        assert(!resolve_hic_five_prime(e,on));
        auto f=end_with(PlacementStatus::Split,{candidate(0,0,70,60,100),candidate(0,0,90,90,300)});
        assert(resolve_hic_five_prime(f,on)&&f.anchor.start==300);
    }
    // Unique, ambiguous and seedless ends are never changed.
    for(auto status:{PlacementStatus::Unique,PlacementStatus::Ambiguous,PlacementStatus::NoSeeds}) {
        EndPlacement e;e.status=status;e.candidates={candidate(0,0,80,80,100)};
        assert(!resolve_hic_five_prime(e,on)&&e.status==status);
    }
    // End to end through the localizer: a Hi-C read whose 5' 80 bases come from
    // the reference and whose 3' 70 bases come from another component.
    NativeLocalizer::Policy policy;policy.k=15;policy.w=7;
    const auto ref=dna(12000,1234),other=dna(4000,3456);
    NativeGraph g;g.names={"a"};g.labels={ref};g.paths={{"CHM13",{0}}};
    auto contexts=native_graph_contexts(g,0);
    PlacementSequence outside;outside.name="other";outside.sequence=other;outside.component=1;outside.space=100;outside.reference_contig=1;
    contexts.push_back(std::move(outside));
    NativeLocalizer localizer(contexts,policy,2);
    const auto targets=target_index({ref,other},policy.k,policy.w);
    FastqRead a{"pair0",ref.substr(4500,80)+other.substr(500,70),std::string(150,'I'),1};
    FastqRead b{"pair0",ref.substr(1000,150),std::string(150,'I'),2};
    const auto resolved=native_pair_evidence(targets,localizer,a,b,"pair0",MoleculeKind::HiC,0,on,true);
    assert(resolved.audit.find("\taccepted\t")!=std::string::npos&&resolved.audit.find("five_prime,")!=std::string::npos);
    assert(!resolved.marker.ends[0].empty()&&!resolved.marker.ends[1].empty());
    assert(resolved.marker.geometry==MoleculeGeometry::Known&&resolved.marker.reference_distance>0);
    // Markers come from the 5' (reference) interval only.
    std::set<size_t> ref_markers;
    minimizers(ref,policy.k,policy.w,[&](auto word,size_t){auto it=targets.ids.find(word);if(it!=targets.ids.end())ref_markers.insert(it->second);});
    for(size_t x:resolved.marker.ends[0])assert(ref_markers.count(x));
    // The rule never applies to short reads.
    const auto sr=native_pair_evidence(targets,localizer,a,b,"pair0",MoleculeKind::ShortRead,0,on,true);
    assert(sr.audit.find("\tunresolved_end\t")!=std::string::npos);
    std::cout<<"PASS Hi-C 5' junction resolution: 5'-most selection, support and competitor gates, ties, unchanged other states, localizer integration, 5'-only markers, short reads unaffected\n";
}
