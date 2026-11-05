#include "NativeLocalizer.hpp"
#include <cassert>
#include <iostream>
#include <random>

using namespace hic;
static std::string dna(size_t n,uint64_t seed) {
    std::mt19937_64 rng(seed);std::string s(n,'A');for(auto&c:s)c="ACGT"[rng()%4];return s;
}
static PlacementSequence context(std::string sequence,uint32_t component=0,uint32_t space=0,uint64_t origin=0,bool reference=true) {
    PlacementSequence s;s.name="context";s.sequence=std::move(sequence);s.component=component;s.space=space;s.origin=origin;
    s.reference_contig=reference?component:UINT32_MAX;return s;
}
static void same_decision(const EndPlacement&a,const EndPlacement&b) {
    assert(a.status==b.status&&a.component==b.component&&a.best_support==b.best_support&&a.competing_support==b.competing_support);
    assert(a.anchor.status==b.anchor.status);
    if(a.anchor.status==AnchorStatus::Known)assert(std::tie(a.anchor.contig,a.anchor.start,a.anchor.word)==std::tie(b.anchor.contig,b.anchor.start,b.anchor.word));
    assert(a.candidates.size()==b.candidates.size()&&a.retained==b.retained);
}
int main() {
    NativeLocalizer::Policy p;p.k=15;p.w=7;
    const auto reference=dna(1400,127),other=dna(1600,161);
    // Build at each thread count. The inputs and all tie decisions are fixed.
    for(int threads:{1,2,4,8}) {
        NativeLocalizer index({context(reference),context(other,1,1)},p,threads);
        for(size_t start:{size_t(30),size_t(270),size_t(940)}) {
            auto read=reference.substr(start,150);
            auto a=index.locate(read),b=index.locate(reverse_complement(read));
            same_decision(a,b);
            assert(a.status==PlacementStatus::Unique&&a.component==0&&a.best_support==150);
            assert(a.anchor.status==AnchorStatus::Known&&a.anchor.start>=start&&a.anchor.start<=start+150-p.k);
            assert(a.candidates.front().reverse!=b.candidates.front().reverse);
            // A selected anchor is an exact reference word, not a fitted start.
            uint64_t exact=UINT64_MAX;words(reference.substr(a.anchor.start,p.k),p.k,[&](auto x,size_t,bool){exact=x;});
            assert(exact==a.anchor.word);
            read[77]=read[77]=='A'?'C':'A';
            same_decision(index.locate(read),index.locate(reverse_complement(read)));
        }
        auto first=index.locate(reference.substr(30,150)),second=index.locate(reference.substr(700,150));
        auto d=selected_anchor_distance(first,second);
        assert(d.status==AnchorDistanceStatus::Known&&d.separation==second.anchor.start-first.anchor.start);
        assert(selected_anchor_distance(second,first).separation==d.separation);
        for(bool a:{false,true})for(bool b:{false,true}) {
            auto r1=reference.substr(30,150),r2=reference.substr(700,150);
            if(a)r1=reverse_complement(r1);
            if(b)r2=reverse_complement(r2);
            assert(selected_anchor_distance(index.locate(r1),index.locate(r2)).separation==d.separation);
        }
        assert(selected_anchor_distance(first,index.locate(other.substr(20,150))).status==AnchorDistanceStatus::Trans);
        assert(selected_anchor_distance(first,index.locate("N")).status==AnchorDistanceStatus::Unknown);
        assert(index.locate("NNN").status==PlacementStatus::NoSeeds);
        auto chimera=reference.substr(30,90)+other.substr(400,90);
        auto c=index.locate(chimera);
        assert(c.status==PlacementStatus::Split&&!c.component_unique()&&c.anchor.status==AnchorStatus::Split);
        same_decision(c,index.locate(reverse_complement(chimera)));
    }
    // Donor copies of one physical locus do not cause repetitive-seed rejection.
    p.maximum_loci=1;
    auto s=context(reference);
    std::vector<PlacementSequence> copies(51,s);
    NativeLocalizer one(copies,p,4);
    auto r=one.locate(reference.substr(300,150));
    assert(r.status==PlacementStatus::Unique&&r.best_support==150);
    uint64_t key=0;minimizers(reference,p.k,p.w,[&](auto x,size_t){key=x;});
    assert(one.physical_loci(key)==1);
    // The same string at another locus is not a donor duplicate.
    auto duplicate=s;duplicate.space=9;
    NativeLocalizer repeated({s,duplicate},p,2);
    assert(repeated.physical_loci(key)==2);
    assert(repeated.locate(reference.substr(300,150)).status==PlacementStatus::NoSeeds);
    // With repetitive seeds retained: same component, ambiguous physical locus.
    p.maximum_loci=64;
    NativeLocalizer ambiguous({s,duplicate},p,2);
    r=ambiguous.locate(reference.substr(300,150));
    assert(r.status==PlacementStatus::Ambiguous&&r.component_unique()&&r.anchor.status==AnchorStatus::Ambiguous);
    same_decision(r,ambiguous.locate(reverse_complement(reference.substr(300,150))));
    // Reverse-oriented context and overlapping chunks agree on physical starts.
    auto reversed=s;reversed.reverse=true;reversed.sequence=reverse_complement(reference);
    NativeLocalizer oriented({reversed},p,2);
    same_decision(oriented.locate(reference.substr(300,150)),NativeLocalizer({s},p,1).locate(reference.substr(300,150)));
    const size_t cut=600,halo=200; // >= read length - 1
    NativeLocalizer chunks({context(reference.substr(0,cut+halo)),context(reference.substr(cut-halo),0,0,cut-halo)},p,2);
    r=chunks.locate(reference.substr(cut-80,150));
    assert(r.status==PlacementStatus::Unique&&r.best_support==150);
    same_decision(r,chunks.locate(reverse_complement(reference.substr(cut-80,150))));
    // No reference projection: valid relation-only localization, NOT distance 0.
    NativeLocalizer novel({context(reference,0,3,0,false)},p,1);
    r=novel.locate(reference.substr(300,150));
    assert(r.status==PlacementStatus::Unique&&r.anchor.status==AnchorStatus::Unprojected);
    assert(selected_anchor_distance(r,r).status==AnchorDistanceStatus::Unknown);
    // Filtering before identity deduplication, including the reverse range.
    SequenceIndex targets;targets.k=3;targets.w=1;targets.ids.emplace(0,0);targets.targets.resize(1);
    targets.targets[0].word=0;targets.targets[0].unique_locus=true;targets.targets[0].count=43;targets.targets[0].occurrence_count=77;
    assert(native_interval_markers(targets,"AAAAAA",{{1,6}})==std::vector<size_t>{0});
    assert(native_interval_markers(targets,"TTTTTT",{{0,5}})==std::vector<size_t>{0});
    assert(native_interval_markers(targets,"AAAAAA",{{1,3}}).empty());
    assert(targets.targets[0].count==43&&targets.targets[0].occurrence_count==77);
    NativeLocalizer::Policy small;small.k=3;small.w=1;small.maximum_loci=1;
    NativeLocalizer poly({context("AAAAAA")},small,1);
    assert(poly.physical_loci(0)==4); // not one collapsed minimizer identity
    assert(poly.locate("AAAAAA").status==PlacementStatus::NoSeeds);
    // Even-k palindromes carry both strands; never invent a unique orientation.
    small.k=2;small.minimum_seeds=1;small.maximum_loci=8;
    NativeLocalizer palindrome({context("AT")},small,1);
    assert(palindrome.locate("AT").status==PlacementStatus::Ambiguous);
    // Internal construction chunks partition windows/starts, not sequence
    // context. Full-length extension crosses every chunk boundary unchanged.
    const auto long_reference=dna(600000,11271);
    NativeLocalizer serial({context(long_reference)},p,1),parallel({context(long_reference)},p,8);
    assert(serial.stats().context_bases==long_reference.size());
    assert(serial.stats().ingest_batches>1);
    assert(serial.stats().selected_words==parallel.stats().selected_words);
    assert(serial.stats().retained_occurrences==parallel.stats().retained_occurrences);
    assert(serial.stats().physical_loci==parallel.stats().physical_loci);
    assert(serial.stats().peak_ingest_scratch_bytes<=NativeLocalizer::ingest_scratch_limit);
    assert(parallel.stats().peak_ingest_scratch_bytes<=NativeLocalizer::ingest_scratch_limit);
    for(size_t boundary:{size_t(1u<<16),size_t(2u<<16),size_t(1u<<18),size_t(2u<<18)}) {
        const auto read=long_reference.substr(boundary-75,150);
        auto a=serial.locate(read),b=parallel.locate(read);
        same_decision(a,b);assert(a.best_support==150&&a.status==PlacementStatus::Unique);
        same_decision(a,parallel.locate(reverse_complement(read)));
    }
    std::unordered_set<uint64_t> whole_selected;
    minimizers(long_reference,p.k,p.w,[&](auto word,size_t){whole_selected.insert(word);});
    std::unordered_map<uint64_t,size_t> expected_counts;
    words(long_reference,p.k,[&](auto word,size_t,bool){if(whole_selected.count(word))++expected_counts[word];});
    for(auto [word,n]:expected_counts)assert(parallel.physical_loci(word)==n);
    bool caught=false;try{one.locate(std::string(p.maximum_read_length+1,'A'));}catch(const std::invalid_argument&){caught=true;}
    assert(caught);
    std::cout<<"PASS native exact-extension localizer: RC/strand and anchor/mate invariance; physical locus versus donor copies; split/trans/unknown states; chunk continuity; filtering before identity collapse; 1/2/4/8 threads\n";
}
