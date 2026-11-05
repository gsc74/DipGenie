// Independent frozen serial construction: test oracle, never an inference path.
#include "SequenceEvidence.hpp"
#include <cassert>
#include <iostream>
#include <random>
#include <chrono>
using namespace hic;
static CommonSegments serial_segments(const NativeGraph&g,size_t k,int threads=1) {
        if(k<2||threads<1)throw std::invalid_argument("segment options");
        const size_t H=g.paths.size(),V=g.labels.size();
        std::vector<bool> common(2*V,true);
        for(auto&p:g.paths) {
            std::vector<unsigned char> count(V,0);std::vector<bool> present(2*V,false);
            for(uint32_t v:p.walk){if(count[v/2]<2)++count[v/2];present[v]=true;}
            for(size_t v=0;v<2*V;++v)common[v]=common[v]&&present[v]&&count[v/2]==1&&g.labels[v/2].size()>=k-1;
        }
        std::vector<std::unordered_map<uint32_t,size_t>> positions(H);
        for(size_t h=0;h<H;++h)for(size_t i=0;i<g.paths[h].walk.size();++i) {
            uint32_t v=g.paths[h].walk[i];if(common[v])positions[h][v]=i;
        }
        CommonSegments result;std::vector<size_t> previous(H,0);bool first=true;
        for(uint32_t v:g.paths[0].walk)if(common[v]) {
            bool valid=true;size_t terminal=0;
            for(size_t h=0;h<H;++h){valid=valid&&(first||positions[h].at(v)>previous[h]);terminal+=positions[h].at(v)+1==g.paths[h].walk.size();}
            if(terminal&&terminal<H&&g.labels[v/2].size()==k-1)valid=false;
            if(valid){result.anchors.push_back(v);for(size_t h=0;h<H;++h)previous[h]=positions[h].at(v);first=false;}
        }
        result.cuts.resize(H);result.sequences.resize(H);
        std::unordered_set<uint32_t> selected(result.anchors.begin(),result.anchors.end());
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(threads) schedule(dynamic,1)
        #endif
        for(size_t h=0;h<H;++h) {
            auto&cuts=result.cuts[h];cuts.push_back(0);size_t offset=0;
            for(uint32_t v:g.paths[h].walk) {
                if(selected.count(v))cuts.push_back(offset+k-1);
                result.sequences[h]+=g.spell_vertex(v);offset+=g.labels[v/2].size();
            }
            cuts.push_back(offset);
        }
        for(auto&cuts:result.cuts)if(cuts.size()!=result.anchors.size()+2||!std::is_sorted(cuts.begin(),cuts.end()))
            throw std::runtime_error("Inconsistent common anchors");
        return result;
    }
static SequenceIndex serial_index(const NativeGraph&g,const CommonSegments&parts,unsigned k=31,unsigned w=25) {
        (void)checked_pair_domain(g.paths.size());
        SequenceIndex out;out.k=k;out.w=w;std::unordered_set<uint64_t> selected;
        for(auto&s:parts.sequences)minimizers(s,k,w,[&](uint64_t x,size_t){selected.insert(x);});
        std::vector<uint64_t> order(selected.begin(),selected.end());std::sort(order.begin(),order.end());
        out.ids.reserve(order.size());out.targets.reserve(order.size());
        for(auto x:order){out.ids.emplace(x,out.targets.size());WordTarget t;t.word=x;out.targets.push_back(std::move(t));}
        std::vector<uint32_t> count(order.size());
        for(size_t h=0;h<g.paths.size();++h) {
            std::fill(count.begin(),count.end(),0);
            // Physical oriented base coordinates distinguish different loci,
            // including k-mers crossing small graph vertices.
            std::vector<uint64_t> coordinates;coordinates.reserve(parts.sequences[h].size());
            for(uint32_t v:g.paths[h].walk) {
                size_t length=g.labels[v/2].size();if(length>UINT32_MAX)throw std::overflow_error("vertex coordinates");
                for(uint32_t p=0;p<length;++p)coordinates.push_back((uint64_t(v)<<32)|p);
            }
            auto flipped=[&](uint64_t position){uint32_t v=position>>32,p=uint32_t(position);return(uint64_t(v^1)<<32)|(g.labels[v/2].size()-1-p);};
            size_t block=0;
            words(parts.sequences[h],k,[&](uint64_t word,size_t p,bool reverse){
                auto found=out.ids.find(word);if(found==out.ids.end())return;
                while(p+k>parts.cuts[h][block+1])++block;
                auto&t=out.targets[found->second];t.carriers.set(h);
                auto it=std::lower_bound(t.occurrences.begin(),t.occurrences.end(),int(block),[](const auto&a,int b){return a.first<b;});
                if(it==t.occurrences.end()||it->first!=int(block))it=t.occurrences.insert(it,{int(block),DonorMask{}});
                it->second.set(h);
                auto&copies=count[found->second];
                if(copies==UINT32_MAX)throw std::overflow_error("Panel string multiplicity");
                if(copies++)t.unique_locus=false;
                uint64_t a=reverse?flipped(coordinates[p+k-1]):coordinates[p];
                uint64_t b=reverse?flipped(coordinates[p]):coordinates[p+k-1];
                if(!t.location_set){t.location_set=true;t.location_start=a;t.location_end=b;}
                else if(t.location_start!=a||t.location_end!=b)t.unique_locus=false;
            });
            for(size_t x=0;x<count.size();++x)if(count[x]) {
                auto&hist=out.targets[x].copy_histogram;
                auto it=std::lower_bound(hist.begin(),hist.end(),count[x],[](auto a,uint32_t n){return a.first<n;});
                if(it==hist.end()||it->first!=count[x])hist.insert(it,{count[x],1});
                else ++it->second;
            }
        }
        return out;
    }

static void same(const SequenceIndex&a,const SequenceIndex&b) {
    assert(a.k==b.k&&a.w==b.w&&a.ids==b.ids&&a.targets.size()==b.targets.size());
    for(size_t i=0;i<a.targets.size();++i) {
        const auto&x=a.targets[i];const auto&y=b.targets[i];
        assert(x.word==y.word&&x.count==y.count&&x.carriers==y.carriers);
        assert(x.location_set==y.location_set&&x.unique_locus==y.unique_locus);
        assert(x.location_start==y.location_start&&x.location_end==y.location_end);
        assert(x.occurrences==y.occurrences&&x.copy_histogram==y.copy_histogram);
    }
}
static void check(const NativeGraph&g,unsigned k,unsigned w) {
    const auto parts=CommonSegments::build(g,std::max(2u,k),1);
    const auto reference_parts=serial_segments(g,std::max(2u,k),1);
    assert(parts.anchors==reference_parts.anchors&&parts.cuts==reference_parts.cuts&&parts.sequences==reference_parts.sequences);
    const auto expected=serial_index(g,parts,k,w);
    for(int t:{1,2,4,8}) {
        const auto parallel_parts=CommonSegments::build(g,std::max(2u,k),t);
        assert(parts.anchors==parallel_parts.anchors&&parts.cuts==parallel_parts.cuts&&parts.sequences==parallel_parts.sequences);
        same(expected,SequenceIndex::build(g,parallel_parts,k,w,t));
    }
}
int main() {
    std::mt19937_64 rng(87231);
    auto sequence=[&](size_t n){std::string s;for(size_t i=0;i<n;++i)s+="ACGTN"[rng()%5];return s;};
    for(int trial=0;trial<36;++trial) {
        NativeGraph g;
        g.labels={"ACGTACGTACGTACGTACGTACGTACGTACGTACGT",sequence(50),"AAAAAAAAAAAAAAAAAAAA","CGCGCGCGCGCGCGCGCGC"};
        for(int i=0;i<12;++i)g.labels.push_back(sequence(1+rng()%37));
        for(int h=0;h<4;++h) {
            NativeGraph::Path p;p.name=std::to_string(h);p.walk.push_back(0);
            for(int i=0;i<25;++i)p.walk.push_back(2+uint32_t(rng()%(2*g.labels.size()-2)));
            p.walk.push_back(2);
            for(size_t i=1;i<p.walk.size();++i)g.edges.insert(NativeGraph::edge(p.walk[i-1],p.walk[i]));
            g.paths.push_back(std::move(p));
        }
        check(g,1,1);check(g,3,1);check(g,7,5);check(g,31,25);
    }
    // More than two discovery chunks: minimizer windows and k-mer endpoints
    // straddle both chunk boundaries. Ambiguity, repeat and reverse tests.
    NativeGraph large;std::string long_sequence((1u<<21)+103,'A');
    for(char&c:long_sequence)c="ACGT"[rng()%4];
    for(size_t p:{size_t(1024),size_t((1u<<20)-15),size_t((1u<<20)+17)})
        long_sequence.replace(p,7,"NNNNNNN");
    large.labels={long_sequence,"ACGTACGTACGTACGTACGTACGTACGTACGTACGT"};
    large.paths={{"forward",{0,2,0}},{"reverse",{1,2,0}},{"short", {2}}};
    check(large,31,25);
    // Review witness H1=(a,v), H2=(b,v,c), |v|=k-1: v ends H1 only, so a cut
    // inside v would not be strictly interior on H1. It is not an anchor.
    {
        NativeGraph w;w.labels={"C","G","AA","T"};w.paths={{"H1",{0,4}},{"H2",{2,4,6}}};
        const auto parts=CommonSegments::build(w,3,1);
        assert(parts.anchors.empty()&&parts.blocks()==1);
        const auto oracle=serial_segments(w,3,1);
        assert(oracle.anchors==parts.anchors&&oracle.cuts==parts.cuts);
        NativeGraph both;both.labels={"C","G","AA"};both.paths={{"H1",{0,4}},{"H2",{2,4}}};
        assert(CommonSegments::build(both,3,1).anchors.size()==1);
    }
    // No complete minimizer windows, including a single-nucleotide panel.
    NativeGraph tiny;tiny.labels={"A"};tiny.paths={{"one",{0}}};
    check(tiny,31,25);
    bool rejected=false;
    auto parts=CommonSegments::build(tiny,31);
    try{(void)SequenceIndex::build(tiny,parts,0,25,4);}catch(const std::invalid_argument&){rejected=true;}
    assert(rejected);
    rejected=false;
    try{(void)SequenceIndex::build(tiny,parts,31,0,4);}catch(const std::invalid_argument&){rejected=true;}
    assert(rejected);
    std::cout<<"PASS serial/parallel sequence-index identity: keys, loci, strands, carrier masks, copy counts, sorted segment occurrences; 1/2/4/8 workers; small nodes, cycles, ambiguity, repeats, chunk boundaries, empty targets\n";
}
