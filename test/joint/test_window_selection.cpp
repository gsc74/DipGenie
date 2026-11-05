#include "TestTemp.hpp"
#include "PanelSelection.hpp"
#include <cassert>
#include <random>
#include <unistd.h>
using namespace hic;
// Synthetic variation graph: reference blocks separated by biallelic SNP
// bubbles. Path 0 is REF; samples carry random ALT sets. Reads are simulated
// error-free from a TARGET diploid genome, never given to selection directly.
struct Toy {
    NativeGraph g;std::vector<std::vector<char>> alt; // per path, per site
    std::vector<uint32_t> ref_block,ref_allele,alt_allele;size_t sites=0;
};
static Toy make_toy(std::mt19937_64&rng,size_t samples,size_t sites,size_t spacing) {
    Toy t;t.sites=sites;
    auto random=[&](size_t n){std::string s(n,'A');for(char&c:s)c="ACGT"[rng()%4];return s;};
    for(size_t i=0;i<=sites;++i){t.ref_block.push_back(uint32_t(t.g.labels.size()));t.g.labels.push_back(random(spacing));}
    for(size_t i=0;i<sites;++i) {
        const char r="ACGT"[rng()%4];char a;do a="ACGT"[rng()%4];while(a==r);
        t.ref_allele.push_back(uint32_t(t.g.labels.size()));t.g.labels.push_back(std::string(1,r));
        t.alt_allele.push_back(uint32_t(t.g.labels.size()));t.g.labels.push_back(std::string(1,a));
    }
    auto add=[&](const std::string&name,const std::vector<char>&alts) {
        NativeGraph::Path p;p.name=name;
        for(size_t i=0;i<sites;++i){p.walk.push_back(2*t.ref_block[i]);p.walk.push_back(2*(alts[i]?t.alt_allele[i]:t.ref_allele[i]));}
        p.walk.push_back(2*t.ref_block[sites]);t.g.paths.push_back(p);t.alt.push_back(alts);
    };
    add("REF",std::vector<char>(sites,0));
    for(size_t s=0;s<samples;++s)for(int h=1;h<=2;++h) {
        std::vector<char> a(sites);for(auto&x:a)x=rng()%3==0;
        add("S"+std::to_string(100+s)+"#"+std::to_string(h)+"#chr",a);
    }
    return t;
}
static void write_reads(std::mt19937_64&rng,const std::vector<std::string>&genome,const std::filesystem::path&dir,size_t pairs) {
    std::ofstream reads[2]={std::ofstream(dir/"R1.fq"),std::ofstream(dir/"R2.fq")};
    for(size_t i=0;i<pairs;++i) {
        const auto&g=genome[rng()%genome.size()];const size_t start=rng()%(g.size()-400);
        const auto fragment=g.substr(start,400);
        const std::string ends[2]={fragment.substr(0,150),reverse_complement(fragment.substr(250,150))};
        for(int side=0;side<2;++side)reads[side]<<'@'<<i<<'/'<<side+1<<'\n'<<ends[side]<<"\n+\n"<<std::string(150,'I')<<'\n';
    }
}
static std::vector<std::string> names(const NativeGraph&g){std::vector<std::string> v;for(auto&p:g.paths)v.push_back(p.name);return v;}
int main() {
    std::mt19937_64 rng(77113);
    std::string pattern_buffer=test_temp_pattern("dipgenie-window");char*pattern=pattern_buffer.data();
    const char*created=mkdtemp(pattern);if(!created)throw std::runtime_error("mkdtemp");
    const std::filesystem::path dir(created);
    // Projection: non-reference vertices take the latest reference coordinate.
    {
        auto t=make_toy(rng,1,3,50);
        std::unordered_map<uint32_t,uint64_t> start;uint64_t p=0;
        for(uint32_t v:t.g.paths[0].walk){start.emplace(v/2,p);p+=t.g.labels[v/2].size();}
        std::vector<char> all(3,1);t.g.paths.push_back(t.g.paths[0]);
        for(size_t i=0;i<3;++i)t.g.paths.back().walk[2*i+1]=2*t.alt_allele[i];
        const auto v=window_select::projection(t.g,t.g.paths.size()-1,start);
        assert(window_select::project(v,10)==10);
        assert(window_select::project(v,50)==50);   // ALT vertex after a 50 bp block
        assert(window_select::project(v,51)==51);   // next reference block resumes exactly
    }
    // 1) Exact recovery: reads from sample S105 select it first.
    auto toy=make_toy(rng,12,400,120);
    auto genome=[&](const std::vector<std::vector<char>>&haps) {
        std::vector<std::string> out;
        for(const auto&a:haps) {
            NativeGraph::Path p;p.name="x";
            for(size_t i=0;i<toy.sites;++i){p.walk.push_back(2*toy.ref_block[i]);p.walk.push_back(2*(a[i]?toy.alt_allele[i]:toy.ref_allele[i]));}
            p.walk.push_back(2*toy.ref_block[toy.sites]);
            NativeGraph g=toy.g;g.paths={p};out.push_back(g.spell(0));
        }
        return out;
    };
    const size_t target=1+2*5;
    write_reads(rng,genome({toy.alt[target],toy.alt[target+1]}),dir,(48000*40)/300);
    std::vector<ReadLibrary> libs={{"sr","SR",(dir/"R1.fq").string(),(dir/"R2.fq").string()}};
    {
        auto g=toy.g;std::vector<size_t> chosen;WindowSelectionOptions o;o.cap=3;o.window_bases=10000;
        const auto s=window_select_panel(g,libs,4,o,&chosen);
        assert((chosen==std::vector<size_t>{0,target,target+1}));
        std::cerr<<"exact: items="<<s.items<<" cost "<<s.initial_cost<<"->"<<s.final_cost<<" mean="<<double(s.mean)<<'\n';
        // Count noise may misclassify a few items; the true pair must remove
        // nearly all read inconsistency of the REF-only panel.
        assert(s.initial_cost>0&&s.final_cost*33<s.initial_cost&&s.steps.size()==1&&s.steps[0].first=="S105");
    }
    // 1b) Low depth (3x): too few reads for any confident hard dosage class, but the
    // summed per-item likelihood still identifies the sample.
    write_reads(rng,genome({toy.alt[target],toy.alt[target+1]}),dir,(48000*3)/300);
    {
        auto g=toy.g;std::vector<size_t> chosen;WindowSelectionOptions o;o.cap=3;o.window_bases=10000;
        const auto s=window_select_panel(g,libs,4,o,&chosen);
        std::cerr<<"low-depth: items="<<s.items<<" cost "<<s.initial_cost<<"->"<<s.final_cost<<" mean="<<double(s.mean)<<'\n';
        assert((chosen==std::vector<size_t>{0,target,target+1}));
        assert(s.final_cost<s.initial_cost);
    }
    // 2) Local mosaic: HAP1 = S102 left half + S108 right half; HAP2 = S110.
    // A pure global-content ranking need not pick all three; the local
    // facility-location objective reaches zero cost with exactly these samples.
    {
        std::vector<char> mosaic(toy.sites);
        const size_t a=1+2*2,b=1+2*8,c=1+2*10;
        for(size_t i=0;i<toy.sites;++i)mosaic[i]=i<toy.sites/2?toy.alt[a][i]:toy.alt[b][i];
        write_reads(rng,genome({mosaic,toy.alt[c]}),dir,(48000*40)/300);
        auto g=toy.g;std::vector<size_t> chosen;WindowSelectionOptions o;o.cap=7;o.window_bases=12000;
        const auto s=window_select_panel(g,libs,4,o,&chosen);
        assert((chosen==std::vector<size_t>{0,a,a+1,b,b+1,c,c+1}));
        std::cerr<<"mosaic: items="<<s.items<<" cost "<<s.initial_cost<<"->"<<s.final_cost<<'\n';
        assert(s.final_cost*33<s.initial_cost); // one window spans the mosaic switch
        // Worker-count parity of the complete selection and audit.
        for(int threads:{1,2,8}) {
            auto h=toy.g;std::vector<size_t> again;std::vector<long double> gains,base;
            auto x=window_select_panel(h,libs,threads,o,&again,&gains);
            assert(again==chosen&&x.final_cost==s.final_cost&&x.items==s.items&&x.steps==s.steps);
        }
        // Integrated wrapper: complete pairs, REF first, audit written.
        auto w=toy.g;const auto out=dir/"select";std::filesystem::create_directory(out);
        select_panel(w,libs,out,2,7,12000);
        assert(w.paths.size()==7&&w.paths[0].name=="REF");
        for(size_t i=1;i<7;i+=2)assert(w.paths[i].name.substr(0,4)==w.paths[i+1].name.substr(0,4));
        assert(std::filesystem::file_size(out/"panel-selection.tsv")>0&&std::filesystem::file_size(out/"window-selection.tsv")>0);
        std::filesystem::remove_all(out);
    }
    // 3) A panel no larger than the cap is unchanged.
    {
        auto g=toy.g;g.paths.resize(7);const auto before=names(g);
        select_panel(g,libs,dir,1,7,12000);assert(names(g)==before);
    }
    std::filesystem::remove_all(dir);
    std::cout<<"PASS local window panel selection: projection, exact sample recovery, mosaic recovery, worker parity, audit\n";
}
