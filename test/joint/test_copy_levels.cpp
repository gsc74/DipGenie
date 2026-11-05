#include "TestTemp.hpp"
#include "BuildModel.hpp"
#include "NativeGraph.hpp"
#include <cassert>
#include <filesystem>
#include <iostream>
#include <random>
#include <unistd.h>
using namespace hic;
using hic_calibration::copy_level_posterior;
static hic_calibration::CompleteDosage model(long double mean,long double shape) {
    hic_calibration::CompleteDosage d;d.calibration={mean,mean+mean*mean/shape,shape};
    d.error_rate=.5L;d.error_zero=.95L;d.dropout=0;return d;
}
int main() {
    // 1) Posterior identities.
    {
        const auto d=model(40,20);
        const std::map<uint32_t,uint64_t> classes{{0,7},{1,5},{2,4},{4,3}};
        for(uint64_t n:{0ULL,3ULL,19ULL,41ULL,77ULL,160ULL,300ULL}) {
            const auto levels=copy_level_posterior(d,n,classes);
            assert(levels.size()==4);
            long double expected_total=0,level_sum=0,top=-INFINITY,norm=0;
            std::map<uint32_t,long double> joint;
            for(auto[i,a]:classes)for(auto[j,b]:classes){const auto v=std::log((long double)a*b)+d.copies_log(n,i+j);joint[i*100+j]=v;top=std::max(top,v);}
            for(auto&[k,v]:joint){v=std::exp(v-top);norm+=v;}
            for(auto&[k,v]:joint)expected_total+=(k/100+k%100)*v/norm;
            for(size_t k=0;k<levels.size();++k) {
                const auto&p=levels[k];
                assert(std::fabs(p[0]+p[1]+p[2]-1)<1e-12L&&p[0]>=0&&p[1]>=0&&p[2]>=0);
                level_sum+=p[1]+2*p[2];
                if(k)assert(levels[k][0]>=levels[k-1][0]-1e-15L); // nested levels
            }
            assert(std::fabs(level_sum-expected_total)<1e-9L); // sum_k E[d_k] = E[i+j]
        }
        // Strong count evidence for total six with classes {0,2,4}: only (2,4).
        const auto p=copy_level_posterior(model(40,1e6),120,{{0,10},{2,10},{4,10}});
        assert(p[0][2]>.99L&&p[1][2]>.99L&&p[2][1]>.99L&&p[3][1]>.99L);
        // One copy level reduces to the ordinary panel-prior three-state posterior.
        const auto q=copy_level_posterior(d,25,{{0,3},{1,2}});
        long double w[3]={9*std::exp(d.copies_log(25,0)),12*std::exp(d.copies_log(25,1)),4*std::exp(d.copies_log(25,2))};
        const long double s=w[0]+w[1]+w[2];
        for(int k=0;k<3;++k)assert(std::fabs(q[0][k]-w[k]/s)<1e-12L);
    }
    // 2) Raw-read tandem repeat: donors carry 2, 3 or 6 units; HG has 3+6.
    // Every repeat string is present in every allele, so presence cannot
    // identify the pair; total copy number 9 uniquely identifies {3,6}.
    std::mt19937_64 rng(5150);
    auto random=[&](size_t n){std::string s(n,'A');for(char&c:s)c="ACGT"[rng()%4];return s;};
    const std::string left=random(3000),right=random(3000),unit=random(47);
    NativeGraph g;g.labels={left,right};
    std::map<int,uint32_t> allele;
    for(int r:{2,3,6}){allele[r]=uint32_t(g.labels.size());std::string s;for(int i=0;i<r;++i)s+=unit;g.labels.push_back(s);}
    for(auto[r,v]:allele){g.edges.insert(NativeGraph::edge(0,2*v));g.edges.insert(NativeGraph::edge(2*v,2));}
    auto path=[&](const std::string&name,int r){g.paths.push_back({name,{0,2*allele[r],2}});};
    path("REF",3);path("A#1",2);path("A#2",2);path("B#1",3);path("B#2",3);path("C#1",6);path("C#2",6);
    std::string pattern_buffer=test_temp_pattern("dipgenie-copy");char*pattern=pattern_buffer.data();
    const char*created=mkdtemp(pattern);assert(created);const std::filesystem::path dir(created);
    {
        std::string haps[2];for(int h=0;h<2;++h){haps[h]=left;for(int i=0;i<(h?6:3);++i)haps[h]+=unit;haps[h]+=right;}
        std::ofstream reads[2]={std::ofstream(dir/"R1.fq"),std::ofstream(dir/"R2.fq")};
        for(size_t i=0;i<6000;++i) {
            const auto&s=haps[i%2];const size_t start=rng()%(s.size()-300);const auto f=s.substr(start,300);
            const std::string e[2]={f.substr(0,120),reverse_complement(f.substr(180,120))};
            for(int side=0;side<2;++side)reads[side]<<'@'<<i<<'/'<<side+1<<'\n'<<e[side]<<"\n+\n"<<std::string(120,'I')<<'\n';
        }
    }
    for(int threads:{1,4}) {
        const auto parts=CommonSegments::build(g,31,threads);
        auto index=SequenceIndex::build(g,parts,31,25,threads);
        const auto raw=read_evidence(index,(dir/"R1.fq").string(),(dir/"R2.fq").string(),{},true,threads);
        const int H=int(g.paths.size());
        auto prepared=build_model(index,raw,H,int(parts.blocks()),2,threads);
        assert(prepared.copy_strings>0&&prepared.copy_level_targets>0);
        {   // Independent oracle: every exact occurrence in every read end.
            std::map<uint64_t,uint64_t> oracle;
            for(int side=0;side<2;++side){std::ifstream in(dir/(side?"R2.fq":"R1.fq"));std::string h,seq,plus,q;
                while(std::getline(in,h)&&std::getline(in,seq)&&std::getline(in,plus)&&std::getline(in,q))
                    words(seq,31,[&](uint64_t w,size_t,bool){++oracle[w];});}
            for(size_t x=0;x<prepared.target_ids.size();++x)if(prepared.copy_level[x])
                assert(prepared.occurrence_counts[x]==oracle[index.targets[prepared.target_ids[x]].word]);
        }
        const auto result=solve_presence(prepared.model,256,threads);
        assert(!result.path.empty()&&prepared.model.score(result.path)==result.score);
        std::multiset<int> units;
        for(int h=0;h<2;++h) {
            std::vector<size_t> donors;
            for(int z:result.path)donors.push_back(size_t(h?z%H:z/H));
            const auto s=parts.mosaic(donors);
            int n=0;for(size_t p=s.find(unit);p!=std::string::npos;p=s.find(unit,p+unit.size()))++n;
            units.insert(n);
        }
        assert((units==std::multiset<int>{3,6}));
    }
    std::filesystem::remove_all(dir);
    std::cout<<"PASS copy-level posterior identities and raw-read tandem-repeat copy-number recovery (presence-indistinguishable alleles)\n";
}
