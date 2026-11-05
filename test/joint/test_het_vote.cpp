#include "BuildModel.hpp"
#include <cassert>
#include <iostream>
#include <random>
using namespace hic;
int main() {
    std::mt19937_64 rng(77);
    for(int trial=0;trial<500;++trial) {
        MoleculeMarkers f;f.id="m"+std::to_string(trial);
        const size_t n=8;std::vector<size_t> remap(n);
        for(size_t i=0;i<n;++i)remap[i]=rng()%5==0?SIZE_MAX:i;
        for(int e=0;e<2;++e)for(size_t i=0;i<1+rng()%4;++i)f.ends[e].push_back(rng()%n);
        for(auto&v:f.ends){std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());}
        std::vector<Score> het(n);for(auto&h:het)h=rng()%3==0?0:Score(rng()%1025);
        std::map<std::pair<size_t,size_t>,Score> uniform,weighted;
        const bool u=molecule_pair_units(f,remap,[&](const auto&p,Score w){uniform[p]+=w;});
        const bool h=molecule_pair_units(f,remap,[&](const auto&p,Score w){weighted[p]+=w;},&het);
        Score total=0,prod=0;
        for(const auto&[p,w]:uniform){total+=w;prod+=het[p.first]*het[p.second];}
        assert(!u||total==1024);
        assert(h==(u&&prod>0));
        if(!h)continue;
        Score sum=0;
        for(const auto&[p,w]:weighted) {
            const Score q=het[p.first]*het[p.second];
            assert(q>0&&uniform.count(p));sum+=w;
            const long double exact=1024.L*q/prod;
            assert(w>=Score(std::floor(exact))&&w<=Score(std::floor(exact))+Score(uniform.size()));
        }
        assert(sum==1024);
        // Determinism.
        std::map<std::pair<size_t,size_t>,Score> again;
        molecule_pair_units(f,remap,[&](const auto&p,Score w){again[p]+=w;},&het);
        assert(again==weighted);
    }
    // A non-informative pair cannot dilute the vote: one het pair takes it all.
    MoleculeMarkers f;f.id="x";f.ends[0]={0,1};f.ends[1]={2,3};
    std::vector<size_t> remap{0,1,2,3};std::vector<Score> het{512,0,1024,0};
    std::map<std::pair<size_t,size_t>,Score> w;
    assert(molecule_pair_units(f,remap,[&](const auto&p,Score v){w[p]+=v;},&het));
    assert(w.size()==1&&w.at({0,2})==1024);
    std::cout<<"PASS heterozygosity-proportional molecule vote: conservation, proportionality, omission, determinism, no dilution\n";
}
