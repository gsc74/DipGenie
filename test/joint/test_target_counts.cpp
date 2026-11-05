#include "TestTemp.hpp"
#include "SequenceEvidence.hpp"
#include <cassert>
#include <filesystem>
#include <iostream>
#include <unistd.h>
using namespace hic;
int main(){
    std::string pattern_buffer=test_temp_pattern("dipgenie-target-count");char*pattern=pattern_buffer.data();
    const char*created=mkdtemp(pattern);assert(created);std::filesystem::path dir(created);
    const std::vector<std::string>seq={"AAAAACGTCCCCCC", "GGGGGACGTAAAAA", "ACGTNACGTACGT", "ACGT"};
    SequenceIndex clean;clean.k=3;clean.w=4;
    // An externally selected panel universe, independent of read windows.
    for(const auto&s:seq)words(s,clean.k,[&](uint64_t word,size_t,bool){
        if(clean.ids.count(word))return;
        clean.ids.emplace(word,clean.targets.size());WordTarget t;t.word=word;t.carriers=1;t.occurrences={{0,1}};clean.targets.push_back(t);
    });
    std::vector<uint64_t>oracle(clean.targets.size()),min_oracle(clean.targets.size());
    std::ofstream files[2]={std::ofstream(dir/"R1.fq"),std::ofstream(dir/"R2.fq")};
    for(size_t i=0;i<seq.size();++i)for(int side=0;side<2;++side){
        const auto s=side?reverse_complement(seq[i]):seq[i];
        files[side]<<"@r"<<i<<'/'<<side+1<<'\n'<<s<<"\n+\n"<<std::string(s.size(),'I')<<'\n';
        std::set<uint64_t>seen,selected;
        words(s,clean.k,[&](uint64_t word,size_t,bool){seen.insert(word);});
        minimizers(s,clean.k,clean.w,[&](uint64_t word,size_t){selected.insert(word);});
        for(auto word:seen)++oracle[clean.ids.at(word)];
        for(auto word:selected)++min_oracle[clean.ids.at(word)];
    }
    for(auto&f:files)f.close();
    uint64_t all=0,selected=0;for(size_t x=0;x<oracle.size();++x){all+=oracle[x];selected+=min_oracle[x];}
    assert(all>selected);
    for(int threads:{1,2,4}){
        auto index=clean;const auto evidence=read_evidence(index,(dir/"R1.fq").string(),(dir/"R2.fq").string(),{},true,threads);
        assert(evidence.target_end_observations==all&&evidence.minimizer_target_end_observations==selected);
        for(size_t x=0;x<oracle.size();++x)assert(index.targets[x].count==oracle[x]);
        // Phase markers must remain selected minimizers, not all count hits.
        for(const auto&m:evidence.markers){
            const size_t i=std::stoul(m.id.substr(1));
            for(int side=0;side<2;++side){
                auto s=side?reverse_complement(seq[i]):seq[i];std::set<size_t>ids;
                minimizers(s,clean.k,clean.w,[&](uint64_t w,size_t){ids.insert(clean.ids.at(w));});
                assert(m.ends[side]==std::vector<size_t>(ids.begin(),ids.end()));
            }
        }
    }
    std::filesystem::remove_all(dir);
    std::cout<<"PASS context-independent exact target counts, one vote per end, reverse complements, Ns, short windows, repeated words, thread parity; phase markers unchanged\n";
}
