#include "ShardedWords.hpp"
#include <cassert>
#include <iostream>
#include <random>

int main() {
    std::mt19937_64 random(621091);
    hic::ShardedWords buffers[3];hic::ShardedWords::Scratch scratch;
    for(size_t shards:{size_t(1),size_t(7),size_t(256)})for(int round=0;round<90;++round) {
        auto owner=[&](uint64_t x){return round%3?size_t((x^(x>>32))%shards):size_t(0);};
        for(auto&buffer:buffers) {
            buffer.clear();std::vector<std::vector<uint64_t>> expected(shards);
            const size_t n=round%9?random()%24000:0;
            for(size_t i=0;i<n;++i) {
                const uint64_t word=i%5?random():i%10?UINT64_MAX:0;
                buffer.push_back(word);expected[owner(word)].push_back(word);
            }
            buffer.finish(shards,owner,scratch);
            for(size_t s=0;s<shards;++s) {
                auto&v=expected[s];std::sort(v.begin(),v.end());
                const auto got=buffer[s];assert(got.size()==v.size());
                for(size_t i=0;i<v.size();++i)assert(got[i]==v[i]);
            }
        }
    }
    bool failed=false;
    try{buffers[0].finish(0,[](uint64_t){return size_t(0);},scratch);}
    catch(const std::invalid_argument&){failed=true;}
    assert(failed);
    std::cout<<"PASS contiguous shards: 810 reused batches; exact multiplicity, empty/skewed shards, boundary words and shared scratch\n";
}
