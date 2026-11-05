#pragma once
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include "IntegerSort.hpp"

namespace hic {
// Exact sorted shards in one allocation, instead of one growing vector per
// shard. A worker reuses partition scratch; a task retains its data capacity
// between input batches. No observation is sampled or deduplicated here.
class ShardedWords {
    std::vector<uint64_t> data_;
    std::vector<size_t> offsets_;
public:
    struct Scratch {std::vector<uint64_t> data;std::vector<size_t> next;};
    struct View {
        const std::vector<uint64_t>& data;
        size_t offset,count;
        size_t size()const{return count;}
        uint64_t operator[](size_t i)const{return data[offset+i];}
    };
    void clear(){data_.clear();offsets_.clear();}
    void push_back(uint64_t word){data_.push_back(word);}
    template<class Owner> void finish(size_t shards,const Owner&owner,Scratch&scratch) {
        if(!shards)throw std::invalid_argument("Positive word shard count required");
        offsets_.assign(shards+1,0);
        for(uint64_t word:data_)++offsets_.at(owner(word)+1);
        for(size_t i=1;i<offsets_.size();++i)offsets_[i]+=offsets_[i-1];
        scratch.next.assign(offsets_.begin(),offsets_.end()-1);
        scratch.data.resize(data_.size());
        for(uint64_t word:data_)scratch.data[scratch.next[owner(word)]++]=word;
        data_.swap(scratch.data);
        for(size_t i=0;i<shards;++i)if(offsets_[i+1]>offsets_[i])
            sort_words(data_.data()+offsets_[i],offsets_[i+1]-offsets_[i]);
    }
    View operator[](size_t shard)const {
        return {data_,offsets_.at(shard),offsets_.at(shard+1)-offsets_[shard]};
    }
};
}
