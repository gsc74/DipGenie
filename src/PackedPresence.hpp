#pragma once
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hic {
// Literal presence on two haplotypes has four states, not 256. Packing is
// lossless and does not merge histories or alter donor/phase assignments.
class PackedPresence {
    size_t size_=0;
    std::vector<uint64_t> words_;
public:
    PackedPresence()=default;
    explicit PackedPresence(size_t n):size_(n),words_(n/32+(n%32!=0),0){}
    size_t size()const{return size_;}
    const std::vector<uint64_t>& words()const{return words_;}
    unsigned char operator[](size_t i)const{return (words_[i/32]>>(2*(i%32)))&3;}
    void add(size_t i,unsigned char presence){words_[i/32]|=uint64_t(presence&3)<<(2*(i%32));}
    bool operator==(const PackedPresence&other)const{return size_==other.size_&&words_==other.words_;}
    // Disjoint ranges into a zero-initialized destination, as constructed by
    // FrontierCopy. Source and destination may have different bit alignment.
    void copy_zero_range(const PackedPresence&source,size_t from,size_t to,size_t length) {
        while(length) {
            const size_t offset=to%32,start=from%32,n=std::min(length,32-offset);
            uint64_t value=source.words_[from/32]>>(2*start);
            if(start+n>32)value|=source.words_[from/32+1]<<(2*(32-start));
            if(n<32)value&=(uint64_t(1)<<(2*n))-1;
            words_[to/32]|=value<<(2*offset);
            from+=n;to+=n;length-=n;
        }
    }
    std::vector<unsigned char> unpack()const {
        std::vector<unsigned char> result(size_);
        for(size_t i=0;i<size_;++i)result[i]=(*this)[i];
        return result;
    }
};
}
