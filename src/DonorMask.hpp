#pragma once
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace hic {
// Inline first word: existing small panels need no per-mask allocation.
// Larger panels retain every donor bit rather than truncating or silently
// selecting a subset. The type stores membership, not path-search history.
class DonorMask {
    uint64_t low_=0;
    std::vector<uint64_t> high_;
public:
    DonorMask(uint64_t low=0):low_(low){}
    void set(size_t bit) {
        if(bit<64){low_|=uint64_t(1)<<bit;return;}
        const size_t word=bit/64-1;
        if(high_.size()<=word)high_.resize(word+1,0);
        high_[word]|=uint64_t(1)<<(bit%64);
    }
    bool test(size_t bit)const {
        if(bit<64)return (low_>>bit)&1;
        const size_t word=bit/64-1;
        return word<high_.size()&&((high_[word]>>(bit%64))&1);
    }
    DonorMask& operator|=(const DonorMask&other) {
        low_|=other.low_;
        if(high_.size()<other.high_.size())high_.resize(other.high_.size(),0);
        for(size_t i=0;i<other.high_.size();++i)high_[i]|=other.high_[i];
        return *this;
    }
    explicit operator bool()const{return low_||!high_.empty();}
    size_t count()const {
        size_t n=__builtin_popcountll(low_);
        for(uint64_t word:high_)n+=__builtin_popcountll(word);
        return n;
    }
    size_t intersection_count(const DonorMask&other)const {
        size_t n=__builtin_popcountll(low_&other.low_);
        for(size_t i=0;i<std::min(high_.size(),other.high_.size());++i)
            n+=__builtin_popcountll(high_[i]&other.high_[i]);
        return n;
    }
    bool within(size_t donors)const {
        if(!donors)return !static_cast<bool>(*this);
        const size_t words=1+(donors-1)/64;
        if(high_.size()+1>words)return false;
        if(donors%64==0||high_.size()+1<words)return true;
        const uint64_t last=high_.empty()?low_:high_.back();
        return (last>>(donors%64))==0;
    }
    static DonorMask full(size_t donors) {
        DonorMask m;
        if(!donors)return m;
        m.low_=donors>=64?UINT64_MAX:(uint64_t(1)<<donors)-1;
        if(donors>64) {
            m.high_.assign((donors-1)/64,UINT64_MAX);
            if(donors%64)m.high_.back()=(uint64_t(1)<<(donors%64))-1;
        }
        return m;
    }
    friend bool operator==(const DonorMask&a,const DonorMask&b) {
        return a.low_==b.low_&&a.high_==b.high_;
    }
    friend bool operator!=(const DonorMask&a,const DonorMask&b){return !(a==b);}
    uint64_t fingerprint()const {
        auto mix=[](uint64_t x){x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
            x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);};
        uint64_t h=mix(low_)^mix(high_.size());
        for(size_t i=0;i<high_.size();++i)h=mix(h^mix(high_[i]+uint64_t(i)));
        return h;
    }
    friend bool operator<(const DonorMask&a,const DonorMask&b) {
        // Numeric ordering preserves the old uint64_t ordering for H<=64.
        if(a.high_.size()!=b.high_.size())return a.high_.size()<b.high_.size();
        for(size_t i=a.high_.size();i>0;--i)
            if(a.high_[i-1]!=b.high_[i-1])return a.high_[i-1]<b.high_[i-1];
        return a.low_<b.low_;
    }
};

inline int checked_pair_domain(size_t donors) {
    if(!donors||donors>size_t(std::numeric_limits<int>::max())/donors)
        throw std::invalid_argument("Joint donor-pair domain exceeds integer index range");
    return int(donors*donors);
}
}
