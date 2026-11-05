#pragma once
#include <algorithm>
#include <cstdint>
#include <vector>
#include "ksort.h"

namespace hic {
namespace integer_sort_detail {
// Members are implicitly inline: reuse the bundled, allocation-free radix
// implementation without exporting another set of global ksort symbols.
struct Radix {
    static uint64_t key(uint64_t x){return x;}
    KRADIX_SORT_INIT(word, uint64_t, key, 8)
};
}
inline void sort_words(uint64_t*first,size_t count) {
    if(count<2)return;
    // Tiny per-read/shard arrays stay sequential comparison sorts. The radix
    // path runs inside the existing pool, never a nested OpenMP team.
    if(count<4096){std::sort(first,first+count);return;}
    uint64_t low=first[0],high=first[0];bool ordered=true;
    for(size_t i=1;i<count;++i) {
        low=std::min(low,first[i]);high=std::max(high,first[i]);
        ordered&=first[i-1]<=first[i];
    }
    if(ordered)return;
    // Skip identical high bytes: repetitive or narrow-range keys must not
    // pay for eight full passes. Since low!=high, clz is always defined.
    const unsigned shift=((63u-unsigned(__builtin_clzll(low^high)))/8u)*8u;
    integer_sort_detail::Radix{}.rs_sort_word(first,first+count,8,int(shift));
}
inline void sort_words(std::vector<uint64_t>&v){sort_words(v.data(),v.size());}
}
