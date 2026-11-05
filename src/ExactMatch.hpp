#pragma once
#include <cstddef>
#include <cstdint>
#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#define DG_EXACT_X86 1
#else
#define DG_EXACT_X86 0
#endif

namespace hic::exact_match {
inline bool valid(char c){return c=='A'||c=='C'||c=='G'||c=='T';}

template<bool Backward> inline size_t scalar(const char*a,const char*b,size_t n) {
    size_t i=0;
    for(;i<n;++i) {
        const size_t p=Backward?n-1-i:i;
        if(!valid(a[p])||a[p]!=b[p])break;
    }
    return i;
}

#if DG_EXACT_X86
template<bool Backward> __attribute__((target("avx2")))
inline size_t avx2(const char*a,const char*b,size_t n) {
    size_t done=0;
    while(n-done>=32) {
        const size_t p=Backward?n-done-32:done;
        const auto x=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(a+p));
        const auto y=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(b+p));
        auto allowed=_mm256_or_si256(_mm256_cmpeq_epi8(x,_mm256_set1_epi8('A')),
                                    _mm256_cmpeq_epi8(x,_mm256_set1_epi8('C')));
        allowed=_mm256_or_si256(allowed,_mm256_or_si256(
            _mm256_cmpeq_epi8(x,_mm256_set1_epi8('G')),
            _mm256_cmpeq_epi8(x,_mm256_set1_epi8('T'))));
        const uint32_t bad=~uint32_t(_mm256_movemask_epi8(
            _mm256_and_si256(allowed,_mm256_cmpeq_epi8(x,y))));
        if(bad)return done+(Backward?size_t(__builtin_clz(bad)):size_t(__builtin_ctz(bad)));
        done+=32;
    }
    return done+(Backward?scalar<true>(a,b,n-done):scalar<false>(a+done,b+done,n-done));
}

template<bool Backward> __attribute__((target("avx512f,avx512bw")))
inline size_t avx512(const char*a,const char*b,size_t n) {
    size_t done=0;
    while(n-done>=64) {
        const size_t p=Backward?n-done-64:done;
        const auto x=_mm512_loadu_si512(a+p),y=_mm512_loadu_si512(b+p);
        const uint64_t allowed=_mm512_cmpeq_epi8_mask(x,_mm512_set1_epi8('A'))|
            _mm512_cmpeq_epi8_mask(x,_mm512_set1_epi8('C'))|
            _mm512_cmpeq_epi8_mask(x,_mm512_set1_epi8('G'))|
            _mm512_cmpeq_epi8_mask(x,_mm512_set1_epi8('T'));
        const uint64_t bad=~(allowed&_mm512_cmpeq_epi8_mask(x,y));
        if(bad)return done+(Backward?size_t(__builtin_clzll(bad)):size_t(__builtin_ctzll(bad)));
        done+=64;
    }
    return done+(Backward?scalar<true>(a,b,n-done):scalar<false>(a+done,b+done,n-done));
}
#endif

using Function=size_t(*)(const char*,const char*,size_t);
template<bool Backward> inline Function dispatch() {
#if DG_EXACT_X86
    if(__builtin_cpu_supports("avx512f")&&__builtin_cpu_supports("avx512bw"))return avx512<Backward>;
    if(__builtin_cpu_supports("avx2"))return avx2<Backward>;
#endif
    return scalar<Backward>;
}

// NativeLocalizer supplies normalized uppercase strings. Ambiguous bases
// stop extension even when both bytes match. Vector loads stay inside [0,n).
// Runtime dispatch preserves the portable binary's scalar fallback.
inline size_t prefix(const char*a,const char*b,size_t n) {
    static const auto fn=dispatch<false>();return fn(a,b,n);
}
inline size_t suffix(const char*a,const char*b,size_t n) {
    static const auto fn=dispatch<true>();return fn(a,b,n);
}
}
