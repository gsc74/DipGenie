#include "ExactMatch.hpp"
#include <algorithm>
#include <cassert>
#include <string>
#include <random>
#include <iostream>
#include <sys/mman.h>
#include <unistd.h>
using namespace hic::exact_match;
static void check(const char*a,const char*b,size_t n) {
    const auto p=scalar<false>(a,b,n),s=scalar<true>(a,b,n);
    assert(prefix(a,b,n)==p&&suffix(a,b,n)==s);
#if DG_EXACT_X86
    if(__builtin_cpu_supports("avx2")) {
        assert(avx2<false>(a,b,n)==p);assert(avx2<true>(a,b,n)==s);
    }
    if(__builtin_cpu_supports("avx512f")&&__builtin_cpu_supports("avx512bw")) {
        assert(avx512<false>(a,b,n)==p);assert(avx512<true>(a,b,n)==s);
    }
#endif
}
int main() {
    std::mt19937_64 rng(8192);
    for(size_t n=0;n<260;++n)for(size_t p=0;p<=n;++p)for(int shape=0;shape<3;++shape) {
        std::string a(n,'A'),b(n,'A');
        if(p<n) {
            if(shape==0)b[p]='C';
            if(shape==1)a[p]=b[p]='N';
            if(shape==2)a[p]=b[p]=char(128);
        }
        check(a.data(),b.data(),n);
    }
    for(int t=0;t<5000;++t) {
        std::string a(rng()%1025,'A');for(auto&c:a)c="ACGT"[rng()%4];auto b=a;
        if(!b.empty())b[rng()%b.size()]="ACGTN"[rng()%5];
        check(a.data(),b.data(),a.size());
    }
    // Inaccessible pages catch overreads at both ends, including zero-length
    // intervals and every vector/tail boundary, without allocator padding.
    const size_t page=size_t(sysconf(_SC_PAGESIZE));
    auto*p=static_cast<char*>(mmap(nullptr,4*page,PROT_READ|PROT_WRITE,
                                 MAP_PRIVATE|MAP_ANONYMOUS,-1,0));
    assert(p!=MAP_FAILED);
    assert(!mprotect(p,page,PROT_NONE));assert(!mprotect(p+3*page,page,PROT_NONE));
    std::fill_n(p+page,2*page,'A');
    for(size_t n=0;n<260;++n) {
        check(p+page,p+page,n);check(p+3*page-n,p+3*page-n,n);
    }
    assert(!munmap(p,4*page));
    std::cout<<"scalar/AVX2/AVX512 exact extension and guarded boundaries passed\n";
}
