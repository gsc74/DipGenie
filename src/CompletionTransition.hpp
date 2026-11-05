#pragma once
#include "InferenceTypes.hpp"
#include "DonorMask.hpp"

namespace hic {
// Separable exact max-plus transform used by checkpointed bounds.
struct CompletionTransition {
    struct Entry {Score score;int switches;};
    static bool better(Entry a,Entry b) {
        return a.score!=b.score?a.score>b.score:a.switches<b.switches;
    }
    // Hamming donor-pair transition: two separable max-plus transforms,
    // O(H^2), rather than enumerating H^4 predecessor/destination pairs.
    // free_swaps: unordered-pair charge = min over the two labelings of the
    // destination, i.e. a Hamming transform followed by a zero-cost exchange.
    static void transition(std::vector<Entry>&v,int H,Score lambda,bool free_swaps=false) {
        for(int a=0;a<H;++a) {
            Entry best=v[a*H];
            for(int b=1;b<H;++b)if(better(v[a*H+b],best))best=v[a*H+b];
            best.score-=lambda;++best.switches;
            for(int b=0;b<H;++b)if(better(best,v[a*H+b]))v[a*H+b]=best;
        }
        for(int b=0;b<H;++b) {
            Entry best=v[b];
            for(int a=1;a<H;++a)if(better(v[a*H+b],best))best=v[a*H+b];
            best.score-=lambda;++best.switches;
            for(int a=0;a<H;++a)if(better(best,v[a*H+b]))v[a*H+b]=best;
        }
        if(free_swaps)for(int a=0;a<H;++a)for(int b=a+1;b<H;++b) {
            auto&x=v[a*H+b];auto&y=v[b*H+a];
            if(better(x,y))y=x;else if(better(y,x))x=y;
        }
    }
};
}
