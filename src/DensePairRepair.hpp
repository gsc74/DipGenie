#pragma once
#include "PresenceDP.hpp"
#include <climits>
#include <cmath>

namespace hic::global_phase {
struct DensePairResult {
    std::vector<int> pairs;
    Score score=impossible;
    size_t trace_bytes=0;
};
namespace dense_detail {
struct Best {Score value=impossible;uint32_t pair=UINT32_MAX;};
inline bool better(const Best&a,const Best&b){return a.value>b.value||(a.value==b.value&&a.pair<b.pair);}
struct Two {
    Best first,second;
    void add(Best a){
        if(a.value==impossible)return;
        if(better(a,first)){second=first;first=a;}
        else if(better(a,second))second=a;
    }
};
template<class Trace> DensePairResult solve(const std::vector<std::vector<Score>>&local,
        int H,int budget,int threads,size_t checkpoint_interval,bool free_swaps=false) {
    const size_t B=local.size(),D=size_t(checked_pair_domain(size_t(H)));
    const int R=int(std::min(size_t(budget),2*(B-1)));
    const size_t resources=size_t(R)+1;
    if(D>SIZE_MAX/resources)throw std::overflow_error("complete pair repair trace");
    const size_t stride=D*resources;
    // Checkpointed traceback: scores are kept every s blocks and the
    // predecessor table is rebuilt one interval at a time during traceback.
    // Rebuilt rows equal the forward rows exactly, so the result is identical
    // to a full table while memory is O((B/s) stride Score + s stride Trace).
    size_t s=checkpoint_interval;
    if(!s){
        const double ratio=double(sizeof(Score))/double(sizeof(Trace));
        s=std::max<size_t>(1,size_t(std::ceil(std::sqrt(double(B)*ratio))));
    }
    s=std::min(s,B);
    if(stride>SIZE_MAX/sizeof(Trace)/s||stride>SIZE_MAX/sizeof(Score)/((B+s-1)/s+2))
        throw std::overflow_error("complete pair repair trace");
    const Trace absent=std::numeric_limits<Trace>::max();
    RowExecutor pool(threads);
    // Computes row `block` from row `block-1`; trace may be null.
    auto advance=[&](size_t block,const std::vector<Score>&previous,std::vector<Score>&next,Trace*trace){
        if(free_swaps) {
            // At-most-resource recurrence. A predecessor has charge 0 if its
            // unordered pair agrees, <=1 if it shares a donor, and <=2 always.
            // The inclusive sets permit row/column/global maxima, O(H^2) per
            // resource row, without H^4 scans or an ordered-charge surrogate.
            // A predecessor included with excess allowance is dominated by
            // its entry at the true charge: previous values are monotone in r.
            // Thus the stored predecessor also traces back at its true charge.
            pool.for_each((size_t(R)+16)/16,[&](size_t task){
                std::vector<Best> rows(H),cols(H);
                for(int r=int(task*16);r<=std::min(R,int((task+1)*16-1));++r) {
                    if(r) {
                        std::fill(rows.begin(),rows.end(),Best{});
                        std::fill(cols.begin(),cols.end(),Best{});
                        const auto*prior=previous.data()+size_t(r-1)*D;
                        for(int a=0;a<H;++a)for(int b=0;b<H;++b) {
                            const uint32_t z=uint32_t(a*H+b);const Best v{prior[z],z};
                            if(better(v,rows[a]))rows[a]=v;
                            if(better(v,cols[b]))cols[b]=v;
                        }
                    }
                    Best any;
                    if(r>=2) {
                        const auto*prior=previous.data()+size_t(r-2)*D;
                        for(size_t z=0;z<D;++z) {
                            const Best v{prior[z],uint32_t(z)};
                            if(better(v,any))any=v;
                        }
                    }
                    for(int a=0;a<H;++a)for(int b=0;b<H;++b) {
                        const uint32_t z=uint32_t(a*H+b),rev=uint32_t(b*H+a);
                        Best best{previous[size_t(r)*D+z],z};
                        auto take=[&](Best v){if(better(v,best))best=v;};
                        take({previous[size_t(r)*D+rev],rev});
                        if(r){take(rows[a]);take(rows[b]);take(cols[a]);take(cols[b]);}
                        if(r>=2)take(any);
                        const size_t index=size_t(r)*D+z;
                        const bool valid=best.value!=impossible&&local[block][z]!=impossible;
                        next[index]=valid?best.value+local[block][z]:impossible;
                        if(trace)trace[index]=valid?Trace(best.pair):absent;
                    }
                }
            });
            return;
        }
        const int limit=int(std::min(size_t(R),2*block));
        const int old_limit=int(std::min(size_t(R),2*(block-1)));
        std::fill(next.begin(),next.end(),impossible);
        // Tasks own exact-resource rows. Each scratch array is O(H), not a
        // worker-sized copy of the chromosome or traceback matrix.
        pool.for_each((size_t(limit)+16)/16,[&](size_t task){
            std::vector<Two> rows(H),cols(H),older_rows(H),exclude_col(H);
            for(int r=int(task*16);r<=std::min(limit,int((task+1)*16-1));++r){
                const bool one=r>0&&r-1<=old_limit,two=r>1&&r-2<=old_limit;
                if(one){
                    std::fill(rows.begin(),rows.end(),Two{});std::fill(cols.begin(),cols.end(),Two{});
                    const auto*row=previous.data()+size_t(r-1)*D;
                    for(int a=0;a<H;++a)for(int b=0;b<H;++b){
                        const uint32_t z=uint32_t(a*H+b);const Best value{row[z],z};rows[a].add(value);cols[b].add(value);
                    }
                }
                if(two){
                    std::fill(older_rows.begin(),older_rows.end(),Two{});
                    const auto*row=previous.data()+size_t(r-2)*D;
                    for(int a=0;a<H;++a)for(int b=0;b<H;++b){
                        const uint32_t z=uint32_t(a*H+b);older_rows[a].add({row[z],z});
                    }
                    // For each forbidden second donor b, take each row's
                    // best alternative to b, then retain two distinct rows.
                    std::fill(exclude_col.begin(),exclude_col.end(),Two{});
                    for(int b=0;b<H;++b)for(int a=0;a<H;++a){
                        const auto&t=older_rows[a];const auto value=t.first.pair!=UINT32_MAX&&int(t.first.pair%H)==b?t.second:t.first;
                        exclude_col[b].add(value);
                    }
                }
                for(int a=0;a<H;++a)for(int b=0;b<H;++b){
                    const uint32_t z=uint32_t(a*H+b);Best best;
                    auto take=[&](Best v){if(v.value!=impossible&&better(v,best))best=v;};
                    if(r<=old_limit)take({previous[size_t(r)*D+z],z});
                    if(one){
                        take(rows[a].first.pair==z?rows[a].second:rows[a].first);
                        take(cols[b].first.pair==z?cols[b].second:cols[b].first);
                    }
                    if(two){const auto&t=exclude_col[b];take(t.first.pair!=UINT32_MAX&&int(t.first.pair/H)==a?t.second:t.first);}
                    const size_t index=size_t(r)*D+z;
                    const bool valid=best.value!=impossible&&local[block][z]!=impossible;
                    next[index]=valid?best.value+local[block][z]:impossible;
                    if(trace)trace[index]=valid?Trace(best.pair):absent;
                }
            }
        });
    };
    std::vector<std::vector<Score>> checkpoints;checkpoints.reserve((B+s-1)/s);
    std::vector<Score> previous(stride,impossible),next(stride,impossible);
    std::copy(local[0].begin(),local[0].end(),previous.begin());
    if(free_swaps)for(int r=1;r<=R;++r)
        std::copy(local[0].begin(),local[0].end(),previous.begin()+size_t(r)*D);
    checkpoints.push_back(previous);
    for(size_t block=1;block<B;++block){
        advance(block,previous,next,nullptr);
        previous.swap(next);
        if(block%s==0&&block+1<B)checkpoints.push_back(previous);
    }
    const size_t trace_bytes=checkpoints.size()*stride*sizeof(Score)+s*stride*sizeof(Trace);
    DensePairResult result;result.trace_bytes=trace_bytes;int resource=-1,choice=-1;
    for(int r=0;r<=R;++r)for(size_t z=0;z<D;++z)
        if(previous[size_t(r)*D+z]>result.score){result.score=previous[size_t(r)*D+z];resource=r;choice=int(z);}
    if(choice<0)return result;
    auto charge=[&](int a,int b){const int direct=(a/H!=b/H)+(a%H!=b%H);
        return free_swaps?std::min(direct,(a/H!=b%H)+(a%H!=b/H)):direct;};
    result.pairs.resize(B);
    result.pairs[B-1]=choice;
    std::vector<Trace> trace;
    for(size_t c=checkpoints.size();c--&&B>1;){
        // Interval covers blocks (c*s, min((c+1)*s, B-1)].
        const size_t first=c*s+1,last=std::min((c+1)*s,B-1);
        if(first>last)continue;
        trace.assign((last-first+1)*stride,absent);
        previous=checkpoints[c];
        for(size_t block=first;block<=last;++block){
            advance(block,previous,next,trace.data()+(block-first)*stride);
            previous.swap(next);
        }
        for(size_t block=last;block>=first;--block){
            const auto prior=trace[(block-first)*stride+size_t(resource)*D+size_t(choice)];
            if(prior==absent||size_t(prior)>=D)throw std::logic_error("complete repair missing traceback");
            resource-=charge(int(prior),choice);if(resource<0)throw std::logic_error("complete repair negative resource");
            choice=int(prior);result.pairs[block-1]=choice;
        }
        std::vector<Score>().swap(checkpoints[c]);
    }
    if(resource!=0)throw std::logic_error("complete repair resource replay");
    int used=0;Score score=0;
    for(size_t b=0;b<B;++b){score+=local[b][result.pairs[b]];if(b)used+=charge(result.pairs[b-1],result.pairs[b]);}
    if(used>budget||score!=result.score)throw std::logic_error("complete repair score replay");
    return result;
}
}
// Exact additive joint-pair repair for ordered OR min-swapped donor charge.
// Ordered charge uses exact-resource rows with strict exclusions. Free swaps
// use at-most-resource rows with inclusive unordered/row/column/global maxima.
// Time O(B (R+1) H^2) for each pass; checkpointed traceback uses
// O((B/s+s) (R+1) H^2) auxiliary memory, excluding the local-score input.
inline DensePairResult dense_pair_repair(const std::vector<std::vector<Score>>&local,
        int H,int budget,int threads,size_t checkpoint_interval=0,bool free_swaps=false){
    if(H<1||budget<0||threads<1||local.empty()||local.size()>size_t(INT_MAX/2))
        throw std::invalid_argument("complete pair repair dimensions");
    const size_t D=size_t(checked_pair_domain(size_t(H)));
    for(const auto&row:local)if(row.size()!=D)throw std::invalid_argument("complete pair repair costs");
    if(D<UINT16_MAX)return dense_detail::solve<uint16_t>(local,H,budget,threads,checkpoint_interval,free_swaps);
    return dense_detail::solve<uint32_t>(local,H,budget,threads,checkpoint_interval,free_swaps);
}
}
