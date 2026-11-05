#pragma once
#include "SequenceEvidence.hpp"
#include "ControlCalibration.hpp"
#include "ReadLibraries.hpp"
#include <array>
#include <cmath>
#include <map>
#include <numeric>

namespace hic {
// Local read-consistency panel preselection (reads + panel only).
//
// Items are canonical k/w minimizer strings of the complete input panel that
// occur once per carrying path, in one projected reference window, and are
// absent from at least one path. Each path position is projected onto the
// reference by the latest reference-path vertex visited by that path; cyclic
// paths and non-reference vertices therefore need no linear coordinate.
// Exact read-end counts (once per end) give, for every item, the
// log-likelihood of its count under dosage d in {0,1,2} (error rate for d=0,
// fitted conserved-control negative-binomial rates for d=1,2). A haplotype
// pair (a,b) costs the summed negative log-likelihood, relative to each
// item's best dosage, of the dosage d=1[a carries]+1[b carries] it implies.
// Every item contributes, so the selection uses weak evidence at low depth
// instead of discarding it. Complete samples are chosen greedily to minimise
// the sum over windows of the best retained pair cost (REF always retained),
// i.e. a monotone facility-location objective on local diploid read
// likelihood. This is panel preselection only: it makes
// no claim about the diploid objective or a variant-accuracy bound.
struct WindowSelectionOptions {
    uint64_t window_bases=25000;
    size_t cap=51;
    long double scale=16; // cost units per nat; per-item cost is capped at 255 units
};
struct WindowSelectionSummary {
    size_t items=0,windows=0,controls=0,ambiguous=0,repeated=0,unclassified=0;
    long double mean=0,shape=0;
    long initial_cost=0,final_cost=0;
    std::vector<std::pair<std::string,long>> steps;
};
namespace window_select {
inline long double nb_log(uint64_t c,long double mean,long double shape) {
    const long double x=static_cast<long double>(c);
    if(std::isinf(shape))return x*std::log(mean)-mean-std::lgamma(x+1);
    return std::lgamma(x+shape)-std::lgamma(shape)-std::lgamma(x+1)+shape*std::log(shape/(shape+mean))+x*std::log(mean/(shape+mean));
}
// Per-dosage cost in integer units: round(scale*(max_k l_k - l_d)), capped at 255.
inline std::array<uint8_t,3> penalties(uint64_t c,long double mean,long double shape,long double scale) {
    const long double error=std::max(0.5L,mean/100);
    const long double l[3]={nb_log(c,error,std::numeric_limits<long double>::infinity()),nb_log(c,mean/2,shape),nb_log(c,mean,shape)};
    const long double top=std::max({l[0],l[1],l[2]});
    std::array<uint8_t,3> p{};
    for(int d=0;d<3;++d)p[d]=uint8_t(std::min<long double>(255,std::llround(scale*(top-l[d]))));
    return p;
}
// Vertex starts on one path with their projected reference coordinate.
struct Projected {uint64_t offset,coordinate;bool reference;};
inline std::vector<Projected> projection(const NativeGraph&g,size_t h,
        const std::unordered_map<uint32_t,uint64_t>&reference_start) {
    std::vector<Projected> v;v.reserve(g.paths[h].walk.size());
    uint64_t offset=0,last=0;
    for(uint32_t vertex:g.paths[h].walk) {
        const uint64_t n=g.labels[vertex/2].size();
        auto it=reference_start.find(vertex/2);
        if(it!=reference_start.end()){v.push_back({offset,it->second,true});last=it->second+n;}
        else v.push_back({offset,last,false});
        offset+=n;
    }
    return v;
}
inline uint64_t project(const std::vector<Projected>&v,uint64_t p) {
    auto it=std::upper_bound(v.begin(),v.end(),p,[](uint64_t x,const Projected&e){return x<e.offset;});
    if(it==v.begin())return 0;
    --it;
    return it->reference?it->coordinate+(p-it->offset):it->coordinate;
}
}
inline WindowSelectionSummary window_select_panel(NativeGraph&g,const std::vector<ReadLibrary>&libraries,
        int threads,const WindowSelectionOptions&o,std::vector<size_t>*chosen_out=nullptr,
        std::vector<long double>*scores_out=nullptr) {
    WindowSelectionSummary s;
    const size_t H=g.paths.size();
    if(H<=o.cap)throw std::logic_error("window selection requires more paths than the cap");
    if(o.cap<3||o.cap%2!=1||H%2!=1||!o.window_bases)throw std::runtime_error("Panel selection requires REF plus complete sample pairs");
    const size_t B=(H+63)/64;
    RowExecutor pool(threads);
    std::unordered_map<uint32_t,uint64_t> reference_start;
    {uint64_t p=0;for(uint32_t v:g.paths[0].walk){reference_start.emplace(v/2,p);p+=g.labels[v/2].size();}}
    // Sharded identity table: word -> entry. window UINT32_MAX = conflicting.
    struct Entry {uint32_t index=0;};
    constexpr size_t shards=256;
    auto owner=[](uint64_t x){x^=x>>31;x*=0x9e3779b97f4a7c15ULL;return size_t(x>>56);};
    std::vector<FlatMap<uint64_t,Entry>> table(shards);
    std::vector<std::vector<uint64_t>> carriers(shards); // B words per entry
    struct WindowEntry {uint32_t window=0;bool repeated=false;};
    std::vector<std::vector<WindowEntry>> window_entry(shards);
    const size_t batch=std::max<size_t>(1,size_t(threads));
    // Pass 1: the dictionary is the union of every path's k/w minimizers.
    for(size_t first=0;first<H;first+=batch) {
        const size_t last=std::min(H,first+batch);
        std::vector<std::vector<std::vector<uint64_t>>> per(last-first,std::vector<std::vector<uint64_t>>(shards));
        pool.for_each(last-first,[&](size_t i) {
            const auto sequence=g.spell(first+i);
            minimizers(sequence,31,25,[&](uint64_t word,size_t){per[i][owner(word)].push_back(word);});
        });
        pool.for_each(shards,[&](size_t sh) {
            auto&map=table[sh];
            for(const auto&path:per)for(uint64_t word:path[sh]) {
                auto [it,inserted]=map.try_emplace(word);
                if(inserted)it->second.index=uint32_t(map.size()-1);
            }
        });
    }
    for(size_t sh=0;sh<shards;++sh){carriers[sh].assign(table[sh].size()*B,0);window_entry[sh].assign(table[sh].size(),{});}
    // Pass 2: carriers, copy number and projected window come from EVERY exact
    // occurrence in each path, not from that path's own minimizer selection.
    // A neighbouring variant can suppress selection without removing the string.
    std::vector<std::vector<char>> seen(shards);
    for(size_t sh=0;sh<shards;++sh)seen[sh].assign(table[sh].size(),0);
    for(size_t first=0;first<H;first+=batch) {
        const size_t last=std::min(H,first+batch);
        struct Hit {uint32_t index;uint32_t window;bool repeated;};
        std::vector<std::vector<std::vector<Hit>>> per(last-first,std::vector<std::vector<Hit>>(shards));
        pool.for_each(last-first,[&](size_t i) {
            const size_t h=first+i;const auto sequence=g.spell(h);
            const auto proj=window_select::projection(g,h,reference_start);
            std::vector<std::vector<std::pair<uint32_t,uint64_t>>> occ(shards); // index, position
            words(sequence,31,[&](uint64_t word,size_t p,bool){
                const size_t sh=owner(word);auto it=table[sh].find(word);
                if(it!=table[sh].end())occ[sh].emplace_back(it->second.index,p);
            });
            for(size_t sh=0;sh<shards;++sh) {
                auto&v=occ[sh];std::sort(v.begin(),v.end());
                for(size_t a=0;a<v.size();) {
                    size_t b=a+1;while(b<v.size()&&v[b].first==v[a].first)++b;
                    const uint64_t w=window_select::project(proj,v[a].second)/o.window_bases;
                    per[i][sh].push_back({v[a].first,uint32_t(std::min<uint64_t>(w,UINT32_MAX-1)),b-a>1});
                    a=b;
                }
            }
        });
        pool.for_each(shards,[&](size_t sh) {
            for(size_t i=0;i<per.size();++i)for(const auto&hit:per[i][sh]) {
                auto&e=window_entry[sh][hit.index];
                if(!seen[sh][hit.index]){seen[sh][hit.index]=1;e.window=hit.window;}
                else if(e.window!=hit.window)e.window=UINT32_MAX;
                e.repeated|=hit.repeated;
                const size_t h=first+i;carriers[sh][size_t(hit.index)*B+h/64]|=uint64_t(1)<<(h%64);
            }
        });
    }
    // Exact once-per-end counts for every tabled word. SR supplies dosage when
    // present, as in the inference dosage policy; otherwise all libraries.
    bool has_sr=false;for(const auto&l:libraries)has_sr|=l.kind=="SR";
    std::vector<std::vector<uint32_t>> counts(shards);
    for(size_t sh=0;sh<shards;++sh)counts[sh].assign(table[sh].size(),0);
    for(const auto&library:libraries) {
        if(has_sr&&library.kind!="SR")continue;
        FastqInput input[2]={FastqInput(library.r1),FastqInput(library.r2)};
        for(;;) {
            std::vector<FastqRead> reads[2];
            pool.for_each(2,[&](size_t side){FastqRead r;reads[side].reserve(16384);while(reads[side].size()<16384&&input[side].next(r))reads[side].push_back(std::move(r));});
            if(reads[0].size()!=reads[1].size())throw std::runtime_error("Window selection mate count mismatch");
            if(reads[0].empty())break;
            const size_t n=reads[0].size(),tasks=(n+255)/256;
            std::vector<std::vector<std::vector<uint32_t>>> found(tasks,std::vector<std::vector<uint32_t>>(shards));
            pool.for_each(tasks,[&](size_t t) {
                for(size_t i=t*256;i<std::min(n,(t+1)*256);++i)for(int side=0;side<2;++side) {
                    std::vector<std::pair<size_t,uint32_t>> hit;
                    words(reads[side][i].sequence,31,[&](uint64_t w,size_t,bool){
                        const size_t sh=owner(w);auto it=table[sh].find(w);
                        if(it!=table[sh].end())hit.emplace_back(sh,it->second.index);
                    });
                    std::sort(hit.begin(),hit.end());hit.erase(std::unique(hit.begin(),hit.end()),hit.end());
                    for(auto[sh,x]:hit)found[t][sh].push_back(x);
                }
            });
            pool.for_each(shards,[&](size_t sh){for(const auto&t:found)for(uint32_t x:t[sh])++counts[sh][x];});
        }
    }
    auto full=[&](const uint64_t*bits) {
        for(size_t b=0;b<B;++b) {
            const uint64_t want=b+1<B||H%64==0?~uint64_t(0):(uint64_t(1)<<(H%64))-1;
            if(bits[b]!=want)return false;
        }
        return true;
    };
    // Conserved single-copy controls. Read-repeat outliers (strings unique in
    // the panel but collapsed repeats in the sample) are excluded by a fixed
    // 3x-median cap before the moment fit; items above it are unclassified.
    std::vector<uint32_t> control_counts;
    for(size_t sh=0;sh<shards;++sh)for(const auto&[w,e]:table[sh])
        if(seen[sh][e.index]&&!window_entry[sh][e.index].repeated&&window_entry[sh][e.index].window!=UINT32_MAX&&full(&carriers[sh][size_t(e.index)*B]))control_counts.push_back(counts[sh][e.index]);
    if(control_counts.empty())throw std::runtime_error("No conserved window-selection controls");
    std::nth_element(control_counts.begin(),control_counts.begin()+control_counts.size()/2,control_counts.end());
    const uint64_t ceiling=3*uint64_t(std::max<uint32_t>(1,control_counts[control_counts.size()/2]));
    hic_calibration::Controls controls;
    for(uint32_t c:control_counts)if(c<=ceiling)controls.add(c);
    const auto fit=hic_calibration::Calibration::fit(controls);
    s.controls=controls.n;s.mean=fit.mean;s.shape=fit.shape;
    // Items sorted by window, then word for determinism.
    struct Item {uint32_t window;uint64_t word;size_t shard;uint32_t index;std::array<uint8_t,3> pen;};
    std::vector<Item> items;
    for(size_t sh=0;sh<shards;++sh)for(const auto&[w,e]:table[sh]) {
        const uint64_t*bits=&carriers[sh][size_t(e.index)*B];
        const auto&we=window_entry[sh][e.index];
        if(full(bits))continue;
        if(we.repeated){++s.repeated;continue;}
        if(we.window==UINT32_MAX){++s.ambiguous;continue;}
        if(counts[sh][e.index]>ceiling){++s.unclassified;continue;}
        const auto pen=window_select::penalties(counts[sh][e.index],fit.mean,fit.shape,o.scale);
        if(pen[0]==pen[1]&&pen[1]==pen[2]){++s.unclassified;continue;} // no dosage information
        items.push_back({we.window,w,sh,e.index,pen});
    }
    std::sort(items.begin(),items.end(),[](const Item&a,const Item&b){return std::tie(a.window,a.word)<std::tie(b.window,b.word);});
    s.items=items.size();
    // Each window starts on a fresh 64-bit word.
    std::vector<size_t> window_begin;size_t bitpos=0;std::vector<size_t> position(items.size());
    for(size_t i=0;i<items.size();++i) {
        if(!i||items[i].window!=items[i-1].window){bitpos=(bitpos+63)/64*64;window_begin.push_back(bitpos/64);}
        position[i]=bitpos++;
    }
    const size_t words_total=(bitpos+63)/64;window_begin.push_back(words_total);
    const size_t W=window_begin.size()-1;s.windows=W;
    // Bit-plane b of P[d] marks items whose dosage-d cost has bit b set.
    std::vector<uint64_t> P[3][8];for(auto&d:P)for(auto&v:d)v.assign(words_total,0);
    std::vector<std::vector<uint64_t>> D(H,std::vector<uint64_t>(words_total));
    pool.for_each(H,[&](size_t h){
        for(size_t i=0;i<items.size();++i)
            if(carriers[items[i].shard][size_t(items[i].index)*B+h/64]>>(h%64)&1)D[h][position[i]/64]|=uint64_t(1)<<(position[i]%64);
    });
    for(size_t i=0;i<items.size();++i)for(int d=0;d<3;++d)for(int b=0;b<8;++b)
        if(items[i].pen[d]>>b&1)P[d][b][position[i]/64]|=uint64_t(1)<<(position[i]%64);
    table.clear();carriers.clear();counts.clear();window_entry.clear();seen.clear();
    auto cost=[&](size_t a,size_t b,size_t w) {
        long c=0;const auto&A=D[a];const auto&Bv=D[b];
        for(size_t x=window_begin[w];x<window_begin[w+1];++x) {
            const uint64_t M[3]={~A[x]&~Bv[x],A[x]^Bv[x],A[x]&Bv[x]};
            for(int d=0;d<3;++d)for(int b=0;b<8;++b)c+=long(__builtin_popcountll(M[d]&P[d][b][x]))<<b;
        }
        return c;
    };
    auto sample=[](const std::string&name){auto p=name.find('#');return p==std::string::npos?name:name.substr(0,p);};
    std::map<std::string,std::vector<size_t>> groups;for(size_t h=1;h<H;++h)groups[sample(g.paths[h].name)].push_back(h);
    std::vector<std::string> names;std::vector<std::vector<size_t>> members;
    for(auto&[n,v]:groups){if(v.size()!=2)throw std::runtime_error("Cannot resolve two panel haplotypes for "+n);names.push_back(n);members.push_back(v);}
    std::vector<size_t> S{0};std::vector<long> best(W);
    pool.for_each(W,[&](size_t w){best[w]=cost(0,0,w);});
    s.initial_cost=std::accumulate(best.begin(),best.end(),0L);
    std::vector<char> taken(names.size());
    std::vector<long double> score(H,0);
    for(size_t step=0;step<(o.cap-1)/2;++step) {
        std::vector<long> gain(names.size(),-1);
        pool.for_each(names.size(),[&](size_t c) {
            if(taken[c])return;
            long total=0;
            for(size_t w=0;w<W;++w) {
                long b=best[w];
                for(size_t y:members[c]){for(size_t x:S)b=std::min(b,cost(x,y,w));for(size_t x:members[c])b=std::min(b,cost(x,y,w));}
                total+=best[w]-b;
            }
            gain[c]=total;
        });
        // Exact ties resolve by sample name order.
        const size_t c=size_t(std::max_element(gain.begin(),gain.end())-gain.begin());
        taken[c]=1;
        pool.for_each(W,[&](size_t w){for(size_t y:members[c]){for(size_t x:S)best[w]=std::min(best[w],cost(x,y,w));for(size_t x:members[c])best[w]=std::min(best[w],cost(x,y,w));}});
        for(size_t y:members[c]){S.push_back(y);score[y]=static_cast<long double>(gain[c]);}
        s.steps.emplace_back(names[c],gain[c]);
    }
    s.final_cost=std::accumulate(best.begin(),best.end(),0L);
    std::sort(S.begin(),S.end());
    if(chosen_out)*chosen_out=S;
    if(scores_out)*scores_out=score;
    return s;
}
}
