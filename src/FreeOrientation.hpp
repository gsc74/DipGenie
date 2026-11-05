#pragma once
#include "GlobalPhaseMove.hpp"
#include <unordered_map>

namespace hic::free_orientation {
// Exact orientation local search for the unordered-pair charge.
//
// Exchanging the two haplotypes' donors on every segment of an interval
// leaves the unordered recombination charge unchanged, so every interval
// flip is feasible. Candidates come from a segment-level surrogate (linked
// factors whose two endpoints are single-segment heterozygous targets); each
// candidate is accepted only if its EXACT objective change - all genotype and
// phase factors, including multi-segment strings - is positive. The incumbent
// therefore never decreases and is replayed by PresenceModel::score at the end.
struct Report {size_t rounds=0,candidates=0,evaluated=0,accepted=0,flipped_segments=0,continuity_flips=0;
    int exchanges_before=0,exchanges_after=0;Score before=0,after=0;};
// Exchanges realized by a labeled pair sequence: ordered minus unordered charge.
inline int exchanges(const PresenceModel&m,const std::vector<int>&path){
    int e=0;for(size_t b=1;b<path.size();++b)e+=m.ordered_charge(path[b-1],path[b])-m.charge(path[b-1],path[b]);return e;
}

class State {
    const PresenceModel&m;
    std::vector<int>&path;
    std::vector<size_t> occ_begin;std::vector<std::pair<uint32_t,DonorMask>> occ;     // per segment
    std::vector<size_t> adj_begin;std::vector<uint32_t> adj;                          // per target -> contacts
    std::vector<int32_t> c0,c1,d0,d1;std::vector<unsigned char> presence,touched_flag;
    std::vector<uint32_t> touched;
    int swapped(int z)const{return (z%m.H)*m.H+z/m.H;}
    Score target_value(size_t x,unsigned char p)const{
        return phase_reward_scale*m.weights.genotype*PresenceModel::genotype_reward(m.targets[x],p);
    }
    Score contact_value(const PresenceContact&e,unsigned char a,unsigned char b)const{
        return phase_factor(a,b,e.relation,e.panel,m.neutral_phase,m.presence_phase)*e.units*m.weights.contact;
    }
public:
    State(const PresenceModel&model,std::vector<int>&p):m(model),path(p) {
        const size_t T=m.targets.size(),B=size_t(m.blocks);
        occ_begin.assign(B+1,0);
        for(const auto&t:m.targets)for(const auto&o:t.occurrences)++occ_begin[size_t(o.first)+1];
        for(size_t b=0;b<B;++b)occ_begin[b+1]+=occ_begin[b];
        occ.resize(occ_begin[B]);auto fill=occ_begin;
        for(size_t x=0;x<T;++x)for(const auto&o:m.targets[x].occurrences)occ[fill[size_t(o.first)]++]={uint32_t(x),o.second};
        adj_begin.assign(T+1,0);
        for(const auto&e:m.contacts){++adj_begin[e.x+1];++adj_begin[e.y+1];}
        for(size_t x=0;x<T;++x)adj_begin[x+1]+=adj_begin[x];
        adj.resize(adj_begin[T]);auto at=adj_begin;
        for(size_t i=0;i<m.contacts.size();++i){adj[at[m.contacts[i].x]++]=uint32_t(i);adj[at[m.contacts[i].y]++]=uint32_t(i);}
        c0.assign(T,0);c1.assign(T,0);d0.assign(T,0);d1.assign(T,0);presence.assign(T,0);touched_flag.assign(T,0);
        for(size_t b=0;b<B;++b)for(size_t i=occ_begin[b];i<occ_begin[b+1];++i){
            const auto s=m.signature(occ[i].second,path[b]);c0[occ[i].first]+=s&1;c1[occ[i].first]+=s>>1;
        }
        for(size_t x=0;x<T;++x)presence[x]=(c0[x]>0)|((c1[x]>0)<<1);
    }
    unsigned char now(size_t x)const{return presence[x];}
    // Exact objective change of exchanging donors on segments [first,last).
    // commit=true applies it.
    Score flip(int first,int last,bool commit) {
        touched.clear();
        for(int b=first;b<last;++b)for(size_t i=occ_begin[b];i<occ_begin[b+1];++i){
            const uint32_t x=occ[i].first;const auto s=m.signature(occ[i].second,path[b]);
            const int a=s&1,c=s>>1;
            if(!touched_flag[x]){touched_flag[x]=1;touched.push_back(x);}
            d0[x]+=c-a;d1[x]+=a-c;
        }
        auto after=[&](size_t x){return (unsigned char)((c0[x]+d0[x]>0)|((c1[x]+d1[x]>0)<<1));};
        Score delta=0;
        for(uint32_t x:touched){
            const auto p=presence[x],q=after(x);if(p==q)continue;
            delta+=target_value(x,q)-target_value(x,p);
            for(size_t k=adj_begin[x];k<adj_begin[x+1];++k){
                const auto&e=m.contacts[adj[k]];const size_t y=e.x==x?e.y:e.x;
                const bool other=touched_flag[y]&&after(y)!=presence[y];
                if(other&&y<x)continue; // count a doubly changed factor once
                const auto py=presence[y],qy=other?after(y):py;
                const auto [ox,oy,nx,ny]=e.x==x?std::tuple{p,py,q,qy}:std::tuple{py,p,qy,q};
                delta+=contact_value(e,nx,ny)-contact_value(e,ox,oy);
            }
        }
        for(uint32_t x:touched){
            if(commit){c0[x]+=d0[x];c1[x]+=d1[x];presence[x]=(c0[x]>0)|((c1[x]>0)<<1);}
            d0[x]=d1[x]=0;touched_flag[x]=0;
        }
        if(commit)for(int b=first;b<last;++b)path[b]=swapped(path[b]);
        return delta;
    }
};

inline Report improve(const PresenceModel&m,Result&result,int rounds=24,size_t per_round=512) {
    if(!m.free_swaps)throw std::invalid_argument("Orientation local search requires the unordered-pair charge");
    if(result.path.empty()||m.score(result.path)!=result.score)throw std::invalid_argument("orientation incumbent");
    Report report;report.before=report.after=result.score;
    const int recombinations=[&]{int r=0;for(size_t b=1;b<result.path.size();++b)r+=m.charge(result.path[b-1],result.path[b]);return r;}();
    State state(m,result.path);
    const size_t T=m.targets.size();
    std::vector<int> segment(T,-1);
    for(size_t x=0;x<T;++x){const auto&o=m.targets[x].occurrences;if(o.size()==1)segment[x]=o.front().first;}
    for(int round=0;round<rounds;++round) {
        ++report.rounds;
        // Segment-level surrogate on currently heterozygous single-segment endpoints.
        std::unordered_map<uint64_t,Score> pooled;std::vector<int> nodes;
        for(const auto&e:m.contacts){
            const int a=segment[e.x],b=segment[e.y];if(a<0||b<0||a==b)continue;
            const int sx=phase_orientation(state.now(e.x)),sy=phase_orientation(state.now(e.y));if(!sx||!sy)continue;
            const auto lo=std::min(a,b),hi=std::max(a,b);
            pooled[(uint64_t(lo)<<32)|uint32_t(hi)]+=e.units*e.relation*sx*sy;
            nodes.push_back(lo);nodes.push_back(hi);
        }
        std::sort(nodes.begin(),nodes.end());nodes.erase(std::unique(nodes.begin(),nodes.end()),nodes.end());
        const int n=int(nodes.size());if(n<2)break;
        auto id=[&](int s){return int(std::lower_bound(nodes.begin(),nodes.end(),s)-nodes.begin());};
        std::vector<Score> single(n);std::vector<std::vector<std::pair<int,Score>>> left(n);
        for(const auto&[key,w]:pooled){if(!w)continue;const int a=id(int(key>>32)),b=id(int(key&0xffffffffu));
            single[a]-=2*w;single[b]-=2*w;left[b].push_back({a,w});}
        global_phase::Maximum tree(n);
        struct Candidate{Score gain;int lo,hi;};std::vector<Candidate> found;
        for(int b=0;b<n;++b){
            tree.activate(b,1,0,n-1);tree.add(b,single[b],1,0,n-1);
            for(auto[a,w]:left[b])tree.add(a,4*w,1,0,n-1);
            if(tree.value[1]>0)found.push_back({tree.value[1],tree.index[1],b});
        }
        std::sort(found.begin(),found.end(),[](const auto&a,const auto&b){return a.gain!=b.gain?a.gain>b.gain:std::tie(a.lo,a.hi)<std::tie(b.lo,b.hi);});
        found.erase(std::unique(found.begin(),found.end(),[](const auto&a,const auto&b){return a.lo==b.lo&&a.hi==b.hi;}),found.end());
        if(found.size()>per_round)found.resize(per_round);
        report.candidates+=found.size();
        size_t accepted=0;
        for(const auto&c:found){
            const int first=nodes[c.lo],last=nodes[c.hi]+1;
            ++report.evaluated;
            const Score delta=state.flip(first,last,false);
            if(delta>0){state.flip(first,last,true);result.score+=delta;++accepted;report.flipped_segments+=size_t(last-first);}
        }
        report.accepted+=accepted;
        if(!accepted)break;
    }
    // Continuity tie-break (lexicographic secondary criterion). Among labelings
    // of equal-or-better score, prefer the fewest exchanges, i.e. keep each
    // continuing donor on its haplotype as the v1 ordered charge would. Flipping
    // [b,e) between consecutive exchange boundaries removes both exchanges and
    // is accepted only when the exact objective change is >= 0.
    report.exchanges_before=exchanges(m,result.path);
    for(int pass=0;pass<8;++pass) {
        std::vector<int> boundary;
        for(int b=1;b<m.blocks;++b)if(m.ordered_charge(result.path[b-1],result.path[b])>m.charge(result.path[b-1],result.path[b]))boundary.push_back(b);
        if(boundary.empty())break;
        size_t changed=0;
        for(size_t i=0;i<boundary.size();++i) {
            const int b=boundary[i],e=i+1<boundary.size()?boundary[i+1]:m.blocks;
            if(m.ordered_charge(result.path[b-1],result.path[b])==m.charge(result.path[b-1],result.path[b]))continue;
            const Score delta=state.flip(b,e,false);
            if(delta>=0){state.flip(b,e,true);result.score+=delta;++changed;++report.continuity_flips;report.flipped_segments+=size_t(e-b);++i;}
        }
        if(!changed)break;
    }
    report.exchanges_after=exchanges(m,result.path);
    int after=0;for(size_t b=1;b<result.path.size();++b)after+=m.charge(result.path[b-1],result.path[b]);
    if(after!=recombinations)throw std::logic_error("Orientation flip changed the unordered charge");
    if(m.score(result.path)!=result.score||result.score>result.upper)throw std::logic_error("Orientation local search replay");
    report.after=result.score;
    return report;
}
}
