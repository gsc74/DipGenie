#pragma once
#include "PresenceDP.hpp"
#include "DensePairRepair.hpp"
#include <numeric>
#include <queue>

namespace hic::global_phase {
// A large-neighbourhood move for the SAME objective. The temporary phase
// proposal may exceed R; only a repaired, fully rescored feasible pair can
// replace the incumbent. No truth, archived path or external phaser is used.
struct Block { int begin=0,end=0,current=0; std::vector<int> choices; };
struct Edge { int a,b; Score weight; };
struct Repair { std::vector<int> choices; Score score=impossible; size_t trace_bytes=0; };

// Exact additive resource DP over the proposed block/donor alternatives.
// One byte stores each predecessor; its resource index follows from the
// actual donor transition. Current-pair choices guarantee a feasible path.
inline Repair repair(const PresenceModel&m,const std::vector<Block>&blocks,
                     const std::vector<std::vector<Score>>&local,int threads) {
    if(blocks.empty()||local.size()!=blocks.size())throw std::invalid_argument("phase repair dimensions");
    bool complete=true;
    for(const auto&block:blocks){
        if(block.choices.size()!=size_t(m.D())){complete=false;break;}
        for(size_t i=0;i<block.choices.size();++i)if(block.choices[i]!=int(i)){complete=false;break;}
        if(!complete)break;
    }
    if(complete){auto result=dense_pair_repair(local,m.H,m.budget,threads,0,m.free_swaps&&!m.ordered_repair);return{std::move(result.pairs),result.score,result.trace_bytes};}
    const int R=int(std::min<size_t>(size_t(m.budget),2*(blocks.size()-1)));
    size_t widest=0,bytes=0;
    std::vector<size_t> offset(blocks.size());
    for(size_t b=0;b<blocks.size();++b) {
        const size_t n=blocks[b].choices.size();
        if(!n||n>=255||local[b].size()!=n)throw std::invalid_argument("phase repair choices");
        offset[b]=bytes;
        if(n>SIZE_MAX/size_t(R+1)||bytes>SIZE_MAX-n*size_t(R+1))throw std::overflow_error("phase repair trace");
        bytes+=n*size_t(R+1);widest=std::max(widest,n);
    }
    std::vector<unsigned char> trace(bytes,255);
    std::vector<Score> previous(size_t(R+1)*widest,impossible),next(previous.size(),impossible);
    for(size_t j=0;j<blocks[0].choices.size();++j)previous[j]=local[0][j];
    RowExecutor pool(threads);
    for(size_t b=1;b<blocks.size();++b) {
        const auto&from=blocks[b-1].choices;const auto&to=blocks[b].choices;
        std::vector<unsigned char> cost(from.size()*to.size());
        for(size_t j=0;j<to.size();++j)for(size_t i=0;i<from.size();++i)
            cost[j*from.size()+i]=m.charge(from[i],to[j]);
        const int limit=int(std::min<size_t>(R,2*b)),old_limit=int(std::min<size_t>(R,2*(b-1)));
        pool.for_each((size_t(limit)+64)/64,[&](size_t task) {
            for(int r=int(task*64);r<=std::min(limit,int((task+1)*64-1));++r)
              for(size_t j=0;j<to.size();++j) {
                Score best=impossible;unsigned char predecessor=255;
                for(size_t i=0;i<from.size();++i) {
                    const int q=r-cost[j*from.size()+i];if(q<0||q>old_limit)continue;
                    const Score v=previous[size_t(q)*widest+i];
                    if(v>best){best=v;predecessor=static_cast<unsigned char>(i);}
                }
                next[size_t(r)*widest+j]=best==impossible?impossible:best+local[b][j];
                trace[offset[b]+size_t(r)*to.size()+j]=predecessor;
              }
        });
        previous.swap(next);
    }
    Repair out;out.trace_bytes=bytes;int resource=-1,choice=-1;
    for(int r=0;r<=R;++r)for(size_t j=0;j<blocks.back().choices.size();++j)
        if(previous[size_t(r)*widest+j]>out.score){out.score=previous[size_t(r)*widest+j];resource=r;choice=int(j);}
    if(choice<0)return out;
    out.choices.resize(blocks.size());
    for(size_t b=blocks.size();b--;) {
        out.choices[b]=choice;if(!b)break;
        const auto p=trace[offset[b]+size_t(resource)*blocks[b].choices.size()+size_t(choice)];
        if(p==255||p>=blocks[b-1].choices.size())throw std::logic_error("phase repair missing trace");
        resource-=m.charge(blocks[b-1].choices[p],blocks[b].choices[choice]);choice=p;
    }
    if(resource!=0)throw std::logic_error("phase repair resource replay");
    Score replay=0;int used=0;
    for(size_t b=0;b<blocks.size();++b) {
        replay+=local[b][out.choices[b]];
        if(b)used+=m.charge(blocks[b-1].choices[out.choices[b-1]],blocks[b].choices[out.choices[b]]);
    }
    if(replay!=out.score||used>m.budget)throw std::logic_error("phase repair score replay");
    return out;
}

inline Score energy(const std::vector<Edge>&edges,const std::vector<int>&s) {
    Score value=0;for(const auto&e:edges)value+=e.weight*s[e.a]*s[e.b];return value;
}
inline std::vector<int> coordinate(const std::vector<Edge>&edges,std::vector<int>s) {
    std::vector<std::vector<std::pair<int,Score>>> adjacent(s.size());
    std::vector<Score> gain(s.size());
    for(const auto&e:edges) {
        adjacent[e.a].push_back({e.b,e.weight});adjacent[e.b].push_back({e.a,e.weight});
        gain[e.a]-=2*e.weight*s[e.a]*s[e.b];gain[e.b]-=2*e.weight*s[e.a]*s[e.b];
    }
    std::priority_queue<std::pair<Score,int>> queue;
    for(size_t i=0;i<s.size();++i)queue.push({gain[i],-int(i)});
    while(!queue.empty()) {
        auto[g,negative]=queue.top();queue.pop();const int a=-negative;
        if(g!=gain[a])continue;
        if(g<=0)break;
        s[a]=-s[a];gain[a]=-gain[a];queue.push({gain[a],-a});
        for(auto[b,w]:adjacent[a]){gain[b]-=4*w*s[a]*s[b];queue.push({gain[b],-b});}
    }
    return s;
}
// Range-add / maximum tree: find the best whole interval flip in
// O((blocks+contacts) log blocks), rather than enumerating all intervals.
struct Maximum {
    int n;std::vector<Score> value,lazy;std::vector<int> index;
    explicit Maximum(int size):n(size),value(4*size,impossible),lazy(4*size),index(4*size,-1){}
    void apply(int v,Score d){if(value[v]!=impossible)value[v]+=d;lazy[v]+=d;}
    void push(int v){if(lazy[v]){apply(v*2,lazy[v]);apply(v*2+1,lazy[v]);lazy[v]=0;}}
    void pull(int v){int c=(value[v*2]>=value[v*2+1])?v*2:v*2+1;value[v]=value[c];index[v]=index[c];}
    void activate(int at,int v,int l,int r){if(l==r){value[v]=0;lazy[v]=0;index[v]=at;return;}push(v);int mid=(l+r)/2;
        if(at<=mid)activate(at,v*2,l,mid);else activate(at,v*2+1,mid+1,r);pull(v);}
    void add(int end,Score d,int v,int l,int r){if(end>=r){apply(v,d);return;}push(v);int mid=(l+r)/2;
        add(end,d,v*2,l,mid);if(end>mid)add(end,d,v*2+1,mid+1,r);pull(v);}
};
inline std::vector<int> intervals(const std::vector<Edge>&edges,std::vector<int>s) {
    const int n=int(s.size());if(!n)return s;
    std::vector<std::vector<std::pair<int,Score>>> left(n);
    for(const auto&e:edges)left[e.b].push_back({e.a,e.weight});
    // A finite proposal effort; acceptance never depends on convergence.
    for(int iteration=0;iteration<32;++iteration) {
        std::vector<Score> single(n);
        for(const auto&e:edges){Score w=-2*e.weight*s[e.a]*s[e.b];single[e.a]+=w;single[e.b]+=w;}
        Maximum tree(n);Score best=0;int lo=-1,hi=-1;
        for(int b=0;b<n;++b) {
            tree.activate(b,1,0,n-1);tree.add(b,single[b],1,0,n-1);
            for(auto[a,w]:left[b])tree.add(a,4*w*s[a]*s[b],1,0,n-1);
            if(tree.value[1]>best){best=tree.value[1];lo=tree.index[1];hi=b;}
        }
        if(lo<0)break;
        const auto before=energy(edges,s);for(int b=lo;b<=hi;++b)s[b]=-s[b];
        if(energy(edges,s)!=before+best)throw std::logic_error("global interval replay");
    }
    return s;
}
inline std::vector<int> orient(const std::vector<Edge>&edges,int n) {
    std::vector<size_t> order(edges.size());std::iota(order.begin(),order.end(),0);
    auto magnitude=[](Score x){return x<0?-x:x;};
    std::sort(order.begin(),order.end(),[&](size_t a,size_t b){const auto x=magnitude(edges[a].weight),y=magnitude(edges[b].weight);return x!=y?x>y:a<b;});
    std::vector<int> parent(n),rank(n);std::iota(parent.begin(),parent.end(),0);
    auto root=[&](int a){while(a!=parent[a]){parent[a]=parent[parent[a]];a=parent[a];}return a;};
    std::vector<std::vector<std::pair<int,int>>> tree(n);
    for(size_t i:order) {
        const auto&e=edges[i];int a=root(e.a),b=root(e.b);if(a==b)continue;
        if(rank[a]<rank[b])std::swap(a,b);
        parent[b]=a;
        if(rank[a]==rank[b])++rank[a];
        const int sign=e.weight>0?1:-1;tree[e.a].push_back({e.b,sign});tree[e.b].push_back({e.a,sign});
    }
    std::vector<int> s(n);
    for(int i=0;i<n;++i)if(!s[i]) {
        s[i]=1;std::vector<int> stack{i};
        while(!stack.empty()){const int a=stack.back();stack.pop_back();for(auto[b,w]:tree[a])if(!s[b]){s[b]=w*s[a];stack.push_back(b);}}
    }
    s=intervals(edges,coordinate(edges,std::move(s)));
    auto other=intervals(edges,coordinate(edges,std::vector<int>(n,1)));
    return energy(edges,other)>energy(edges,s)?other:s;
}

struct Report {size_t span=0,blocks=0,edges=0,trace_bytes=0;int proposed_recombinations=0,repaired_recombinations=0;Score before=0,unconstrained=0,repaired=0,after=0;bool accepted=false;};
inline Report move(const PresenceModel&m,Result&result,size_t span,int threads) {
    if(!span||result.path.empty()||m.score(result.path)!=result.score)throw std::invalid_argument("global move incumbent");
    Report report;report.span=span;report.before=result.score;report.after=result.score;
    std::vector<Block> blocks;std::vector<int> segment_block(m.blocks);
    for(int b=0;b<m.blocks;++b) {
        if(blocks.empty()||blocks.back().current!=result.path[b]||size_t(b-blocks.back().begin)>=span)
            blocks.push_back({b,b+1,result.path[b],{}});
        blocks.back().end=b+1;segment_block[b]=int(blocks.size()-1);
    }
    report.blocks=blocks.size();
    const size_t T=m.targets.size();
    std::vector<int> owner(T,-1);std::vector<unsigned char> present(T),desired(T);
    std::vector<DonorMask> masks(T);
    std::vector<std::vector<size_t>> owned(blocks.size());
    for(size_t x=0;x<T;++x) {
        const auto&occ=m.targets[x].occurrences;
        for(const auto&[b,mask]:occ)present[x]|=m.signature(mask,result.path[b]);
        const int a=segment_block[occ.front().first],b=segment_block[occ.back().first];
        if(a==b){owner[x]=a;owned[a].push_back(x);for(const auto&entry:occ)masks[x]|=entry.second;}
    }
    std::map<std::pair<int,int>,Score> pooled;
    for(const auto&e:m.contacts) {
        int a=owner[e.x],b=owner[e.y];if(a<0||b<0||a==b)continue;
        const int sign=phase_orientation(present[e.x])*phase_orientation(present[e.y])*e.relation;
        if(!sign)continue;
        if(a>b)std::swap(a,b);
        pooled[{a,b}]+=e.units*sign;
    }
    std::vector<Edge> edges;for(const auto&[key,w]:pooled)if(w)edges.push_back({key.first,key.second,w});
    report.edges=edges.size();if(edges.empty())return report;
    auto orientation=orient(edges,int(blocks.size()));
    desired=present;
    for(size_t x=0;x<T;++x)if(owner[x]>=0&&orientation[owner[x]]<0)
        desired[x]=((present[x]&1)<<1)|((present[x]&2)>>1);
    std::vector<int> proposed=result.path;
    auto swapped=[&](int z){return (z%m.H)*m.H+z/m.H;};
    for(size_t b=0;b<blocks.size();++b)for(int j=blocks[b].begin;j<blocks[b].end;++j)
        if(orientation[b]<0)proposed[j]=swapped(proposed[j]);
    for(size_t b=1;b<proposed.size();++b)report.proposed_recombinations+=m.charge(proposed[b-1],proposed[b]);
    // Evaluate the temporary proposal honestly, without changing m.budget.
    // Repeated strings may span blocks: recompute their actual presence.
    std::vector<unsigned char> pp(T);
    for(size_t x=0;x<T;++x){for(const auto&[b,mask]:m.targets[x].occurrences)pp[x]|=m.signature(mask,proposed[b]);
        report.unconstrained+=phase_reward_scale*m.weights.genotype*PresenceModel::genotype_reward(m.targets[x],pp[x]);}
    for(const auto&e:m.contacts)report.unconstrained+=e.units*m.weights.contact*phase_factor(pp[e.x],pp[e.y],e.relation,e.panel,m.neutral_phase,m.presence_phase);
    if(report.proposed_recombinations<=m.budget&&report.unconstrained>result.score) {
        if(m.score(proposed)!=report.unconstrained)throw std::logic_error("global proposal replay");
        result.path=proposed;result.score=report.unconstrained;
    }
    // Both coarse and fine repairs retain every ordered donor pair from
    // the fixed panel. Ranking each haplotype separately and keeping its
    // top two donors can discard the best feasible joint phase repair.
    // The dense resource DP evaluates these pairs in O(B R H^2), without
    // enumerating H^4 transitions or changing the recombination budget.
    RowExecutor pool(threads);
    pool.for_each(blocks.size(),[&](size_t b) {
        auto&c=blocks[b].choices;
        c.resize(size_t(m.D()));std::iota(c.begin(),c.end(),0);
    });
    // Twice the additive proposal score avoids fractional contact weights.
    // Cross-block contacts freeze the other endpoint to its proposed phase.
    // This is a search surrogate only: acceptance always uses m.score().
    std::vector<Score> field(T);
    std::vector<std::vector<size_t>> internal(blocks.size());
    for(size_t i=0;i<m.contacts.size();++i) {
        const auto&e=m.contacts[i];const int a=owner[e.x],b=owner[e.y];
        if(a>=0&&a==b){internal[a].push_back(i);continue;}
        const Score w=e.units*m.weights.contact*e.relation;
        if(a>=0)field[e.x]+=(b<0?2:1)*w*phase_orientation(desired[e.y]);
        if(b>=0)field[e.y]+=(a<0?2:1)*w*phase_orientation(desired[e.x]);
    }
    std::vector<std::vector<Score>> local(blocks.size());
    pool.for_each(blocks.size(),[&](size_t b) {
        local[b].resize(blocks[b].choices.size());
        for(size_t i=0;i<blocks[b].choices.size();++i) {
            const int z=blocks[b].choices[i];Score value=0;
            for(size_t x:owned[b]) {
                const auto p=m.signature(masks[x],z);
                value+=2*phase_reward_scale*m.weights.genotype*PresenceModel::genotype_reward(m.targets[x],p)+field[x]*phase_orientation(p);
            }
            for(size_t c:internal[b]) {const auto&e=m.contacts[c];value+=2*e.units*m.weights.contact*phase_factor(m.signature(masks[e.x],z),m.signature(masks[e.y],z),e.relation,e.panel,m.neutral_phase,m.presence_phase);}
            local[b][i]=value;
        }
    });
    auto fixed=repair(m,blocks,local,threads);report.trace_bytes=fixed.trace_bytes;
    if(fixed.choices.empty())throw std::logic_error("Incumbent missing from global repair domain");
    std::vector<int> candidate(m.blocks);
    for(size_t b=0;b<blocks.size();++b)for(int j=blocks[b].begin;j<blocks[b].end;++j)candidate[j]=blocks[b].choices[fixed.choices[b]];
    for(size_t b=1;b<candidate.size();++b)report.repaired_recombinations+=m.charge(candidate[b-1],candidate[b]);
    const Score candidate_score=m.score(candidate);
    report.repaired=candidate_score;
    if(candidate_score==impossible)throw std::logic_error("Global repair violates budget");
    if(candidate_score>result.score){result.path=std::move(candidate);result.score=candidate_score;}
    report.after=result.score;report.accepted=result.score>report.before;
    if(m.score(result.path)!=result.score||result.score>result.upper)throw std::logic_error("Global move final certificate replay");
    return report;
}
}
