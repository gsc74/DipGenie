#pragma once
#include "InferenceTypes.hpp"
#include "RowExecutor.hpp"
#include "PackedCompletion.hpp"
#include "DonorMask.hpp"
#include "PhaseFactor.hpp"
#include "JointContactCompletion.hpp"
#include "FrontierCopy.hpp"
#include "FlatHash.hpp"
#include <unordered_map>
#include <functional>
#include <set>
#include <tuple>
#include <memory>
#include <array>

namespace hic {
// Literal target presence, not occurrence count. A repeated target can occur
// in several segments and is scored ONCE, when its last occurrence is fixed.
struct PresenceTarget {
    bool hom=false;
    std::vector<std::pair<int,DonorMask>> occurrences; // segment, donor bitset
    Score reward_units=1; // fixed nonnegative dosage confidence; 1 recovers v1
    Score other_reward_units=0; // reward for the other nonzero dosage, not absence
    Score absent_reward_units=0;
};
struct PresenceContact {size_t x,y;Score units;int relation=1;bool panel=false;};
struct PresenceModel {
    int H=0,blocks=0,budget=0;
    std::vector<PresenceTarget> targets;
    std::vector<PresenceContact> contacts;
    NormalizedWeights weights{1,0};
    // false: v1 ordered charge (each haplotype's donor changes; the NP-hard
    // special case). true: diploid charge of the UNORDERED donor pair, i.e. the
    // fewer switches over both haplotype labelings of the next pair. Exchanging
    // the two haplotypes' donors is then free; phase across segments is decided
    // only by the linked read/Hi-C/panel phase factors, not by the budget.
    bool free_swaps=false;
    bool neutral_phase=false; // see phase_factor(); false = v1 factor 1+r*x*y
    bool presence_phase=false; // linked-presence factor; overrides neutral_phase
    // All-pair budget repair under ordered charge (a feasible restriction of
    // min-swapped charge): it never inserts score-neutral donor-pair swaps.
    bool ordered_repair=false;
    int D()const{return checked_pair_domain(size_t(H));}
    int ordered_charge(int a,int b)const{return(a/H!=b/H)+(a%H!=b%H);}
    int charge(int a,int b)const{
        const int c=ordered_charge(a,b);
        if(!free_swaps||!c)return c;
        return std::min(c,ordered_charge(a,(b%H)*H+b/H));
    }
    unsigned char signature(const DonorMask&mask,int z)const{return mask.test(size_t(z/H))+2*mask.test(size_t(z%H));}
    static Score agreement(unsigned char a,unsigned char b,int relation=1){
        return (a==1||a==2)&&(b==1||b==2)&&((a==b)==(relation==1));
    }
    static Score coverage(unsigned kind,unsigned char a){return kind==2?a==0:kind==1?a==3:(a==1||a==2);}
    static Score genotype_reward(const PresenceTarget&t,unsigned char a){
        return a?(coverage(t.hom,a)?t.reward_units:t.other_reward_units):t.absent_reward_units;
    }
    static Score target_upper(const PresenceTarget&t){return std::max({t.reward_units,t.other_reward_units,t.absent_reward_units});}
    // Full-graph bound U_0 = U_G + U_C: every coverage and phase factor at its
    // own maximum. Both parts are valid on any feasible set, including
    // exchange-feasible pairs, since no factor exceeds its maximum.
    Score genotype_bound()const {
        __int128 v=0;
        for(const auto&t:targets)v=score_sum(v,score_product(score_product(phase_reward_scale,target_upper(t)),weights.genotype));
        if(v>std::numeric_limits<Score>::max()/4)throw std::overflow_error("presence bound");
        return Score(v);
    }
    Score phase_bound()const {
        __int128 v=0;
        for(auto&e:contacts)v=score_sum(v,score_product(score_product(phase_reward_scale,e.units),weights.contact));
        if(v>std::numeric_limits<Score>::max()/4)throw std::overflow_error("presence bound");
        return Score(v);
    }
    Score universal_bound()const {
        const __int128 v=score_sum(genotype_bound(),phase_bound());
        if(v>std::numeric_limits<Score>::max()/4)throw std::overflow_error("presence bound");
        return Score(v);
    }
    void check()const {
        if(H<1||blocks<1||budget<0||weights.genotype<0||weights.contact<0)
            throw std::invalid_argument("presence dimensions/weights");
        (void)D();
        for(auto&t:targets) {
            if(t.reward_units<0||t.other_reward_units<0||t.absent_reward_units<0)throw std::invalid_argument("negative target reward");
            if(t.occurrences.empty())throw std::invalid_argument("target has no graph occurrence");
            int last=-1;
            for(const auto&[b,mask]:t.occurrences) {
                if(b<=last||b>=blocks||!mask||!mask.within(size_t(H)))throw std::invalid_argument("target occurrence");
                last=b;
            }
        }
        for(auto&e:contacts)if(e.x>=targets.size()||e.y>=targets.size()||e.x==e.y||e.units<0||(e.relation!=1&&e.relation!=-1))
            throw std::invalid_argument("contact endpoints/weight");
        (void)universal_bound();
    }
    Score score(const std::vector<int>&path)const {
        if(path.size()!=size_t(blocks))throw std::invalid_argument("incomplete path");
        int r=0;for(size_t b=0;b<path.size();++b){if(path[b]<0||path[b]>=D())throw std::invalid_argument("donor pair");if(b)r+=charge(path[b-1],path[b]);}
        if(r>budget)return impossible;
        std::vector<unsigned char> presence(targets.size());Score result=0;
        for(size_t x=0;x<targets.size();++x) {
            for(const auto&[b,mask]:targets[x].occurrences)presence[x]|=signature(mask,path[b]);
            result+=phase_reward_scale*genotype_reward(targets[x],presence[x])*weights.genotype;
        }
        for(auto&e:contacts)result+=phase_factor(presence[e.x],presence[e.y],e.relation,e.panel,neutral_phase,presence_phase)*e.units*weights.contact;
        return result;
    }
};

// Optional read-only search diagnostic. A witness is never inserted into the
// frontier, scored as an initializer, or used to change candidate ordering.
struct PresenceSearchAudit {
    struct Record {
        int completed=0,recombinations=0;
        Score witness_prefix=0,best_equivalent=impossible,priority=impossible,cutoff=impossible;
        size_t rank=0,candidates=0,retained=0,distinct_retained=0;
        bool generated=false,survived=false,equivalent_retained=false;
    };
    std::vector<int> witness;
    std::function<void(const Record&)> emit;
};

inline Result solve_presence(const PresenceModel&m,size_t beam,int threads=1,
                             std::function<void(int,const Result&)> progress={},
                             bool completion=false,const PresenceSearchAudit*audit=nullptr) {
    m.check();if(threads<1)throw std::invalid_argument("threads");
    if(audit&&(audit->witness.size()!=size_t(m.blocks)||m.score(audit->witness)==impossible||!audit->emit))
        throw std::invalid_argument("Search diagnostic requires a feasible complete witness and callback");
    // A segment with no target occurrence cannot change any evidence factor.
    // Delay switches across such segments to the next informative boundary.
    // This also avoids spending the beam on arbitrary donor labels before
    // the first observation. All donor transitions in this domain are valid.
    std::vector<int> occupied(m.blocks),order;
    for(auto&t:m.targets)for(const auto&[b,mask]:t.occurrences){(void)mask;occupied[b]=1;}
    for(int b=0;b<m.blocks;++b)if(occupied[b])order.push_back(b);
    if(order.empty()){Result r;r.path.assign(m.blocks,0);r.score=r.upper=0;if(progress)progress(m.blocks,r);return r;}
    if(order.size()<size_t(m.blocks)) {
        PresenceModel compact=m;compact.blocks=int(order.size());std::vector<int> remap(m.blocks,-1);
        for(size_t i=0;i<order.size();++i)remap[order[i]]=int(i);
        for(auto&t:compact.targets)for(auto&entry:t.occurrences)entry.first=remap[entry.first];
        PresenceSearchAudit compressed_audit;
        if(audit) {
            for(int b:order)compressed_audit.witness.push_back(audit->witness[b]);
            if(compact.score(compressed_audit.witness)!=m.score(audit->witness))
                throw std::logic_error("Witness empty-segment compression changed score");
            compressed_audit.emit=[&](const PresenceSearchAudit::Record&value){
                auto expanded=value;expanded.completed=order[value.completed-1]+1;audit->emit(expanded);
            };
        }
        auto r=solve_presence(compact,beam,threads,[&](int b,const Result&state){
            // The final informative segment also completes the empty suffix.
            // Keep progress in the original segment domain, including 100%.
            if(progress)progress(b==compact.blocks?m.blocks:order[b-1]+1,state);
        },completion,audit?&compressed_audit:nullptr);
        if(r.path.empty())return r;
        std::vector<int> expanded(m.blocks);size_t i=0;
        for(int b=0;b<m.blocks;++b){while(i+1<order.size()&&b>=order[i+1])++i;expanded[b]=r.path[i];}
        r.path=std::move(expanded);
        if(m.score(r.path)!=r.score)throw std::logic_error("Empty-segment compression changed score");
        return r;
    }
    const size_t T=m.targets.size();const int D=m.D();
    // Equality of the complete occurrence pattern implies equal presence on
    // every candidate in this segment domain. Alias only stored state: keep
    // each original coverage/contact factor and the original full-graph bound.
    std::vector<size_t> alias(T);
    std::vector<std::vector<std::pair<int,DonorMask>>> patterns;
    {
        // Store each distinct pattern once, not both in a tree key and in
        // patterns. Fingerprints only find a collision chain; exact pattern
        // equality retains the same aliases and first-occurrence ID order.
        FlatMap<uint64_t,size_t> heads;
        std::vector<size_t> next;
        for(size_t x=0;x<T;++x) {
            const auto&pattern=m.targets[x].occurrences;uint64_t hash=pattern.size();
            for(const auto&[block,mask]:pattern)
                hash=(hash^mask.fingerprint()^uint64_t(block))*0x9e3779b97f4a7c15ULL;
            auto [head,inserted]=heads.try_emplace(hash,SIZE_MAX);(void)inserted;
            size_t id=head->second;
            while(id!=SIZE_MAX&&patterns[id]!=pattern)id=next[id];
            if(id==SIZE_MAX) {
                id=patterns.size();next.push_back(head->second);head->second=id;
                patterns.push_back(pattern);
            }
            alias[x]=id;
        }
    }
    const size_t N=patterns.size();
    std::vector<int> first(N),last(N),use(N),old_position(N,-1);
    std::vector<std::vector<size_t>> born(m.blocks),close(m.blocks),contacts(m.blocks);
    std::vector<std::vector<std::pair<size_t,DonorMask>>> emit(m.blocks);
    for(size_t x=0;x<N;++x) {
        const auto& p=patterns[x];first[x]=p.front().first;last[x]=use[x]=p.back().first;
        born[first[x]].push_back(x);
        for(const auto&[b,mask]:p)emit[b].push_back({x,mask});
    }
    for(size_t x=0;x<T;++x)close[last[alias[x]]].push_back(x);
    // Fold equal factors by addition, not averaging: every original molecule
    // keeps its full weight. A within-class contact rewards heterozygosity.
    std::vector<PresenceContact> folded;
    {
        FlatMap<std::tuple<size_t,size_t,int,bool>,Score> pooled;
        for(const auto&e:m.contacts) {
            size_t x=alias[e.x],y=alias[e.y];if(x>y)std::swap(x,y);
            pooled[{x,y,e.relation,e.panel}]+=e.units;
        }
        for(const auto&[key,units]:pooled) {
            auto[x,y,relation,panel]=key;folded.push_back({x,y,units,relation,panel});
        }
        std::sort(folded.begin(),folded.end(),[](const auto&a,const auto&b){
            return std::tie(a.x,a.y,a.relation,a.panel)<std::tie(b.x,b.y,b.relation,b.panel);});
    }
    for(size_t i=0;i<folded.size();++i) {
        const auto&e=folded[i];const size_t x=e.x,y=e.y;
        int b=std::max(last[x],last[y]);contacts[b].push_back(i);
        use[x]=std::max(use[x],b);use[y]=std::max(use[y],b);
    }
    const auto phase_events=joint_contact_events(folded,first,last,m.blocks,m.weights.contact,completion&&beam&&!m.neutral_phase&&!m.presence_phase);
    std::vector<Score> suffix(m.blocks+1);
    for(int b=m.blocks-1;b>=0;--b) {
        suffix[b]=suffix[b+1];
        for(size_t x:close[b])suffix[b]+=phase_reward_scale*PresenceModel::target_upper(m.targets[x])*m.weights.genotype;
        for(size_t i:contacts[b])suffix[b]+=phase_reward_scale*folded[i].units*m.weights.contact;
    }
    PackedCompletion lookahead;
    RowExecutor pool(threads);
    if(completion&&beam) {
        PackedScoreRows local{size_t(D)};
        std::vector<Score> other(m.blocks+1);
        constexpr int tile=64;
        const int batch=std::min(threads,128)*tile;
        for(int first_block=0;first_block<m.blocks;first_block+=std::min(batch,m.blocks-first_block)) {
          const int count=std::min(batch,m.blocks-first_block);
          std::vector<PackedScoreRows> partial(size_t((count+tile-1)/tile),PackedScoreRows{size_t(D)});
          pool.for_each(partial.size(),[&](size_t task) {
           std::vector<Score> row(size_t(D),0);
           const int start=first_block+int(task)*tile,end=std::min(first_block+count,start+tile);
           for(int b=start;b<end;++b) {
            std::fill(row.begin(),row.end(),0);
            Score local_upper=0;
            std::map<std::pair<unsigned,DonorMask>,Score> groups;
            for(size_t x:close[b])if(first[alias[x]]==b) {
                const auto&t=m.targets[x];const auto&mask=patterns[alias[x]][0].second;
                groups[{t.hom,mask}]+=phase_reward_scale*t.reward_units*m.weights.genotype;
                groups[{!t.hom,mask}]+=phase_reward_scale*t.other_reward_units*m.weights.genotype;
                groups[{2,mask}]+=phase_reward_scale*t.absent_reward_units*m.weights.genotype;
                local_upper+=phase_reward_scale*PresenceModel::target_upper(t)*m.weights.genotype;
            }
            for(auto&[key,weight]:groups)for(int z=0;z<D;++z)
                row[z]+=weight*PresenceModel::coverage(key.first,m.signature(key.second,z));
            for(size_t i:contacts[b]) {
                const auto&e=folded[i];
                if(first[e.x]!=b||first[e.y]!=b)continue;
                const Score weight=e.units*m.weights.contact;local_upper+=phase_reward_scale*weight;
                for(int z=0;z<D;++z)row[z]+=weight*phase_factor(
                    m.signature(patterns[e.x][0].second,z),m.signature(patterns[e.y][0].second,z),e.relation,e.panel,m.neutral_phase,m.presence_phase);
            }
            other[b]=suffix[b]-suffix[b+1]-local_upper;
            if(other[b]<0)throw std::logic_error("Negative nonlocal completion bound");
            partial[task].append(row);
           }
          });
          for(auto&part:partial)local.append_rows(std::move(part));
        }
        for(int b=m.blocks-1;b>=0;--b)other[b]+=other[b+1];
        lookahead=PackedCompletion(m.H,m.blocks,m.budget,std::move(local),std::move(other),0,threads,m.free_swaps);
    }
    struct Trace {size_t parent;int choice;};
    struct Live {int r,last;Score score;size_t trace;PackedPresence values;};
    struct Row {
        int r=0;std::shared_ptr<const PackedPresence> before;
        uint64_t frontier_hash=0;
        size_t first_source=SIZE_MAX,last_source=SIZE_MAX,next_key=SIZE_MAX;
        Score*score=nullptr;size_t*parent=nullptr;
        const Score*completion=nullptr;
        size_t generated=0;
    };
    struct Cell {Score score;size_t parent,row;int z;Score priority=0;};
    auto better=[](const Cell&a,const Cell&b){
        if(a.priority!=b.priority)return a.priority>b.priority;
        return a.score!=b.score?a.score>b.score:std::make_tuple(a.parent,a.z,a.row)<std::make_tuple(b.parent,b.z,b.row);
    };
    std::vector<Trace> traces{{0,-1}};
    std::vector<Live> live{{0,-1,0,0,{}}};std::vector<size_t> previous;
    // Reuse the large destination matrices across segments. Allocation does
    // not initialize fundamental values; independent row workers initialize
    // and first-touch their own entries before evaluating predecessors.
    std::unique_ptr<Score[]> row_scores;
    std::unique_ptr<size_t[]> row_parents;
    size_t row_capacity=0;
    std::vector<Score> completion_rows;
    // Small resource domains admit direct lookup. Cap scratch independently
    // of the requested budget; chromosome-size or larger budgets fall back
    // to sorted sparse resources below.
    std::vector<size_t> resource_slots(m.budget<=65536?size_t(m.budget)+1:0,SIZE_MAX);
    // A reusable, at-most-half-full table locates (budget, fingerprint)
    // buckets. Full packed-state equality still decides every merge. Rows
    // are appended in the original source/delta order, never table order.
    struct RowKey {uint64_t epoch=0,hash=0;int resource=0;size_t first=SIZE_MAX;};
    struct RowSource {size_t state,next;};
    std::vector<RowKey> row_keys;
    std::vector<RowSource> row_sources;
    std::vector<DonorMask> current(N);Result result;Score omitted=impossible;
    std::vector<unsigned char> witness_presence(audit?N:0);
    Score witness_prefix=0;int witness_r=0;
    result.completion_penalty=lookahead.penalty;
    result.completion_relaxed_switches=lookahead.relaxed_switches;
    for(int b=0;b<m.blocks;++b) {
        // One serial reconstruction per tile; worker queries remain read-only.
        if(completion&&beam)lookahead.prepare(b);
        for(size_t i=0;i<previous.size();++i)old_position[previous[i]]=int(i);
        if(b)for(const auto&[x,mask]:emit[b-1]){(void)mask;current[x]=0;}
        for(const auto&[x,mask]:emit[b])current[x]=mask;
        // Sorted, disjoint lists: maintain the contact frontier in linear
        // time instead of sorting all long-range assignments every segment.
        std::vector<size_t> active;active.reserve(previous.size()+born[b].size());
        auto birth=born[b].begin();
        for(size_t x:previous)if(use[x]>b) {
            while(birth!=born[b].end()&&*birth<x) {
                if(use[*birth]>b)active.push_back(*birth);
                ++birth;
            }
            active.push_back(x);
        }
        for(;birth!=born[b].end();++birth)if(use[*birth]>b)active.push_back(*birth);
        result.max_frontier=std::max(result.max_frontier,active.size());
        std::vector<int> previous_positions;previous_positions.reserve(active.size());
        std::vector<std::pair<size_t,DonorMask>> active_emissions;
        for(size_t i=0;i<active.size();++i) {
            const size_t x=active[i];previous_positions.push_back(old_position[x]);
            if(current[x])active_emissions.push_back({i,current[x]});
        }
        const FrontierCopy frontier_copy(previous_positions);
        // Targets confined to this segment have no dependence on history.
        // Aggregate equal donor masks before computing the joint unary row.
        std::map<std::pair<unsigned,DonorMask>,Score> local_groups;
        std::map<std::pair<size_t,unsigned>,Score> dynamic_groups;
        for(size_t x:close[b]) {
            const auto&t=m.targets[x];
            if(first[alias[x]]==b){
                local_groups[{t.hom,current[alias[x]]}]+=t.reward_units;
                local_groups[{!t.hom,current[alias[x]]}]+=t.other_reward_units;
                local_groups[{2,current[alias[x]]}]+=t.absent_reward_units;
            }else{
                dynamic_groups[{alias[x],t.hom}]+=t.reward_units;
                dynamic_groups[{alias[x],!t.hom}]+=t.other_reward_units;
                dynamic_groups[{alias[x],2}]+=t.absent_reward_units;
            }
        }
        std::vector<Score> unary(D);
        for(auto&[key,count]:local_groups)for(int z=0;z<D;++z)
            unary[z]+=phase_reward_scale*count*m.weights.genotype*PresenceModel::coverage(key.first,m.signature(key.second,z));
        // A factor sees only four presence states per endpoint. Compile its
        // exact integer reward once per segment, not once per predecessor
        // and donor pair. Predecessors select a small lookup row; successors
        // then need only a byte lookup and addition. This preserves every
        // factor, including neutral absent/homozygous phase observations.
        struct TransitionFactor {
            int first=-1,second=-1;
            size_t stride=4;
            std::vector<unsigned char> emitted;
            std::vector<Score> rewards;
        };
        std::vector<TransitionFactor> factors;
        std::vector<JointContactEvent> past_phase_events;
        for(auto it=dynamic_groups.begin();it!=dynamic_groups.end();) {
            const size_t x=it->first.first;
            std::array<Score,3> weights{};
            do {
                weights[it->first.second]+=phase_reward_scale*it->second*m.weights.genotype;
                ++it;
            } while(it!=dynamic_groups.end()&&it->first.first==x);
            TransitionFactor f;f.first=old_position[x];f.emitted.resize(D);f.rewards.resize(16);
            for(int z=0;z<D;++z)f.emitted[z]=m.signature(current[x],z);
            for(unsigned old=0;old<4;++old)for(unsigned now=0;now<4;++now) {
                const unsigned presence=old|now;
                f.rewards[4*old+now]=weights[0]*PresenceModel::coverage(false,presence)+
                                     weights[1]*PresenceModel::coverage(true,presence)+
                                     weights[2]*PresenceModel::coverage(2,presence);
            }
            factors.push_back(std::move(f));
        }
        for(size_t i:contacts[b]) {
            const auto&e=folded[i];const Score weight=e.units*m.weights.contact;
            // A contact whose two endpoints are born here is independent
            // of the predecessor history and belongs in the shared unary row.
            if(first[e.x]==b&&first[e.y]==b) {
                for(int z=0;z<D;++z)unary[z]+=weight*phase_factor(
                    m.signature(current[e.x],z),m.signature(current[e.y],z),e.relation,e.panel,m.neutral_phase,m.presence_phase);
                continue;
            }
            TransitionFactor f;f.first=old_position[e.x];f.second=old_position[e.y];
            f.stride=16;f.emitted.resize(D);f.rewards.resize(256);
            for(int z=0;z<D;++z)f.emitted[z]=m.signature(current[e.x],z)|(m.signature(current[e.y],z)<<2);
            for(unsigned old=0;old<16;++old)for(unsigned now=0;now<16;++now) {
                const unsigned presence=old|now;
                f.rewards[16*old+now]=weight*phase_factor(presence&3,presence>>2,e.relation,e.panel,m.neutral_phase,m.presence_phase);
            }
            factors.push_back(std::move(f));
        }
        for(const auto&e:phase_events[b]) {
            if(first[e.x]<b&&first[e.y]<b) {
                if(old_position[e.x]<0||old_position[e.y]<0)throw std::logic_error("Missing pending phase endpoint");
                past_phase_events.push_back(e);
                continue;
            }
            if(first[e.x]==b&&first[e.y]==b) {
                for(int z=0;z<D;++z)unary[z]+=e.value(m.signature(current[e.x],z),m.signature(current[e.y],z));
                continue;
            }
            TransitionFactor f;f.first=old_position[e.x];f.second=old_position[e.y];
            f.stride=16;f.emitted.resize(D);f.rewards.resize(256);
            for(int z=0;z<D;++z)f.emitted[z]=m.signature(current[e.x],z)|(m.signature(current[e.y],z)<<2);
            for(unsigned old=0;old<16;++old)for(unsigned now=0;now<16;++now) {
                const unsigned presence=old|now;
                f.rewards[16*old+now]=e.value(presence&3,presence>>2);
            }
            factors.push_back(std::move(f));
        }
        // Position-sensitive XOR fingerprint permits sparse OR-emission
        // updates. Hashes never establish equality: full packed keys do.
        auto word_hash=[](uint64_t word,size_t position) {
            uint64_t x=word+0x9e3779b97f4a7c15ULL*(position+1);
            x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
            x=(x^(x>>27))*0x94d049bb133111ebULL;
            return x^(x>>31);
        };
        std::vector<std::shared_ptr<const PackedPresence>> gathered(live.size());
        std::vector<uint64_t> hashes(live.size());
        pool.for_each(live.size(),[&](size_t i) {
            auto values=std::make_shared<PackedPresence>(frontier_copy.gather(live[i].values));
            uint64_t hash=values->size();
            for(size_t w=0;w<values->words().size();++w)
                hash^=word_hash(values->words()[w],w);
            gathered[i]=std::move(values);hashes[i]=hash;
        });
        // Hashing only accelerates lookup: full equality resolves collisions.
        // Budget variants share an immutable frontier instead of storing it
        // again in the map key and in every destination row.
        if(live.size()>SIZE_MAX/6)throw std::overflow_error("Row grouping capacity");
        size_t key_capacity=8;
        while(key_capacity<6*live.size()) {
            if(key_capacity>SIZE_MAX/2)throw std::overflow_error("Row key capacity");
            key_capacity*=2;
        }
        if(row_keys.size()<key_capacity)row_keys.resize(key_capacity);
        const uint64_t key_epoch=uint64_t(b)+1;
        row_sources.clear();row_sources.reserve(3*live.size());
        std::vector<Row> rows;
        rows.reserve(3*live.size());
        for(size_t i=0;i<live.size();++i) {
            const auto&s=live[i];const auto&before=gathered[i];
            for(int delta=0;delta<=(b?2:0)&&s.r+delta<=m.budget;++delta) {
                const int resource=s.r+delta;
                size_t slot=size_t(word_hash(hashes[i],size_t(resource)))&(row_keys.size()-1);
                while(row_keys[slot].epoch==key_epoch&&
                      (row_keys[slot].resource!=resource||row_keys[slot].hash!=hashes[i]))
                    slot=(slot+1)&(row_keys.size()-1);
                auto&key=row_keys[slot];
                if(key.epoch!=key_epoch)key={key_epoch,hashes[i],resource,SIZE_MAX};
                size_t ri=SIZE_MAX;
                for(size_t candidate=key.first;candidate!=SIZE_MAX;candidate=rows[candidate].next_key)
                    if(*rows[candidate].before==*before){ri=candidate;break;}
                if(ri==SIZE_MAX) {
                    ri=rows.size();Row row;row.r=resource;row.before=before;row.frontier_hash=hashes[i];
                    row.next_key=key.first;key.first=ri;
                    rows.push_back(std::move(row));
                }
                auto&row=rows[ri];const size_t link=row_sources.size();
                row_sources.push_back({i,SIZE_MAX});
                if(row.first_source==SIZE_MAX)row.first_source=link;
                else row_sources[row.last_source].next=link;
                row.last_source=link;
            }
        }
        if(rows.size()>row_capacity) {
            const size_t extra=std::max(size_t(1),row_capacity/2);
            if(row_capacity>SIZE_MAX-extra)throw std::overflow_error("Destination row capacity");
            row_capacity=std::max(rows.size(),row_capacity+extra);
            if(beam) {
                if(beam<=SIZE_MAX/3)row_capacity=std::min(row_capacity,3*beam);
            }
            if(row_capacity>SIZE_MAX/size_t(D)/std::max(sizeof(Score),sizeof(size_t)))
                throw std::overflow_error("Destination score storage");
            const size_t entries=row_capacity*size_t(D);
            row_scores.reset();row_parents.reset();
            row_scores.reset(new Score[entries]);row_parents.reset(new size_t[entries]);
        }
        pool.for_each(rows.size(),[&](size_t ri) {
            auto&row=rows[ri];
            row.score=row_scores.get()+ri*size_t(D);row.parent=row_parents.get()+ri*size_t(D);
            std::fill_n(row.score,D,impossible);std::fill_n(row.parent,D,size_t(0));
            struct Lookup {const unsigned char*emitted;const Score*rewards;};
            thread_local std::vector<Lookup> lookup;
            lookup.resize(factors.size());
            for(size_t link=row.first_source;link!=SIZE_MAX;link=row_sources[link].next) {
                const size_t si=row_sources[link].state;
                const auto&s=live[si];
                Score starting=s.score;
                for(const auto&e:past_phase_events)starting+=e.value(s.values[old_position[e.x]],s.values[old_position[e.y]]);
                for(size_t i=0;i<factors.size();++i) {
                    const auto&f=factors[i];
                    const unsigned old=(f.first<0?0:s.values[f.first])|
                                       ((f.second<0?0:s.values[f.second])<<2);
                    lookup[i]={f.emitted.data(),f.rewards.data()+old*f.stride};
                }
                auto visit=[&](int z) {
                    Score value=starting+unary[z];
                    for(const auto&f:lookup)value+=f.rewards[f.emitted[z]];
                    ++row.generated;
                    if(value>row.score[z]||(value==row.score[z]&&s.trace<row.parent[z])){row.score[z]=value;row.parent[z]=s.trace;}
                };
                if(!b) {
                    // Fix only the global exchange symmetry, never relabel
                    // the haplotypes at a later segment.
                    for(int a=0;a<m.H;++a)for(int c=a;c<m.H;++c)visit(a*m.H+c);
                } else {
                    // Enumerate exactly the successors in this resource
                    // layer, instead of scanning H^2 cells three times and
                    // rejecting the wrong switch counts in each scan.
                    const int a=s.last/m.H,c=s.last%m.H,delta=row.r-s.r;
                    // Each successor is visited exactly once, in the layer of
                    // its (ordered or unordered) charge. Non-swapped successors
                    // are visited first, so exact ties keep the orientation.
                    const bool swaps=m.free_swaps&&a!=c;
                    if(delta==0){visit(s.last);if(swaps)visit(c*m.H+a);}
                    else if(delta==1) {
                        for(int h=0;h<m.H;++h)if(h!=c)visit(a*m.H+h);
                        for(int h=0;h<m.H;++h)if(h!=a)visit(h*m.H+c);
                        if(swaps) {
                            // (h,a) and (c,h); (a,a) and (c,c) are already above.
                            for(int h=0;h<m.H;++h)if(h!=c&&h!=a)visit(h*m.H+a);
                            for(int h=0;h<m.H;++h)if(h!=a&&h!=c)visit(c*m.H+h);
                        }
                    } else if(delta==2) {
                        for(int h=0;h<m.H;++h)if(h!=a&&!(m.free_swaps&&h==c))
                            for(int k=0;k<m.H;++k)if(k!=c&&!(m.free_swaps&&k==a))visit(h*m.H+k);
                    } else throw std::logic_error("Invalid joint transition charge");
                }
            }
        });
        if(completion&&beam) {
            std::vector<int> resources;resources.reserve(rows.size());
            if(!resource_slots.empty()) {
                for(const auto&row:rows)if(resource_slots[size_t(row.r)]==SIZE_MAX) {
                    resource_slots[size_t(row.r)]=resources.size();resources.push_back(row.r);
                }
            } else {
                for(const auto&row:rows)resources.push_back(row.r);
                std::sort(resources.begin(),resources.end());
                resources.erase(std::unique(resources.begin(),resources.end()),resources.end());
            }
            if(resources.size()>SIZE_MAX/size_t(D))throw std::overflow_error("Completion row cache");
            completion_rows.resize(resources.size()*size_t(D));
            pool.for_each(resources.size(),[&](size_t i) {
                lookahead.remaining_row(b,m.budget-resources[i],suffix[b+1],completion_rows.data()+i*size_t(D));
            });
            for(auto&row:rows) {
                const size_t i=resource_slots.empty()?size_t(std::lower_bound(resources.begin(),resources.end(),row.r)-resources.begin()):resource_slots[size_t(row.r)];
                row.completion=completion_rows.data()+i*size_t(D);
            }
            if(!resource_slots.empty())for(int r:resources)resource_slots[size_t(r)]=SIZE_MAX;
        }
        auto make_cell=[&](size_t i,int z) {
            Cell cell{rows[i].score[z],rows[i].parent[z],i,z};
            cell.priority=cell.score+(rows[i].completion?rows[i].completion[z]:suffix[b+1]);
            return cell;
        };
        Cell witness_cell{impossible,0,0,-1};PresenceSearchAudit::Record observation;
        if(audit) {
            const int z=audit->witness[b];
            if(b)witness_r+=m.charge(audit->witness[b-1],z);
            for(const auto&[x,mask]:emit[b])witness_presence[x]|=m.signature(mask,z);
            for(size_t x:close[b])witness_prefix+=phase_reward_scale*m.weights.genotype*
                PresenceModel::genotype_reward(m.targets[x],witness_presence[alias[x]]);
            for(size_t i:contacts[b]){const auto&e=folded[i];witness_prefix+=e.units*m.weights.contact*
                phase_factor(witness_presence[e.x],witness_presence[e.y],e.relation,e.panel,m.neutral_phase,m.presence_phase);}
            for(const auto&e:phase_events[b])witness_prefix+=e.value(witness_presence[e.x],witness_presence[e.y]);
            observation.completed=b+1;observation.recombinations=witness_r;observation.witness_prefix=witness_prefix;
            // Global exchange of the two haplotype labels is the same
            // physical solution and must not be reported as a search loss.
            for(int flip=0;flip<2;++flip)for(size_t i=0;i<rows.size();++i) {
                const int oriented=flip?(z%m.H)*m.H+z/m.H:z;
                const auto&row=rows[i];if(row.r!=witness_r||row.score[oriented]==impossible)continue;
                bool equal=true;
                for(size_t j=0;j<active.size();++j) {
                    const auto p=witness_presence[active[j]];
                    const auto expected=flip?((p&1)<<1)|((p&2)>>1):p;
                    if(((*row.before)[j]|m.signature(current[active[j]],oriented))!=expected){equal=false;break;}
                }
                if(equal){auto cell=make_cell(i,oriented);if(witness_cell.z<0||better(cell,witness_cell))witness_cell=cell;}
            }
            observation.generated=witness_cell.z>=0;
            if(observation.generated){observation.best_equivalent=witness_cell.score;observation.priority=witness_cell.priority;observation.rank=1;}
        }
        struct RankedRow {
            std::vector<Cell> best_cells;
            Score omitted=impossible;
            size_t candidates=0,above_witness=0;
        };
        std::vector<RankedRow> ranked;
        bool emission_can_merge=false;
        if(!active_emissions.empty()&&rows.size()>1) {
            // OR-emissions can erase differences only at their touched
            // presence slots. If (resource, untouched slots) differs for
            // every row, no destination donor pair can merge two rows.
            // A hash collision merely selects the slower exact merge below;
            // a hash is NEVER used to establish state equality.
            std::vector<std::pair<size_t,uint64_t>> erased;
            for(const auto&[i,mask]:active_emissions) {
                (void)mask;
                const uint64_t bits=uint64_t(3)<<(2*(i%32));
                if(!erased.empty()&&erased.back().first==i/32)erased.back().second|=bits;
                else erased.push_back({i/32,bits});
            }
            std::vector<std::pair<int,uint64_t>> keys;keys.reserve(rows.size());
            for(const auto&row:rows) {
                uint64_t hash=row.frontier_hash;
                for(const auto&[w,bits]:erased) {
                    const auto old=row.before->words()[w];
                    hash^=word_hash(old,w)^word_hash(old&~bits,w);
                }
                keys.emplace_back(row.r,hash);
            }
            std::sort(keys.begin(),keys.end());
            emission_can_merge=std::adjacent_find(keys.begin(),keys.end())!=keys.end();
        }
        if(!emission_can_merge) {
        // Distinct (r,before,z) cells remain distinct after emission. Retain
        // contiguous row scans rather than transposed H^2-by-row scans.
        ranked.resize(rows.size());
        // Row maxima are distinct generated states. Their K-th best is
        // a safe cutoff: no candidate worse than it can enter global top K.
        // This avoids constructing/sorting full Cell records for candidates
        // that are already excluded. It is exact, including all tie fields.
        Cell cutoff{impossible,0,0,-1};
        if(beam&&rows.size()>=beam&&rows.size()>1) {
            std::vector<Cell> maxima(rows.size(),Cell{impossible,0,0,-1});
            pool.for_each(rows.size(),[&](size_t i) {
                auto&best=maxima[i];
                for(int z=0;z<D;++z)if(rows[i].score[z]!=impossible) {
                    const auto cell=make_cell(i,z);
                    if(best.z<0||better(cell,best))best=cell;
                }
            });
            maxima.erase(std::remove_if(maxima.begin(),maxima.end(),[](const Cell&c){return c.z<0;}),maxima.end());
            if(maxima.size()>=beam) {
                std::nth_element(maxima.begin(),maxima.begin()+beam-1,maxima.end(),better);
                cutoff=maxima[beam-1];
            }
        }
        pool.for_each(rows.size(),[&](size_t i) {
            auto&out=ranked[i];auto&top=out.best_cells;
            // Partition once instead of maintaining a heap for each accepted
            // cell. Scratch is local to a worker; retain only K cells per row.
            thread_local std::vector<Cell> cells;
            cells.clear();cells.reserve(size_t(D));
            for(int z=0;z<D;++z)if(rows[i].score[z]!=impossible) {
                const Cell cell=make_cell(i,z);++out.candidates;
                if(audit&&observation.generated&&better(cell,witness_cell))++out.above_witness;
                if(cutoff.z>=0&&better(cutoff,cell)) {
                    out.omitted=std::max(out.omitted,cell.priority);continue;
                }
                cells.push_back(cell);
            }
            if(beam&&cells.size()>beam) {
                std::nth_element(cells.begin(),cells.begin()+beam,cells.end(),better);
                for(size_t j=beam;j<cells.size();++j)out.omitted=std::max(out.omitted,cells[j].priority);
                top.assign(cells.begin(),cells.begin()+beam);
            } else top.assign(cells.begin(),cells.end());
        });
        } else {
            // OR can collapse different incoming frontiers. Merge those
            // cells BEFORE any top-K truncation. A shard owns one destination
            // donor pair, so distinct shards cannot contain the same state.
            ranked.resize(size_t(D));
            Cell cutoff{impossible,0,0,-1};
            if(beam&&size_t(D)>=beam&&!audit) {
                std::vector<Cell> maxima(size_t(D),Cell{impossible,0,0,-1});
                pool.for_each(size_t(D),[&](size_t z) {
                    auto&best=maxima[z];
                    for(size_t i=0;i<rows.size();++i)if(rows[i].score[z]!=impossible) {
                        const auto c=make_cell(i,int(z));
                        if(best.z<0||better(c,best))best=c;
                    }
                });
                maxima.erase(std::remove_if(maxima.begin(),maxima.end(),[](const Cell&c){return c.z<0;}),maxima.end());
                if(maxima.size()>=beam) {
                    std::nth_element(maxima.begin(),maxima.begin()+beam-1,maxima.end(),better);
                    cutoff=maxima[beam-1];
                }
            }
            pool.for_each(size_t(D),[&](size_t z) {
                auto&out=ranked[z];
                // Emission masks occupy only words touched at this segment;
                // do not allocate a complete frontier for every candidate.
                thread_local std::vector<std::pair<size_t,uint64_t>> emissions;
                emissions.clear();
                for(const auto&[i,mask]:active_emissions) {
                    const uint64_t bits=uint64_t(m.signature(mask,int(z)))<<(2*(i%32));
                    if(!bits)continue;
                    if(!emissions.empty()&&emissions.back().first==i/32)emissions.back().second|=bits;
                    else emissions.push_back({i/32,bits});
                }
                struct Keyed {Cell cell;int resource;uint64_t hash;};
                thread_local std::vector<Keyed> cells;
                cells.clear();cells.reserve(rows.size());
                for(size_t i=0;i<rows.size();++i)if(rows[i].score[z]!=impossible) {
                    const auto c=make_cell(i,int(z));
                    if(cutoff.z>=0&&better(cutoff,c)) {
                        ++out.candidates; // safe omissions; may overcount merged low-priority cells
                        out.omitted=std::max(out.omitted,c.priority);continue;
                    }
                    uint64_t hash=rows[i].frontier_hash;
                    for(const auto&[w,bits]:emissions) {
                        const auto old=rows[i].before->words()[w];
                        hash^=word_hash(old,w)^word_hash(old|bits,w);
                    }
                    cells.push_back({c,rows[i].r,hash});
                }
                std::sort(cells.begin(),cells.end(),[](const Keyed&a,const Keyed&b) {
                    return std::tie(a.resource,a.hash,a.cell.row)<std::tie(b.resource,b.hash,b.cell.row);
                });
                auto equal_after=[&](size_t a,size_t b) {
                    const auto&x=rows[a].before->words();const auto&y=rows[b].before->words();
                    size_t e=0;
                    for(size_t w=0;w<x.size();++w) {
                        const uint64_t bits=e<emissions.size()&&emissions[e].first==w?emissions[e++].second:0;
                        if((x[w]|bits)!=(y[w]|bits))return false;
                    }
                    return true;
                };
                auto&top=out.best_cells;
                for(size_t begin=0;begin<cells.size();) {
                    size_t end=begin+1;
                    while(end<cells.size()&&cells[end].resource==cells[begin].resource&&cells[end].hash==cells[begin].hash)++end;
                    const size_t first=top.size();
                    for(size_t i=begin;i<end;++i) {
                        size_t j=first;
                        for(;j<top.size();++j)if(equal_after(cells[i].cell.row,top[j].row))break;
                        if(j==top.size())top.push_back(cells[i].cell);
                        else if(better(cells[i].cell,top[j]))top[j]=cells[i].cell;
                    }
                    begin=end;
                }
                out.candidates+=top.size();
                if(audit&&observation.generated)for(const auto&c:top)
                    if(better(c,witness_cell))++out.above_witness;
                if(beam&&top.size()>beam) {
                    std::nth_element(top.begin(),top.begin()+beam,top.end(),better);
                    for(size_t i=beam;i<top.size();++i)out.omitted=std::max(out.omitted,top[i].priority);
                    top.resize(beam);
                }
            });
        }
        // A global top-K state must be in its shard's top K. Keep a strict
        // total width, independent of R. Resource winners outside that set
        // are approximate omissions and contribute their completion bounds.
        std::vector<Cell> keep;
        size_t candidates=0;
        for(const auto&row:rows)result.generated+=row.generated;
        for(size_t i=0;i<ranked.size();++i) {
            const auto&rank=ranked[i];
            candidates+=rank.candidates;
            omitted=std::max(omitted,rank.omitted);
            if(audit){observation.candidates+=rank.candidates;observation.rank+=rank.above_witness;}
        }
        if(beam) {
            // Hierarchical exact top-K selection among generated candidates.
            // Every omitted candidate contributes a completion bound.
            constexpr size_t group_rows=32;
            struct Group {std::vector<Cell> cells;Score omitted=impossible;};
            std::vector<Group> groups((ranked.size()+group_rows-1)/group_rows);
            auto select=[&](std::vector<Cell>&cells,Score&lost) {
                if(cells.size()<=beam)return;
                std::nth_element(cells.begin(),cells.begin()+beam,cells.end(),better);
                for(size_t j=beam;j<cells.size();++j)
                    lost=std::max(lost,cells[j].priority);
                cells.resize(beam);
            };
            pool.for_each(groups.size(),[&](size_t g) {
                thread_local std::vector<Cell> cells;cells.clear();
                const size_t end=std::min(ranked.size(),(g+1)*group_rows);
                size_t size=0;for(size_t i=g*group_rows;i<end;++i)size+=ranked[i].best_cells.size();
                cells.reserve(size);
                for(size_t i=g*group_rows;i<end;++i)
                    cells.insert(cells.end(),ranked[i].best_cells.begin(),ranked[i].best_cells.end());
                select(cells,groups[g].omitted);
                groups[g].cells.assign(cells.begin(),cells.end());
            });
            for(auto&group:groups) {
                omitted=std::max(omitted,group.omitted);
                keep.insert(keep.end(),group.cells.begin(),group.cells.end());
            }
            select(keep,omitted);
        } else for(const auto&rank:ranked)keep.insert(keep.end(),rank.best_cells.begin(),rank.best_cells.end());
        result.discarded+=candidates-keep.size();
        std::sort(keep.begin(),keep.end(),better);live.clear();live.resize(keep.size());
        const size_t trace_begin=traces.size();
        if(keep.size()>SIZE_MAX-trace_begin)throw std::overflow_error("Trace storage");
        traces.resize(trace_begin+keep.size());
        pool.for_each(keep.size(),[&](size_t i) {
            const auto&cell=keep[i];const auto&row=rows[cell.row];
            Live s{row.r,cell.z,cell.score,trace_begin+i,*row.before};
            for(const auto&[i,mask]:active_emissions)s.values.add(i,m.signature(mask,cell.z));
            traces[trace_begin+i]={cell.parent,cell.z};live[i]=std::move(s);
        });
        if(audit) {
            std::set<std::tuple<int,int,std::vector<unsigned char>>> distinct;
            for(const auto&s:live) {
                distinct.emplace(s.r,s.last,s.values.unpack());
                if(s.r!=witness_r)continue;
                for(int flip=0;flip<2;++flip) {
                    const int z=audit->witness[b],oriented=flip?(z%m.H)*m.H+z/m.H:z;
                    if(s.last!=oriented)continue;
                    bool equal=true;
                    for(size_t j=0;j<active.size();++j) {
                        const auto p=witness_presence[active[j]];
                        const auto expected=flip?((p&1)<<1)|((p&2)>>1):p;
                        if(s.values[j]!=expected){equal=false;break;}
                    }
                    if(equal) {
                        observation.equivalent_retained=true;
                        if(s.score>=witness_prefix)observation.survived=true;
                    }
                }
            }
            observation.retained=live.size();observation.distinct_retained=distinct.size();
            if(!keep.empty())observation.cutoff=keep.back().priority;
            audit->emit(observation);
        }
        result.peak=std::max(result.peak,live.size());
        for(size_t x:previous)old_position[x]=-1;
        previous=std::move(active);
        if(live.empty())return result;
        if(progress&&(b%1000==0||b+1==m.blocks))progress(b+1,result);
    }
    for(size_t t=live.front().trace;t;t=traces[t].parent)result.path.push_back(traces[t].choice);
    std::reverse(result.path.begin(),result.path.end());result.score=m.score(result.path);
    if(result.score!=live.front().score)throw std::logic_error("frontier score differs from literal whole-pair score");
    result.upper=std::max(result.score,omitted);return result;
}

// Reverse only the order of segment decisions, not graph edges or DNA.
// OR-presence, undirected phase factors and the symmetric donor-switch
// charge are unchanged after reversing the returned decision sequence.
inline PresenceModel reverse_decisions(const PresenceModel&m) {
    m.check();PresenceModel reversed=m;
    for(auto&t:reversed.targets) {
        for(auto&entry:t.occurrences)entry.first=m.blocks-1-entry.first;
        std::reverse(t.occurrences.begin(),t.occurrences.end());
    }
    reversed.check();return reversed;
}

struct BidirectionalResult {
    Result best;
    Score forward_score=impossible,reverse_score=impossible;
    bool selected_reverse=false;
    int passes=0;
};
inline BidirectionalResult solve_presence_bidirectional(const PresenceModel&m,size_t beam,int threads=1,
        std::function<void(const char*,int,const Result&)> progress={}) {
    BidirectionalResult output;
    auto forward=solve_presence(m,beam,threads,[&](int b,const Result&r){if(progress)progress("forward",b,r);},true);
    output.passes=1;output.forward_score=forward.path.empty()?impossible:forward.score;
    // An exact segment-domain certificate makes another traversal redundant.
    if(!forward.path.empty()&&forward.upper==forward.score){output.best=std::move(forward);return output;}
    auto reversed=reverse_decisions(m);
    auto backward=solve_presence(reversed,beam,threads,[&](int b,const Result&r){if(progress)progress("reverse",b,r);},true);
    output.passes=2;output.reverse_score=backward.path.empty()?impossible:backward.score;
    std::reverse(backward.path.begin(),backward.path.end());
    if(!backward.path.empty()) {
        if(m.score(backward.path)!=backward.score)throw std::logic_error("Decision reversal changed objective");
        // Canonical global labels; never exchange haplotypes at local boundaries.
        auto swapped=backward.path;for(auto&z:swapped)z=(z%m.H)*m.H+z/m.H;
        if(swapped<backward.path)backward.path=std::move(swapped);
    }
    const bool forward_ok=!forward.path.empty(),backward_ok=!backward.path.empty();
    const auto upper=forward_ok&&backward_ok?std::min(forward.upper,backward.upper):
                     forward_ok?forward.upper:backward.upper;
    const auto generated=forward.generated+backward.generated,discarded=forward.discarded+backward.discarded;
    const auto peak=std::max(forward.peak,backward.peak),frontier=std::max(forward.max_frontier,backward.max_frontier);
    output.selected_reverse=backward_ok&&(!forward_ok||backward.score>forward.score);
    output.best=output.selected_reverse?std::move(backward):std::move(forward);
    output.best.upper=upper;output.best.generated=generated;output.best.discarded=discarded;
    output.best.peak=peak;output.best.max_frontier=frontier;
    if(!output.best.path.empty()&&output.best.score>output.best.upper)throw std::logic_error("Invalid bidirectional bound");
    return output;
}
}
