#include "PresenceDP.hpp"
#include <cassert>
#include <random>
#include <iostream>
using namespace hic;
inline Result serial_solve_presence(const PresenceModel&m,size_t beam,int threads=1,
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
    if(order.empty()){Result r;r.path.assign(m.blocks,0);r.score=r.upper=0;return r;}
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
        auto r=serial_solve_presence(compact,beam,threads,[&](int b,const Result&state){if(progress)progress(order[b-1]+1,state);},completion,audit?&compressed_audit:nullptr);
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
        std::map<std::vector<std::pair<int,DonorMask>>,size_t> lookup;
        for(size_t x=0;x<T;++x) {
            auto [it,inserted]=lookup.emplace(m.targets[x].occurrences,patterns.size());
            if(inserted)patterns.push_back(m.targets[x].occurrences);
            alias[x]=it->second;
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
        std::map<std::tuple<size_t,size_t,int,bool>,Score> pooled;
        for(const auto&e:m.contacts) {
            size_t x=alias[e.x],y=alias[e.y];if(x>y)std::swap(x,y);
            pooled[{x,y,e.relation,e.panel}]+=e.units;
        }
        for(const auto&[key,units]:pooled) {
            auto[x,y,relation,panel]=key;folded.push_back({x,y,units,relation,panel});
        }
    }
    for(size_t i=0;i<folded.size();++i) {
        const auto&e=folded[i];const size_t x=e.x,y=e.y;
        int b=std::max(last[x],last[y]);contacts[b].push_back(i);
        use[x]=std::max(use[x],b);use[y]=std::max(use[y],b);
    }
    const auto phase_events=joint_contact_events(folded,first,last,m.blocks,m.weights.contact,completion&&beam);
    std::vector<Score> suffix(m.blocks+1);
    for(int b=m.blocks-1;b>=0;--b) {
        suffix[b]=suffix[b+1];
        for(size_t x:close[b])suffix[b]+=phase_reward_scale*PresenceModel::target_upper(m.targets[x])*m.weights.genotype;
        for(size_t i:contacts[b])suffix[b]+=phase_reward_scale*folded[i].units*m.weights.contact;
    }
    PackedCompletion lookahead;
    if(completion&&beam) {
        PackedScoreRows local{size_t(D)};
        std::vector<Score> row(size_t(D),0),other(m.blocks+1);
        for(int b=0;b<m.blocks;++b) {
            std::fill(row.begin(),row.end(),0);
            Score local_upper=0;
            std::map<std::pair<bool,DonorMask>,Score> groups;
            for(size_t x:close[b])if(first[alias[x]]==b) {
                const auto&t=m.targets[x];const auto&mask=patterns[alias[x]][0].second;
                groups[{t.hom,mask}]+=phase_reward_scale*t.reward_units*m.weights.genotype;
                groups[{!t.hom,mask}]+=phase_reward_scale*t.other_reward_units*m.weights.genotype;
                local_upper+=phase_reward_scale*PresenceModel::target_upper(t)*m.weights.genotype;
            }
            for(auto&[key,weight]:groups)for(int z=0;z<D;++z)
                row[z]+=weight*PresenceModel::coverage(key.first,m.signature(key.second,z));
            for(size_t i:contacts[b]) {
                const auto&e=folded[i];
                if(first[e.x]!=b||first[e.y]!=b)continue;
                const Score weight=e.units*m.weights.contact;local_upper+=phase_reward_scale*weight;
                for(int z=0;z<D;++z)row[z]+=weight*phase_factor(
                    m.signature(patterns[e.x][0].second,z),m.signature(patterns[e.y][0].second,z),e.relation,e.panel);
            }
            other[b]=suffix[b]-suffix[b+1]-local_upper;
            if(other[b]<0)throw std::logic_error("Negative nonlocal completion bound");
            local.append(row);
        }
        for(int b=m.blocks-1;b>=0;--b)other[b]+=other[b+1];
        lookahead=PackedCompletion(m.H,m.blocks,m.budget,std::move(local),std::move(other));
    }
    struct Trace {size_t parent;int choice;};
    struct Live {int r,last;Score score;size_t trace;std::vector<unsigned char> values;};
    struct Row {
        int r=0;std::vector<unsigned char> before;
        std::vector<size_t> sources;
        std::vector<Score> score;std::vector<size_t> parent;
        size_t generated=0;
    };
    struct Cell {Score score;size_t parent,row;int z;Score priority=0;};
    auto better=[](const Cell&a,const Cell&b){
        if(a.priority!=b.priority)return a.priority>b.priority;
        return a.score!=b.score?a.score>b.score:std::make_tuple(a.parent,a.z,a.row)<std::make_tuple(b.parent,b.z,b.row);
    };
    std::vector<Trace> traces{{0,-1}};
    std::vector<Live> live{{0,-1,0,0,{}}};std::vector<size_t> previous;
    std::vector<DonorMask> current(N);Result result;Score omitted=impossible;
    std::vector<unsigned char> witness_presence(audit?N:0);
    Score witness_prefix=0;int witness_r=0;
    result.completion_penalty=lookahead.penalty;
    result.completion_relaxed_switches=lookahead.relaxed_switches;
    RowExecutor pool(threads);
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
        std::map<std::pair<bool,DonorMask>,Score> local_groups;
        std::map<std::pair<size_t,bool>,Score> dynamic_groups;
        for(size_t x:close[b]) {
            const auto&t=m.targets[x];
            if(first[alias[x]]==b){
                local_groups[{t.hom,current[alias[x]]}]+=t.reward_units;
                local_groups[{!t.hom,current[alias[x]]}]+=t.other_reward_units;
            }else{
                dynamic_groups[{alias[x],t.hom}]+=t.reward_units;
                dynamic_groups[{alias[x],!t.hom}]+=t.other_reward_units;
            }
        }
        struct CoverageGroup {size_t x;bool hom;Score count;};
        std::vector<CoverageGroup> dynamic_targets;
        for(const auto&[key,count]:dynamic_groups)dynamic_targets.push_back({key.first,key.second,count});
        std::vector<Score> unary(D);
        for(auto&[key,count]:local_groups)for(int z=0;z<D;++z)
            unary[z]+=phase_reward_scale*count*m.weights.genotype*PresenceModel::coverage(key.first,m.signature(key.second,z));
        std::map<std::pair<int,std::vector<unsigned char>>,size_t> row_index;
        std::vector<Row> rows;
        for(size_t i=0;i<live.size();++i) {
            auto&s=live[i];auto before=frontier_copy.gather(s.values);
            for(int delta=0;delta<=(b?2:0)&&s.r+delta<=m.budget;++delta) {
                auto key=std::make_pair(s.r+delta,before);auto it=row_index.emplace(std::move(key),rows.size());
                if(it.second) {Row row;row.r=s.r+delta;row.before=before;row.score.assign(D,impossible);row.parent.assign(D,0);rows.push_back(std::move(row));}
                rows[it.first->second].sources.push_back(i);
            }
        }
        pool.for_each(rows.size(),[&](size_t ri) {
            auto&row=rows[ri];
            for(size_t si:row.sources) {
                const auto&s=live[si];
                auto mask=[&](size_t x,int z){return static_cast<unsigned char>((old_position[x]>=0?s.values[old_position[x]]:0)|m.signature(current[x],z));};
                for(int z=0;z<D;++z) {
                    // The objective, transitions and shared budget are
                    // invariant under swapping BOTH complete haplotypes.
                    // Fix this global label symmetry only at initialization;
                    // ordering donors independently at later segments would
                    // instead erase real phase and is not permitted.
                    if(!b&&z/m.H>z%m.H)continue;
                    if(s.r+(b?m.charge(s.last,z):0)!=row.r)continue;
                    Score value=s.score+unary[z];
                    for(const auto&g:dynamic_targets)value+=phase_reward_scale*g.count*m.weights.genotype*PresenceModel::coverage(g.hom,mask(g.x,z));
                    for(size_t i:contacts[b]){const auto&e=folded[i];value+=e.units*m.weights.contact*phase_factor(mask(e.x,z),mask(e.y,z),e.relation,e.panel);}
                    for(const auto&e:phase_events[b])value+=e.value(mask(e.x,z),mask(e.y,z));
                    ++row.generated;
                    if(value>row.score[z]||(value==row.score[z]&&s.trace<row.parent[z])){row.score[z]=value;row.parent[z]=s.trace;}
                }
            }
        });
        auto make_cell=[&](size_t i,int z) {
            Cell cell{rows[i].score[z],rows[i].parent[z],i,z};
            cell.priority=cell.score+(completion&&beam?
                lookahead.remaining(b,z,m.budget-rows[i].r,suffix[b+1]):suffix[b+1]);
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
                phase_factor(witness_presence[e.x],witness_presence[e.y],e.relation,e.panel);}
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
                    if((row.before[j]|m.signature(current[active[j]],oriented))!=expected){equal=false;break;}
                }
                if(equal){auto cell=make_cell(i,oriented);if(witness_cell.z<0||better(cell,witness_cell))witness_cell=cell;}
            }
            observation.generated=witness_cell.z>=0;
            if(observation.generated){observation.best_equivalent=witness_cell.score;observation.priority=witness_cell.priority;observation.rank=1;}
        }
        // Independent reference: materialize the FULL post-emission key in
        // an ordered map, max-merge, then select the distinct-state top K.
        // No packed-key hashing, row maxima or per-shard truncation is used.
        std::map<std::tuple<int,int,std::vector<unsigned char>>,Cell> states;
        for(size_t i=0;i<rows.size();++i) {
            const auto&row=rows[i];result.generated+=row.generated;
            for(int z=0;z<D;++z)if(row.score[z]!=impossible) {
                auto after=row.before;
                for(const auto&[j,mask]:active_emissions)after[j]|=m.signature(mask,z);
                const auto cell=make_cell(i,z);
                auto [it,fresh]=states.emplace(std::make_tuple(row.r,z,std::move(after)),cell);
                if(!fresh&&better(cell,it->second))it->second=cell;
            }
        }
        std::vector<Cell> keep;
        auto discard=[&](const Cell&cell) {
            omitted=std::max(omitted,cell.priority);++result.discarded;
        };
        for(const auto&[key,cell]:states) {
                (void)key;
                if(audit){++observation.candidates;if(observation.generated&&better(cell,witness_cell))++observation.rank;}
                if(!beam){keep.push_back(cell);continue;}
                if(keep.size()<beam){keep.push_back(cell);std::push_heap(keep.begin(),keep.end(),better);}
                else if(better(cell,keep.front())) {
                    discard(keep.front());
                    std::pop_heap(keep.begin(),keep.end(),better);keep.back()=cell;std::push_heap(keep.begin(),keep.end(),better);
                } else discard(cell);
        }
        std::sort(keep.begin(),keep.end(),better);live.clear();
        for(auto&cell:keep) {
            auto&row=rows[cell.row];Live s{row.r,cell.z,cell.score,traces.size(),row.before};
            for(const auto&[i,mask]:active_emissions)s.values[i]|=m.signature(mask,cell.z);
            traces.push_back({cell.parent,cell.z});live.push_back(std::move(s));
        }
        if(audit) {
            std::set<std::tuple<int,int,std::vector<unsigned char>>> distinct;
            for(const auto&s:live) {
                distinct.emplace(s.r,s.last,s.values);
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


int main() {
    std::mt19937 rng(34171);
    {
        // C is 01 versus 10 after segment 1. Both emit 11 in segment 2 at
        // the same charge/donor pair. Keeping the two best pre-emission
        // cells loses the next distinct state; dedup after top-K is too late.
        PresenceModel m;m.H=2;m.blocks=4;m.budget=4;m.free_swaps=true;m.weights={1,1};
        DonorMask donor0;donor0.set(0);
        m.targets={
            {false,{{0,donor0}},10},
            {true,{{1,donor0},{2,donor0}},10},
            {true,{{2,donor0}},20},
            {false,{{3,donor0}},1}};
        m.contacts={{1,3,3,1,false}};
        const auto expected=serial_solve_presence(m,2);
        for(int threads:{1,2,4,8}) {
            size_t layers=0;
            PresenceSearchAudit audit{expected.path,[&](const auto&r){
                assert(r.retained==r.distinct_retained);
                if(r.completed==3)assert(r.retained==2);
                ++layers;
            }};
            const auto actual=solve_presence(m,2,threads,{},false,&audit);
            assert(layers==4&&actual.path==expected.path&&actual.score==expected.score);
        }
    }
    // Every bit alignment, empty/full words, disjoint runs and newly born
    // frontier entries must agree with the original byte representation.
    for(size_t n=0;n<1100;++n) {
        std::vector<unsigned char> bytes(n);PackedPresence packed(n);
        for(size_t i=0;i<n;++i){bytes[i]=rng()%4;packed.add(i,bytes[i]);}
        assert(packed.unpack()==bytes);
        std::vector<int> positions;
        for(size_t i=0;i<n;++i) {
            if(rng()%4)positions.push_back(int(i));
            if(rng()%8==0)positions.push_back(-1);
        }
        const FrontierCopy copy(positions);
        assert(copy.gather(bytes)==copy.gather(packed).unpack());
        assert(packed.words().size()==n/32+(n%32!=0));
    }
    for(size_t from=0;from<64;++from)for(size_t to=0;to<64;++to)
        for(size_t length:{0,1,31,32,33,63,64,65}) {
            PackedPresence source(from+length),dest(to+length+1);
            for(size_t i=0;i<source.size();++i){source.add(i,i%4);source.add(i,(i/3)%4);}
            dest.copy_zero_range(source,from,to,length);
            for(size_t i=0;i<dest.size();++i)
                assert(dest[i]==(i>=to&&i<to+length?source[from+i-to]:0));
        }
    for(int trial=0;trial<72;++trial) {
        PresenceModel model;model.H=2+rng()%4;model.blocks=4+rng()%14;model.budget=rng()%8;
        if(trial>=60){model.H=7;model.blocks=12;model.budget=20;}
        model.weights={1+Score(rng()%4),1+Score(rng()%4)};
        const size_t targets=15+rng()%30;
        for(size_t x=0;x<targets;++x) {
            std::map<int,DonorMask> occurrences;
            for(int n=0,limit=1+int(rng()%4);n<limit;++n) {
                const int b=rng()%model.blocks;
                DonorMask mask;
                for(int h=0;h<model.H;++h)if(rng()%2)mask.set(size_t(h));
                if(!mask)mask.set(size_t(rng()%model.H));
                occurrences[b]=mask;
            }
            PresenceTarget target;
            target.hom=bool(rng()%2);target.reward_units=1+rng()%7;target.other_reward_units=rng()%5;
            for(auto&entry:occurrences)target.occurrences.push_back(entry);
            model.targets.push_back(std::move(target));
        }
        for(size_t n=0;n<targets;++n) {
            const size_t x=rng()%targets,y=rng()%targets;
            if(x!=y)model.contacts.push_back({x,y,1+Score(rng()%7),rng()%2?1:-1,bool(rng()%2)});
        }
        model.check();
        for(size_t beam:{1,3,17,128}) {
            const auto expected=serial_solve_presence(model,beam,1,{},true);
            for(int threads:{1,4,8}) {
                const auto actual=solve_presence(model,beam,threads,{},true);
                assert(expected.path==actual.path&&expected.score==actual.score&&actual.upper>=expected.upper);
                assert(expected.generated==actual.generated);
                assert(expected.peak==actual.peak&&expected.max_frontier==actual.max_frontier);
                assert(expected.completion_penalty==actual.completion_penalty);
                assert(expected.completion_relaxed_switches==actual.completion_relaxed_switches);
            }
            if(trial<8&&beam==3) {
                std::vector<PresenceSearchAudit::Record> a,b;
                PresenceSearchAudit reference{expected.path,[&](const auto&r){a.push_back(r);}};
                PresenceSearchAudit measured{expected.path,[&](const auto&r){b.push_back(r);}};
                const auto before=serial_solve_presence(model,beam,1,{},true,&reference);
                const auto after=solve_presence(model,beam,4,{},true,&measured);
                assert(before.path==after.path&&a.size()==b.size());
                for(size_t i=0;i<a.size();++i) {
                    auto fields=[](const PresenceSearchAudit::Record&r){
                        return std::make_tuple(r.completed,r.recombinations,r.witness_prefix,
                            r.best_equivalent,r.priority,r.cutoff,r.rank,r.candidates,r.retained,
                            r.distinct_retained,r.generated,r.survived,r.equivalent_retained);
                    };
                    assert(fields(a[i])==fields(b[i]));
                }
            }
        }
    }
    // Donor masks spanning a machine word, many tied cells, and a beam
    // smaller than the donor-pair domain exercise hierarchical pruning.
    for(int H:{1,33,65}) {
        PresenceModel model;model.H=H;model.blocks=4;model.budget=4;model.weights={3,2};
        for(int x=0;x<8;++x) {
            DonorMask mask;for(int h=0;h<H;++h)if((h+x)%3)mask.set(size_t(h));
            if(!mask)mask.set(0);
            PresenceTarget target;target.hom=x%2;target.reward_units=3;target.other_reward_units=1;
            target.occurrences={{x%2,mask},{2+x%2,mask}};model.targets.push_back(target);
        }
        model.contacts={{0,1,7,1,false},{2,7,5,-1,true}};
        for(size_t beam:{1,33,512}) {
            const auto expected=serial_solve_presence(model,beam,1,{},true);
            const auto actual=solve_presence(model,beam,8,{},true);
            assert(expected.path==actual.path&&expected.score==actual.score&&actual.upper>=expected.upper);
            assert(expected.generated==actual.generated);
            assert(expected.peak==actual.peak&&expected.max_frontier==actual.max_frontier);
        }
    }
    // Lookup compilation must retain wide exact weights, target aliases,
    // contacts closing at different times, and every neutral/phase relation.
    // The unpruned cases also compare against the original direct scorer.
    for(int exponent:{0,65,90}) {
        PresenceModel model;model.H=2;model.blocks=3;model.budget=4;model.weights={3,2};
        const Score weight=Score(1)<<exponent;
        for(int x=0;x<6;++x) {
            PresenceTarget t;t.hom=x%2;t.reward_units=3*weight;t.other_reward_units=weight;
            DonorMask initial;initial.set(size_t(x%2));t.occurrences={{x%3,initial}};
            if(x<2){DonorMask final;final.set(size_t(1-x));t.occurrences.push_back({2,final});}
            model.targets.push_back(t);
        }
        model.targets.push_back(model.targets[0]);
        model.contacts={{0,6,7*weight,1,false},{0,1,3*weight,-1,true},
                        {2,5,weight,1,false},{3,4,2*weight,-1,false}};
        for(size_t beam:{0,1,7})for(bool completion:{false,true}) {
            const auto expected=serial_solve_presence(model,beam,1,{},completion);
            for(int threads:{1,4}) {
                const auto actual=solve_presence(model,beam,threads,{},completion);
                assert(expected.path==actual.path&&expected.score==actual.score&&actual.upper>=expected.upper);
                assert(expected.generated==actual.generated);
                assert(expected.peak==actual.peak&&expected.max_frontier==actual.max_frontier);
            }
        }
    }
    std::cout<<"PASS distinct-state scalar/parallel beam parity: emission-collapse regression, paths, rewards, conservative bounds and ties on randomized joint-contact models\n";
}
