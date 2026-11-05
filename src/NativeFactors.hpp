#pragma once
#include "BuildModel.hpp"
#include <cmath>

namespace hic::native_factors {
constexpr uint32_t scale=1024;
constexpr uint32_t pseudocount=scale/2; // one half of a physical-pair budget
inline constexpr const char* version="source-library-exact-anchor-v1";
enum class Source : unsigned char {SR,HiC,Panel};
struct Class {
    Source source=Source::Panel;
    uint32_t library=UINT32_MAX;
    MoleculeGeometry geometry=MoleculeGeometry::Unknown;
    uint8_t bin=0,qc=0;
    auto key()const{return std::tie(source,library,geometry,bin,qc);}
    bool operator<(const Class&b)const{return key()<b.key();}
    bool operator==(const Class&b)const{return key()==b.key();}
};
inline uint8_t distance_bin(uint64_t d) {
    constexpr uint64_t edges[]={300,1000,3000,10000,30000,100000,300000,1000000,3000000};
    for(uint8_t i=0;i<9;++i)if(d<edges[i])return i;
    return 9;
}
inline Class classify(const MoleculeMarkers&m) {
    if((m.geometry!=MoleculeGeometry::Known&&m.geometry!=MoleculeGeometry::Unknown&&m.geometry!=MoleculeGeometry::Trans)||m.library==UINT32_MAX||
       (m.localization_qc!=1&&m.localization_qc!=2))throw std::invalid_argument("Missing native molecule class");
    if((m.geometry==MoleculeGeometry::Known&&m.reference_distance<0)||
       (m.geometry!=MoleculeGeometry::Known&&m.reference_distance!=unknown_contact_distance))
        throw std::invalid_argument("Inconsistent native geometry");
    return {m.kind==MoleculeKind::ShortRead?Source::SR:Source::HiC,m.library,m.geometry,
            m.geometry==MoleculeGeometry::Known?distance_bin(uint64_t(m.reference_distance)):uint8_t(0),m.localization_qc};
}
inline std::string name(const Class&c) {
    if(c.source==Source::Panel)return "panel";
    return std::string(c.source==Source::SR?"SR":"HIC")+"/library="+std::to_string(c.library)+"/"+
        (c.geometry==MoleculeGeometry::Known?"anchor_bin="+std::to_string(c.bin):
         c.geometry==MoleculeGeometry::Trans?"trans":"unknown")+"/qc="+std::to_string(c.qc);
}
struct ClassFactors {Class identity;std::vector<PresenceContact> factors;};
struct Snapshot {
    std::vector<ClassFactors> classes;
    size_t molecules=0;
    Score raw_units=0,read_units=0;
};
// Classify BEFORE pooling. Each physical pair has exactly one class shared by
// all its eligible target pairs. Allocation and the outer het discount retain
// the baseline order: allocate 1024 units, sum compatible observations, then
// floor(pooled * t_x * t_y / 1024^2). No source confidence has been applied yet.
inline Snapshot prepare(const PreparedModel&prepared,const SequenceIndex&index,const RawEvidence&raw,
                        bool stage2,bool het_markers,int threads) {
    if(threads<1)throw std::invalid_argument("Native factor workers");
    const auto&m=prepared.model;
    if(prepared.target_ids.size()!=m.targets.size())throw std::invalid_argument("Native target remap");
    const auto&remap=prepared.linkage_remap;
    if(remap.size()!=index.targets.size())throw std::invalid_argument("Native linkage remap missing");
    std::vector<Score> het(m.targets.size());
    for(size_t x=0;x<m.targets.size();++x) {
        const auto&t=m.targets[x];het[x]=t.hom?t.other_reward_units:t.reward_units;
        if(het[x]<0||het[x]>scale)throw std::invalid_argument("Native heterozygosity units");
    }
    struct Key {
        Class identity;size_t x,y;
        bool operator==(const Key&b)const{return identity==b.identity&&x==b.x&&y==b.y;}
        bool operator<(const Key&b)const{return std::tie(identity,x,y)<std::tie(b.identity,b.x,b.y);}
    };
    struct Hash {
        size_t operator()(const Key&k)const {
            uint64_t h=k.x*0x9e3779b97f4a7c15ULL+k.y;
            h^=uint64_t(k.identity.library)*0xbf58476d1ce4e5b9ULL;
            h^=(uint64_t(k.identity.source)<<24)|(uint64_t(k.identity.geometry)<<16)|(uint64_t(k.identity.bin)<<8)|k.identity.qc;
            h^=h>>30;h*=0xbf58476d1ce4e5b9ULL;h^=h>>27;return h;
        }
    };
    struct Contribution {Key key;Score units;};
    struct Task {std::vector<std::vector<Contribution>> bins;size_t molecules=0;};
    RowExecutor pool(threads);const size_t owners=std::min(size_t(256),size_t(threads)*4);
    std::vector<FlatMap<Key,Score,Hash>> sums(owners);
    Snapshot result;constexpr size_t task_size=64;
    const size_t batch_size=owners*task_size;
    std::vector<Task> tasks(owners);
    for(auto&task:tasks)task.bins.resize(owners);
    for(size_t begin=0;begin<raw.markers.size();begin+=batch_size) {
        const size_t end=std::min(raw.markers.size(),begin+batch_size);
        tasks.resize((end-begin+task_size-1)/task_size);
        pool.for_each(tasks.size(),[&](size_t task) {
            auto&out=tasks[task];out.molecules=0;
            for(auto&bin:out.bins)bin.clear();
            for(size_t i=begin+task*task_size;i<std::min(end,begin+(task+1)*task_size);++i) {
                const auto&original=raw.markers[i];const auto c=classify(original);
                MoleculeMarkers filtered;filtered.id=original.id;
                for(int e=0;e<2;++e)for(size_t x:original.ends[e]) {
                    const size_t target=remap.at(x);if(target==SIZE_MAX)continue;
                    const auto&t=m.targets[target];const auto hom=t.hom?t.reward_units:t.other_reward_units;
                    if(het_markers&&(!het[target]||het[target]<hom||het[target]<t.absent_reward_units))continue;
                    filtered.ends[e].push_back(x);
                }
                if(molecule_pair_units(filtered,remap,[&](const auto&pair,Score units) {
                    Key key{c,pair.first,pair.second};out.bins[Hash{}(key)%owners].push_back({key,units});
                },stage2?&het:nullptr))++out.molecules;
            }
            for(auto&bin:out.bins)std::sort(bin.begin(),bin.end(),[](const auto&a,const auto&b){return a.key<b.key;});
        });
        pool.for_each(owners,[&](size_t owner) {
            auto&dest=sums[owner];
            for(const auto&task:tasks) {
                const auto&bin=task.bins[owner];
                for(size_t a=0;a<bin.size();) {
                    size_t b=a+1;Score units=bin[a].units;
                    while(b<bin.size()&&bin[b].key==bin[a].key)
                        units=score_sum(units,bin[b++].units);
                    auto&value=dest[bin[a].key];
                    value=score_sum(value,units);a=b;
                }
            }
        });
        for(const auto&task:tasks) {
            if(task.molecules>SIZE_MAX-result.molecules)throw std::overflow_error("Native molecule count");
            result.molecules+=task.molecules;
        }
    }
    result.raw_units=score_product(Score(result.molecules),scale);
    std::vector<Contribution> merged;
    for(auto&shard:sums) {
        for(const auto&[key,units]:shard) {
            const Score weighted=score_product(score_product(units,het[key.x]),het[key.y])/(Score(scale)*scale);
            if(weighted)merged.push_back({key,weighted});
        }
        std::decay_t<decltype(shard)>{}.swap(shard);
    }
    std::sort(merged.begin(),merged.end(),[](const auto&a,const auto&b){return a.key<b.key;});
    for(const auto&e:merged) {
        if(result.classes.empty()||!(result.classes.back().identity==e.key.identity))result.classes.push_back({e.key.identity,{}});
        result.classes.back().factors.push_back({e.key.x,e.key.y,e.units,1,false});
        result.read_units=score_sum(result.read_units,e.units);
    }
    if(result.read_units>result.raw_units)throw std::logic_error("Native factor mass exceeds physical observations");
    ClassFactors panel;
    for(const auto&e:m.contacts)if(e.panel)panel.factors.push_back(e);
    if(!panel.factors.empty())result.classes.push_back(std::move(panel));
    return result;
}
struct Estimate {
    struct Entry {Score agree=0,oppose=0;long double epsilon=.5L,llr=0;uint32_t confidence=0;};
    std::map<Class,Entry> classes;
    long double maximum_llr=0;
};
// Estimate in pre-confidence factor units (1024 = one full physical-pair
// budget), with a fixed 512-unit pseudocount on each side. Panel participates
// in the SAME maximum as SR and Hi-C. Sparse/unseen classes are not pooled:
// zero informative evidence means epsilon=1/2 and confidence zero, including
// the source-specific unknown-distance classes. Freeze this result for a stage.
inline Estimate estimate(const Snapshot&snapshot,const PresenceModel&m,const std::vector<int>&path) {
    if(path.size()!=size_t(m.blocks))throw std::invalid_argument("Native calibration path length");
    for(int z:path)if(z<0||z>=m.D())throw std::invalid_argument("Native calibration donor pair");
    std::vector<unsigned char> presence(m.targets.size());
    for(size_t x=0;x<m.targets.size();++x)for(const auto&[b,mask]:m.targets[x].occurrences)
        presence[x]|=m.signature(mask,path.at(size_t(b)));
    Estimate result;
    for(const auto&group:snapshot.classes) {
        auto&c=result.classes[group.identity];
        for(const auto&e:group.factors) {
            const int orientation=e.relation*phase_orientation(presence.at(e.x))*phase_orientation(presence.at(e.y));
            if(orientation>0)c.agree=score_sum(c.agree,e.units);
            if(orientation<0)c.oppose=score_sum(c.oppose,e.units);
        }
    }
    for(auto&[identity,c]:result.classes) {
        (void)identity;
        const long double a=static_cast<long double>(c.agree)+pseudocount;
        const long double o=static_cast<long double>(c.oppose)+pseudocount;
        c.epsilon=o/(a+o);c.llr=std::max(0.L,std::log(a)-std::log(o));
        if(!std::isfinite(c.llr)||!std::isfinite(c.epsilon))throw std::overflow_error("Native confidence estimate");
        result.maximum_llr=std::max(result.maximum_llr,c.llr);
    }
    for(auto&[identity,c]:result.classes) {
        (void)identity;
        c.confidence=result.maximum_llr>0?uint32_t(std::llround(scale*c.llr/result.maximum_llr)):0;
        if(c.confidence>scale)throw std::logic_error("Native confidence outside unit interval");
    }
    return result;
}
// Apply immutable class confidence BEFORE compatible factors are pooled.
// Keep the common confidence denominator instead of rounding each small
// factor down to zero. With an estimate, stored units are u*c and denote
// u*c/scale pre-confidence units. W has the same denominator, which cancels
// exactly in the normalized objective. A null estimate retains units u.
// Opposite signs and panel/read provenance remain distinct.
inline Score apply(const Snapshot&snapshot,PresenceModel&m,const Estimate*estimate=nullptr) {
    using Key=std::tuple<size_t,size_t,int,bool>;
    std::map<Key,Score> pooled;
    for(const auto&group:snapshot.classes) {
        uint32_t confidence=scale;
        if(estimate) {
            auto found=estimate->classes.find(group.identity);
            confidence=found==estimate->classes.end()?0:found->second.confidence;
            if(confidence>scale)throw std::invalid_argument("Native confidence exceeds scale");
        }
        for(const auto&e:group.factors) {
            if(e.units<0||e.x>=m.targets.size()||e.y>=m.targets.size()||e.x==e.y||(e.relation!=1&&e.relation!=-1))
                throw std::invalid_argument("Invalid native factor");
            const Score units=estimate?score_product(e.units,confidence):e.units;
            if(units){auto&slot=pooled[{e.x,e.y,e.relation,e.panel}];slot=score_sum(slot,units);}
        }
    }
    m.contacts.clear();m.contacts.reserve(pooled.size());Score total=0;
    for(const auto&[key,units]:pooled) {
        const auto[x,y,relation,panel]=key;m.contacts.push_back({x,y,units,relation,panel});
        total=score_sum(total,units);
    }
    return total;
}
// Rebuild normalization after EVERY snapshot change. Bounds and incumbents
// must then be recomputed by the caller; no score from the old snapshot is a
// valid certificate for these weights.
inline void normalize(PreparedModel&p,Score phase_weight=1) {
    p.read_contact_units=0;p.panel.units=0;
    for(const auto&e:p.model.contacts) {
        auto&units=e.panel?p.panel.units:p.read_contact_units;
        units=score_sum(units,e.units);
    }
    p.contact_units=score_sum(p.read_contact_units,p.panel.units);
    p.model.weights=normalized_weights(p.contact_units,p.het_reward_units,phase_weight);
    p.model.check();
}
inline void write_snapshot(std::ostream&out,const Snapshot&s) {
    out<<"class\ttarget1\ttarget2\tweight_units\trelation\tsource\n";
    for(const auto&g:s.classes)for(const auto&e:g.factors)
        out<<name(g.identity)<<'\t'<<e.x<<'\t'<<e.y<<'\t'<<e.units<<'\t'<<e.relation<<'\t'<<(e.panel?"panel":"molecule")<<'\n';
    if(!out)throw std::runtime_error("Native factor snapshot write failed");
}
} // namespace hic::native_factors
