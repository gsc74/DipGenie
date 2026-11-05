// Positional seeding, exact extension, interval-union voting and diagonal
// reduction adapted from hifiasm hic.cpp at
// ec9a8b222d149d25b7355e83765698640d59b189. No assembler-owned objects or packed
// coordinates are imported. See third_party/HIFIASM-LICENSE.txt (MIT).
#pragma once
#include "SequenceEvidence.hpp"
#include "ExactMatch.hpp"
#include <tuple>

namespace hic {

// An oriented graph interval spelled in a continuous context. `end` is its
// end-exclusive context coordinate; offset is the lower forward-node offset.
// copy=0 denotes a node visited at most once in each panel path. For repeated
// visits the adapter supplies an explicit, conservative occurrence identity.
struct PlacementRun {
    uint64_t end=0,node=0,copy=0,offset=0;
    bool reverse=false;
};
struct PlacementSlice {
    uint64_t node=0,copy=0,begin=0,end=0;
    bool reverse=false;
    auto key()const{return std::tie(node,copy,begin,end,reverse);}
    bool operator<(const PlacementSlice&x)const{return key()<x.key();}
    bool operator==(const PlacementSlice&x)const{return key()==x.key();}
};
// A continuous, oriented sequence context. Copies/chunks of the SAME physical
// coordinate space share space and origin. Distinct repeats MUST use distinct
// spaces or coordinates. A graph adapter must supply this identity; equal DNA
// alone is never evidence that two genomic loci are the same.
// Coordinates increase along the physical space, even for reversed contexts.
struct PlacementSequence {
    std::string name,sequence;
    uint32_t component=0,space=0;
    uint64_t origin=0;
    bool reverse=false;
    uint32_t reference_contig=UINT32_MAX; // only designated reference contexts
    std::vector<PlacementRun> runs; // optional: piecewise graph-coordinate map
};
inline void reverse_placement_slices(std::vector<PlacementSlice>&v) {
    std::reverse(v.begin(),v.end());for(auto&x:v)x.reverse=!x.reverse;
}
inline std::vector<PlacementSlice> placement_slices(const PlacementSequence&s,uint64_t begin,uint64_t end) {
    if(begin>=end||end>s.sequence.size()||s.runs.empty())throw std::invalid_argument("Graph placement interval");
    auto it=std::upper_bound(s.runs.begin(),s.runs.end(),begin,[](uint64_t p,const PlacementRun&r){return p<r.end;});
    std::vector<PlacementSlice> result;
    while(begin<end) {
        if(it==s.runs.end())throw std::logic_error("Uncovered graph placement interval");
        const uint64_t run_begin=it==s.runs.begin()?0:(it-1)->end;
        const uint64_t finish=std::min(end,it->end),a=begin-run_begin,b=finish-run_begin;
        const uint64_t lo=it->offset+(it->reverse?it->end-run_begin-b:a);
        const uint64_t hi=it->offset+(it->reverse?it->end-run_begin-a:b);
        PlacementSlice slice{it->node,it->copy,lo,hi,it->reverse};
        // A chunk boundary inside the same physical interval is not a new locus.
        if(!result.empty()&&result.back().node==slice.node&&result.back().copy==slice.copy&&
           result.back().reverse==slice.reverse&&
           (slice.reverse?result.back().begin==slice.end:result.back().end==slice.begin)) {
            result.back().begin=std::min(result.back().begin,slice.begin);
            result.back().end=std::max(result.back().end,slice.end);
        } else result.push_back(slice);
        begin=finish;++it;
    }
    return result;
}
struct MatchInterval {
    uint32_t begin=0,end=0; // half-open, in the original read
    bool operator==(const MatchInterval&x)const{return begin==x.begin&&end==x.end;}
};
enum class PlacementStatus {NoSeeds,InsufficientSupport,Unique,Ambiguous,Split};
enum class AnchorStatus {Known,NoSeeds,InsufficientSupport,Ambiguous,Split,Unprojected};
struct SelectedAnchor {
    AnchorStatus status=AnchorStatus::NoSeeds;
    uint32_t contig=UINT32_MAX;
    uint64_t start=0,word=0; // exact reference k-mer START, never a read start
};
struct PlacementSeed {
    uint64_t word=0,position=0;
    uint32_t query=0,context=0;
    bool reverse=false;
};
struct PlacementCandidate {
    uint32_t component=0,space=0,support=0;
    bool reverse=false;
    int64_t diagonal_low=0,diagonal_high=0;
    std::vector<MatchInterval> intervals;
    std::vector<PlacementSeed> seeds;
    SelectedAnchor anchor;
    // Canonical-query-oriented physical context, for max-merging aliases from
    // donor paths. Distinct allele routes or repeated visits are not merged.
    std::vector<PlacementSlice> graph_context;
    uint32_t context_clip_left=0,context_clip_right=0;
    int64_t context_jitter=0;
};
inline bool same_placement_alias(const PlacementCandidate&a,const PlacementCandidate&b) {
    return a.component==b.component&&a.context_clip_left==b.context_clip_left&&
        a.context_clip_right==b.context_clip_right&&a.context_jitter==b.context_jitter&&
        a.graph_context==b.graph_context;
}
struct PlacementAliasHash {
    size_t operator()(const PlacementCandidate&c)const {
        uint64_t h=0x9e3779b97f4a7c15ULL;
        auto add=[&](uint64_t x) {
            x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;
            x*=0x94d049bb133111ebULL;x^=x>>31;
            h^=x+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);
        };
        add(c.component);add(c.context_clip_left);add(c.context_clip_right);
        add(uint64_t(c.context_jitter));
        for(const auto&s:c.graph_context) {
            add(s.node);add(s.copy);add(s.begin);add(s.end);add(s.reverse);
        }
        return size_t(h);
    }
};
// Only the already-owned candidate stores its graph footprint. Hashes locate
// a collision chain; full equality still decides aliasing. Preserve input
// order, maximum support and the original selected-anchor tie rule exactly.
template<class Hash=PlacementAliasHash>
inline std::vector<PlacementCandidate> merge_placement_aliases(
        std::vector<PlacementCandidate> candidates,Hash hash={}) {
    FlatMap<size_t,size_t> heads;
    heads.reserve(candidates.size());
    std::vector<size_t> next;next.reserve(candidates.size());
    std::vector<PlacementCandidate> distinct;distinct.reserve(candidates.size());
    for(auto&c:candidates) {
        if(c.graph_context.empty()) {
            distinct.push_back(std::move(c));next.push_back(SIZE_MAX);continue;
        }
        auto [head,inserted]=heads.emplace(hash(c),SIZE_MAX);
        (void)inserted;
        size_t at=head->second;
        while(at!=SIZE_MAX&&!same_placement_alias(c,distinct[at]))at=next[at];
        if(at==SIZE_MAX) {
            next.push_back(head->second);head->second=distinct.size();
            distinct.push_back(std::move(c));continue;
        }
        auto&kept=distinct[at];
        if(c.support>kept.support)kept=std::move(c);
        else if(c.support==kept.support&&c.anchor.status==AnchorStatus::Known&&
                (kept.anchor.status!=AnchorStatus::Known||
                 std::tie(c.anchor.start,c.anchor.word)<std::tie(kept.anchor.start,kept.anchor.word)))
            kept.anchor=c.anchor;
    }
    return distinct;
}
struct EndPlacement {
    PlacementStatus status=PlacementStatus::NoSeeds;
    uint32_t component=UINT32_MAX;
    uint32_t best_support=0,competing_support=0;
    uint64_t skipped_repetitive_seeds=0;
    std::vector<PlacementCandidate> candidates;
    std::vector<size_t> retained;
    SelectedAnchor anchor;
    // Component-unique is deliberately weaker than a unique locus/anchor.
    bool component_unique()const{return component!=UINT32_MAX;}
};
enum class AnchorDistanceStatus {Known,Unknown,Trans};
struct AnchorDistance {
    AnchorDistanceStatus status=AnchorDistanceStatus::Unknown;
    uint64_t separation=0; // meaningful ONLY when status == Known
    AnchorStatus ends[2]={AnchorStatus::NoSeeds,AnchorStatus::NoSeeds};
};
inline AnchorDistance selected_anchor_distance(const EndPlacement&a,const EndPlacement&b) {
    AnchorDistance d;d.ends[0]=a.anchor.status;d.ends[1]=b.anchor.status;
    if(d.ends[0]!=AnchorStatus::Known||d.ends[1]!=AnchorStatus::Known)return d;
    if(a.anchor.contig!=b.anchor.contig){d.status=AnchorDistanceStatus::Trans;return d;}
    d.status=AnchorDistanceStatus::Known;
    d.separation=a.anchor.start>b.anchor.start?a.anchor.start-b.anchor.start:b.anchor.start-a.anchor.start;
    return d;
}

// Preserve selected identities but recover all their occurrences before
// interval filtering. The first clipped copy cannot hide an admitted copy.
// This routine never modifies raw exact per-end dosage or occurrence counts.
inline std::vector<size_t> native_interval_markers(const SequenceIndex&index,
    const std::string&read,const std::vector<MatchInterval>&intervals) {
    std::unordered_set<uint64_t> selected;
    minimizers(read,index.k,index.w,[&](uint64_t x,size_t){selected.insert(x);});
    std::vector<size_t> markers;
    words(read,index.k,[&](uint64_t word,size_t p,bool){
        if(!selected.count(word))return;
        auto it=index.ids.find(word);
        if(it==index.ids.end()||!index.targets[it->second].unique_locus)return;
        for(const auto&range:intervals)if(p>=range.begin&&p+index.k<=range.end){markers.push_back(it->second);break;}
    });
    std::sort(markers.begin(),markers.end());markers.erase(std::unique(markers.begin(),markers.end()),markers.end());
    return markers;
}

class NativeLocalizer {
public:
    struct Policy {
        unsigned k=31,w=25;
        uint32_t maximum_loci=64,minimum_seeds=2;
        // Integer percentages, fixed before inference; not calibrated MAPQ.
        uint32_t minimum_support_percent=50,competitor_percent=80;
        uint32_t maximum_read_length=1024;
    };
    static constexpr const char* version="hifiasm-locus-longest-interval-v2";
    // Index construction changes storage/scheduling only, not evidence policy.
    // The bound covers batch hit vectors, not the persistent context/index or
    // per-word physical-locus sorting. Their sizes are reported separately.
    static constexpr size_t ingest_scratch_limit=256u<<20;
    struct IndexStats {
        uint64_t context_bases=0,selected_words=0,retained_occurrences=0;
        uint64_t physical_loci=0,discarded_words=0;
        size_t peak_ingest_scratch_bytes=0,ingest_batches=0;
    };
private:
    struct Occurrence {
        // Pack the new locus ID into the old 24-byte footprint: chromosome
        // indexes contain billions of donor occurrences.
        uint64_t position=0;
        uint32_t context=0,projection=UINT32_MAX,locus=0; // locus is local to this word
        bool reverse=false;
        Occurrence()=default;
        Occurrence(uint32_t c,uint64_t p,bool r):position(p),context(c),reverse(r){}
    };
    static_assert(sizeof(Occurrence)==24,"Native occurrence storage grew");
    struct Word {std::vector<Occurrence> occurrences;uint64_t loci=0;std::vector<SelectedAnchor> projections;};
    struct Hit {
        uint32_t component,space;bool reverse;
        int64_t diagonal;MatchInterval interval;PlacementSeed seed;
        SelectedAnchor anchor;
    };
    Policy policy_;
    std::vector<PlacementSequence> contexts_;
    FlatMap<uint64_t,size_t> index_;
    std::vector<Word> words_;
    IndexStats stats_;

    uint64_t physical_start(const Occurrence&o)const {
        const auto&s=contexts_[o.context];
        return s.origin+(s.reverse?s.sequence.size()-o.position-policy_.k:o.position);
    }
    static uint32_t union_intervals(std::vector<MatchInterval>&v) {
        std::sort(v.begin(),v.end(),[](auto a,auto b){return std::tie(a.begin,a.end)<std::tie(b.begin,b.end);});
        size_t n=0;uint64_t support=0;
        for(auto x:v) {
            if(n&&x.begin<=v[n-1].end)v[n-1].end=std::max(v[n-1].end,x.end);
            else v[n++]=x;
        }
        v.resize(n);for(auto x:v)support+=x.end-x.begin;
        if(support>UINT32_MAX)throw std::overflow_error("Native query support");
        return uint32_t(support);
    }
    static uint32_t overlap(const std::vector<MatchInterval>&a,const std::vector<MatchInterval>&b) {
        size_t i=0,j=0;uint32_t n=0;
        while(i<a.size()&&j<b.size()) {
            const auto lo=std::max(a[i].begin,b[j].begin),hi=std::min(a[i].end,b[j].end);
            if(lo<hi)n+=hi-lo;
            if(a[i].end<b[j].end)++i;else ++j;
        }
        return n;
    }
public:
    // Immutable after construction: only local/per-thread scratch is used by
    // locate(). Continuous full sequences avoid any node/chunk boundary reset.
    NativeLocalizer(std::vector<PlacementSequence> contexts,Policy policy,int threads=1)
        :policy_(policy),contexts_(std::move(contexts)) {
        if(policy.k<1||policy.k>31||!policy.w||!policy.maximum_loci||!policy.minimum_seeds||threads<1||
           policy.minimum_support_percent>100||policy.competitor_percent>100||!policy.competitor_percent||
           policy.maximum_read_length<policy.k||contexts_.size()>UINT32_MAX)
            throw std::invalid_argument("Native localization policy");
        std::map<uint32_t,uint32_t> owners,reference_owners;
        std::map<uint32_t,bool> graph_spaces;
        for(auto&s:contexts_) {
            if(s.origin>uint64_t(INT64_MAX)-policy.maximum_read_length||
               s.sequence.size()>uint64_t(INT64_MAX)-policy.maximum_read_length-s.origin)
                throw std::overflow_error("Native physical coordinate");
            if(s.sequence.size()>UINT64_MAX-stats_.context_bases)
                throw std::overflow_error("Native total context length");
            stats_.context_bases+=s.sequence.size();
            auto [it,inserted]=owners.emplace(s.space,s.component);
            if(!inserted&&it->second!=s.component)throw std::invalid_argument("Physical space has multiple components");
            const auto [space,added_space]=graph_spaces.emplace(s.space,!s.runs.empty());
            if(!added_space&&(space->second||!s.runs.empty()))
                throw std::invalid_argument("Graph contexts require separate spelled coordinate spaces");
            if(s.reference_contig!=UINT32_MAX) {
                auto [ref,added]=reference_owners.emplace(s.space,s.reference_contig);
                if(!added&&ref->second!=s.reference_contig)throw std::invalid_argument("Physical space has multiple reference contigs");
            }
            for(char&c:s.sequence)c=char(std::toupper(static_cast<unsigned char>(c)));
            if(!s.runs.empty()) {
                // Graph contexts cluster in their own spelled coordinates;
                // orientation lives in the explicit runs, never in a guessed
                // affine projection through insertions/deletions.
                if(s.reverse||s.origin)throw std::invalid_argument("Graph context must use its own spelled coordinates");
                uint64_t previous=0;
                for(const auto&r:s.runs) {
                    if(r.end<=previous||r.end>s.sequence.size()||r.offset>UINT64_MAX-(r.end-previous))
                        throw std::invalid_argument("Malformed graph context runs");
                    previous=r.end;
                }
                if(previous!=s.sequence.size())throw std::invalid_argument("Incomplete graph context runs");
            }
        }
        RowExecutor pool(threads);
        struct Chunk {uint32_t context;size_t begin,count;};
        constexpr size_t chunk_bases=1u<<16;
        auto chunks=[&](size_t right_context) {
            std::vector<Chunk> tasks;
            for(size_t c=0;c<contexts_.size();++c) {
                const size_t length=contexts_[c].sequence.size();
                const size_t count=length>right_context?length-right_context:0;
                for(size_t begin=0;begin<count;begin+=chunk_bases)
                    tasks.push_back({uint32_t(c),begin,std::min(chunk_bases,count-begin)});
            }
            return tasks;
        };
        // Identity discovery then ALL occurrences: repeated equal minimizers
        // must not disappear before multiplicity filtering (e.g. poly-A).
        const auto windows=chunks(size_t(policy.k)+policy.w-2);
        const size_t discovery_batch=std::min(size_t(threads)*2,
            ingest_scratch_limit/(sizeof(uint64_t)*chunk_bases));
        for(size_t first=0;first<windows.size();first+=discovery_batch) {
            std::vector<std::vector<uint64_t>> selected(std::min(discovery_batch,windows.size()-first));
            pool.for_each(selected.size(),[&](size_t i){
                const auto&t=windows[first+i];auto&v=selected[i];v.reserve(t.count);
                minimizers(std::string_view(contexts_[t.context].sequence).substr(t.begin,t.count+policy.k+policy.w-2),
                           policy.k,policy.w,[&](uint64_t x,size_t){v.push_back(x);});
                sort_words(v);v.erase(std::unique(v.begin(),v.end()),v.end());
            });
            size_t scratch=0;
            for(const auto&v:selected){scratch+=v.capacity()*sizeof(uint64_t);for(auto x:v)index_.try_emplace(x,index_.size());}
            stats_.peak_ingest_scratch_bytes=std::max(stats_.peak_ingest_scratch_bytes,scratch);
        }
        words_.resize(index_.size());stats_.selected_words=index_.size();
        if(words_.empty())return;
        const auto starts=chunks(policy.k-1);
        struct Record {size_t id;Occurrence occurrence;};
        struct BatchChunk {std::vector<Record> data;std::vector<size_t> offsets;};
        const size_t batch_limit=std::min(size_t(threads)*2,
            ingest_scratch_limit/(2*sizeof(Record)*chunk_bases));
        const size_t shards=std::min(words_.size(),size_t(threads)*4);
        const size_t per_shard=(words_.size()+shards-1)/shards;
        const auto&lookup=index_;
        for(size_t first=0;first<starts.size();first+=batch_limit) {
            std::vector<BatchChunk> batch(std::min(batch_limit,starts.size()-first));
            pool.for_each(batch.size(),[&](size_t i){
                const auto&t=starts[first+i];std::vector<Record> found;found.reserve(t.count);
                words(std::string_view(contexts_[t.context].sequence).substr(t.begin,t.count+policy.k-1),policy.k,[&](uint64_t x,size_t p,bool reverse){
                    const auto at=lookup.find(x);
                    if(at!=lookup.end())found.push_back({at->second,{t.context,t.begin+p,reverse}});
                });
                auto&dest=batch[i];dest.offsets.assign(shards+1,0);
                for(const auto&r:found)++dest.offsets[r.id/per_shard+1];
                for(size_t s=0;s<shards;++s)dest.offsets[s+1]+=dest.offsets[s];
                auto cursor=dest.offsets;dest.data.resize(found.size());
                for(const auto&r:found)dest.data[cursor[r.id/per_shard]++]=r;
            });
            // Disjoint word owners: no per-hit lock/atomic or serial merge.
            // Batch/chunk order preserves the previous insertion order.
            pool.for_each(shards,[&](size_t owner){
                for(const auto&part:batch)for(size_t i=part.offsets[owner];i<part.offsets[owner+1];++i) {
                    const auto&r=part.data[i];words_[r.id].occurrences.push_back(r.occurrence);
                }
            });
            ++stats_.ingest_batches;
            // Conservative simultaneous found+partitioned-vector capacity.
            size_t scratch=0;
            for(size_t i=0;i<batch.size();++i)
                scratch+=(starts[first+i].count+batch[i].data.capacity())*sizeof(Record);
            stats_.peak_ingest_scratch_bytes=std::max(stats_.peak_ingest_scratch_bytes,scratch);
        }
        // Count physical loci, not donor copies or overlapping index chunks.
        // Keep the contexts for extension; union voting prevents duplicate mass.
        pool.for_each(words_.size(),[&](size_t i) {
            auto&w=words_[i];
            struct Key {
                bool graph=false;uint32_t component=0,space=0;uint64_t start=0;
                std::vector<PlacementSlice> slices;
                auto tuple()const{return std::tie(graph,component,space,start,slices);}
                bool operator<(const Key&b)const{return tuple()<b.tuple();}
                bool operator==(const Key&b)const{return tuple()==b.tuple();}
            };
            struct Decorated {Occurrence occurrence;Key key;};
            std::vector<Decorated> ordered;ordered.reserve(w.occurrences.size());
            for(const auto&o:w.occurrences) {
                const auto&s=contexts_[o.context];Key key;key.component=s.component;
                if(s.runs.empty()){key.space=s.space;key.start=physical_start(o);}
                else {
                    key.graph=true;key.slices=placement_slices(s,o.position,o.position+policy.k);
                    auto reversed=key.slices;reverse_placement_slices(reversed);
                    if(reversed<key.slices)key.slices=std::move(reversed);
                }
                ordered.push_back({o,std::move(key)});
            }
            std::sort(ordered.begin(),ordered.end(),[](const auto&a,const auto&b){
                if(!(a.key==b.key))return a.key<b.key;
                return std::tie(a.occurrence.context,a.occurrence.position)<std::tie(b.occurrence.context,b.occurrence.position);
            });
            for(size_t a=0;a<ordered.size();) {
                size_t b=a+1;while(b<ordered.size()&&ordered[b].key==ordered[a].key)++b;
                ++w.loci;
                if(w.loci>UINT32_MAX)throw std::overflow_error("Native word locus identity");
                SelectedAnchor projection;projection.status=AnchorStatus::Unprojected;
                bool ambiguous=false;
                for(size_t j=a;j<b;++j) {
                    const auto&o=ordered[j].occurrence;const auto&s=contexts_[o.context];
                    if(s.reference_contig==UINT32_MAX)continue;
                    const uint64_t p=physical_start(o);
                    if(projection.status==AnchorStatus::Known&&
                       (projection.contig!=s.reference_contig||projection.start!=p))ambiguous=true;
                    projection={AnchorStatus::Known,s.reference_contig,p,0};
                }
                uint32_t id=UINT32_MAX;
                if(projection.status==AnchorStatus::Known&&!ambiguous) {
                    if(w.projections.size()>=UINT32_MAX)throw std::overflow_error("Native reference projections");
                    id=uint32_t(w.projections.size());w.projections.push_back(projection);
                }
                for(size_t j=a;j<b;++j) {
                    ordered[j].occurrence.projection=id;
                    ordered[j].occurrence.locus=uint32_t(w.loci-1);
                    w.occurrences[j]=ordered[j].occurrence;
                }
                a=b;
            }
            if(w.loci>policy.maximum_loci){std::vector<Occurrence>().swap(w.occurrences);std::vector<SelectedAnchor>().swap(w.projections);}
        });
        for(const auto&w:words_) {
            if(w.loci>UINT64_MAX-stats_.physical_loci||w.occurrences.size()>UINT64_MAX-stats_.retained_occurrences)
                throw std::overflow_error("Native index totals");
            stats_.physical_loci+=w.loci;stats_.retained_occurrences+=w.occurrences.size();
            stats_.discarded_words+=w.loci>policy.maximum_loci;
        }
    }
    const Policy&policy()const{return policy_;}
    const IndexStats&stats()const{return stats_;}
    const std::vector<PlacementSequence>&contexts()const{return contexts_;}
    uint64_t physical_loci(uint64_t word)const {
        auto it=index_.find(word);return it==index_.end()?0:words_[it->second].loci;
    }
    EndPlacement locate(std::string_view original)const {
        if(original.size()>policy_.maximum_read_length)throw std::invalid_argument("Read exceeds native context policy");
        EndPlacement out;
        if(original.size()<policy_.k)return out;
        std::string query(original);for(char&c:query)c=char(std::toupper(static_cast<unsigned char>(c)));
        auto rc=reverse_complement(query);const bool reversed=rc<query;
        if(reversed)query=std::move(rc);
        const uint32_t length=uint32_t(query.size());
        std::vector<Hit> hits;
        // Multiple seeds often lie in the same maximal exact interval. Reuse
        // only the sequence-comparison result, not the seed/locus/anchor vote.
        // Context identity and q->reference mapping must both agree; physical
        // donor aliases and repeated loci are NOT interchangeable cache keys.
        struct ExtensionKey {
            uint32_t context;int64_t offset;bool reverse;
            bool operator==(const ExtensionKey&b)const {
                return context==b.context&&offset==b.offset&&reverse==b.reverse;
            }
        };
        struct ExtensionHash {
            size_t operator()(const ExtensionKey&k)const {
                uint64_t x=uint64_t(k.offset)^(uint64_t(k.context)*0x9e3779b97f4a7c15ULL);
                x^=uint64_t(k.reverse)*0xbf58476d1ce4e5b9ULL;
                x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;
                return size_t(x^(x>>31));
            }
        };
        FlatMap<ExtensionKey,MatchInterval,ExtensionHash> extensions;
        extensions.reserve(std::min(contexts_.size(),size_t(2048))*2);
        // Scan occurrences, not deduplicated identities. Canonicalizing the
        // complete read makes cluster boundaries and all tie decisions RC
        // invariant even when two placements have exactly equal support.
        words(query,policy_.k,[&](uint64_t word,size_t p,bool query_reverse){
            auto it=index_.find(word);if(it==index_.end())return;
            const auto&entry=words_[it->second];
            if(entry.loci>policy_.maximum_loci){++out.skipped_repetitive_seeds;return;}
            // hifiasm retains longest exact extensions before vote pooling.
            // A panel needs an additional distinction: donor routes sharing
            // this physical seed compete, but different loci do not. Otherwise
            // common flanks can give a wrong SNP route 249/250 union support,
            // or one extra matching base can falsely resolve a true paralog.
            // Occurrences are already grouped by physical seed locus. Keep
            // ALL maximum-length ties, without donor-frequency tie breaking.
            size_t group_begin=hits.size();uint32_t group=UINT32_MAX,longest=0;
            for(const auto&o:entry.occurrences) {
                if(o.locus!=group){group=o.locus;group_begin=hits.size();longest=0;}
                const auto&s=contexts_[o.context];
                const bool palindromic=!(policy_.k&1)&&query.substr(p,policy_.k)==reverse_complement(query.substr(p,policy_.k));
                for(unsigned orientation=0;orientation<(palindromic?2u:1u);++orientation) {
                    const bool reverse=palindromic?bool(orientation):query_reverse!=o.reverse;
                    const bool physical_reverse=reverse!=s.reverse;
                    const uint64_t start=physical_start(o);
                    const uint32_t qp=physical_reverse?length-uint32_t(p)-policy_.k:uint32_t(p);
                    const int64_t diagonal=int64_t(start)-int64_t(qp);
                    uint32_t begin=uint32_t(p),end=uint32_t(p)+policy_.k;
                    auto equal=[&](uint32_t q,int64_t pos){
                        if(pos<0||uint64_t(pos)>=s.sequence.size())return false;
                        const int a=base_code(query[q]),b=base_code(s.sequence[size_t(pos)]);
                        return a>=0&&b>=0&&a==(reverse?3-b:b);
                    };
                    const int64_t seed_base=int64_t(o.position)+(reverse?policy_.k-1:0);
                    auto target=[&](uint32_t q){return seed_base+(reverse?int64_t(p)-q:int64_t(q)-int64_t(p));};
                    const ExtensionKey extension{o.context,seed_base+(reverse?int64_t(p):-int64_t(p)),reverse};
                    auto [cached,inserted]=extensions.try_emplace(extension,MatchInterval{});
                    if(!inserted&&cached->second.begin<=begin&&cached->second.end>=end) {
                        begin=cached->second.begin;end=cached->second.end;
                    } else {
                        if(!reverse) {
                            const size_t left=std::min(size_t(begin),size_t(o.position));
                            begin-=uint32_t(exact_match::suffix(query.data()+p-left,s.sequence.data()+o.position-left,left));
                            const size_t right=std::min(size_t(length-end),s.sequence.size()-size_t(o.position)-policy_.k);
                            end+=uint32_t(exact_match::prefix(query.data()+end,s.sequence.data()+o.position+policy_.k,right));
                        } else {
                            while(begin&&equal(begin-1,target(begin-1)))--begin;
                            while(end<length&&equal(end,target(end)))++end;
                        }
                        cached->second={begin,end};
                    }
                    const uint32_t exact_length=end-begin;
                    if(exact_length<longest)continue;
                    if(exact_length>longest){hits.resize(group_begin);longest=exact_length;}
                    SelectedAnchor anchor;anchor.status=AnchorStatus::Unprojected;
                    if(o.projection!=UINT32_MAX){anchor=entry.projections.at(o.projection);anchor.word=word;}
                    hits.push_back({s.component,s.space,physical_reverse,diagonal,{begin,end},
                                    {word,start,uint32_t(p),o.context,physical_reverse},anchor});
                }
            }
        });
        if(hits.empty())return out;
        auto key=[](const Hit&h){return std::make_tuple(h.component,h.space,h.reverse,h.diagonal);};
        std::sort(hits.begin(),hits.end(),[&](const auto&a,const auto&b){
            return std::make_tuple(key(a),a.seed.query,a.seed.word,a.seed.context)<
                   std::make_tuple(key(b),b.seed.query,b.seed.word,b.seed.context);
        });
        const int64_t tolerance=int64_t(length/100)+1;
        std::vector<size_t> seed_order;
        for(size_t i=0;i<hits.size();) {
            size_t j=i+1;const auto&first=hits[i];
            while(j<hits.size()&&hits[j].component==first.component&&hits[j].space==first.space&&
                  hits[j].reverse==first.reverse&&hits[j].diagonal-first.diagonal<=tolerance)++j;
            PlacementCandidate c;c.component=first.component;c.space=first.space;c.reverse=first.reverse;
            c.diagonal_low=first.diagonal;c.diagonal_high=hits[j-1].diagonal;
            c.anchor.status=AnchorStatus::Unprojected;
            // Deduplicate seed identities without one tree allocation per
            // seed. The last sort restores first-occurrence order, including
            // duplicates across nearby diagonals, exactly as the former set.
            seed_order.clear();seed_order.reserve(j-i);
            c.intervals.reserve(j-i);
            for(size_t h=i;h<j;++h) {
                const auto&hit=hits[h];c.intervals.push_back(hit.interval);
                seed_order.push_back(h);
                if(hit.anchor.status==AnchorStatus::Known&&
                   (c.anchor.status!=AnchorStatus::Known||
                    std::tie(hit.anchor.start,hit.anchor.word)<std::tie(c.anchor.start,c.anchor.word)))
                    c.anchor=hit.anchor;
            }
            std::sort(seed_order.begin(),seed_order.end(),[&](size_t a,size_t b) {
                return std::tie(hits[a].seed.query,hits[a].seed.word,a)<
                       std::tie(hits[b].seed.query,hits[b].seed.word,b);
            });
            seed_order.erase(std::unique(seed_order.begin(),seed_order.end(),[&](size_t a,size_t b) {
                return std::tie(hits[a].seed.query,hits[a].seed.word)==
                       std::tie(hits[b].seed.query,hits[b].seed.word);
            }),seed_order.end());
            std::sort(seed_order.begin(),seed_order.end());
            c.seeds.reserve(seed_order.size());
            for(size_t h:seed_order)c.seeds.push_back(hits[h].seed);
            c.support=union_intervals(c.intervals);
            const auto&source=contexts_[first.seed.context];
            if(!source.runs.empty()) {
                const int64_t lo=c.diagonal_low,hi=c.diagonal_high+length;
                const uint64_t begin=uint64_t(std::max(int64_t(0),lo));
                const uint64_t end=uint64_t(std::min(int64_t(source.sequence.size()),hi));
                c.graph_context=placement_slices(source,begin,end);
                c.context_clip_left=lo<0?uint32_t(-lo):0;
                c.context_clip_right=hi>int64_t(source.sequence.size())?uint32_t(hi-int64_t(source.sequence.size())):0;
                c.context_jitter=c.diagonal_high-c.diagonal_low;
                if(c.reverse){reverse_placement_slices(c.graph_context);std::swap(c.context_clip_left,c.context_clip_right);}
            }
            if(c.seeds.size()>=policy_.minimum_seeds)out.candidates.push_back(std::move(c));
            i=j;
        }
        if(out.candidates.empty()){out.status=PlacementStatus::InsufficientSupport;out.anchor.status=AnchorStatus::InsufficientSupport;return out;}
        // Context aliases contribute a maximum, never a sum of donor votes.
        // Physical occurrence IDs keep repeated cyclic visits distinct.
        out.candidates=merge_placement_aliases(std::move(out.candidates));
        std::sort(out.candidates.begin(),out.candidates.end(),[](const auto&a,const auto&b){
            if(a.support!=b.support)return a.support>b.support;
            return std::tie(a.component,a.space,a.reverse,a.diagonal_low)<std::tie(b.component,b.space,b.reverse,b.diagonal_low);
        });
        out.best_support=out.candidates.front().support;
        if(out.candidates.size()>1)out.competing_support=out.candidates[1].support;
        if(uint64_t(out.best_support)*100<uint64_t(length)*policy_.minimum_support_percent) {
            out.status=PlacementStatus::InsufficientSupport;out.anchor.status=AnchorStatus::InsufficientSupport;return out;
        }
        bool split=false;
        for(size_t i=0;i<out.candidates.size();++i) {
            const auto&c=out.candidates[i];
            if(i&&uint64_t(overlap(c.intervals,out.candidates[0].intervals))*100<
                  uint64_t(std::min(c.support,out.best_support))*15)split=true;
            if(uint64_t(c.support)*100>=uint64_t(out.best_support)*policy_.competitor_percent)out.retained.push_back(i);
        }
        if(split) {out.status=PlacementStatus::Split;out.anchor.status=AnchorStatus::Split;}
        else {
            out.component=out.candidates[0].component;
            for(size_t i:out.retained)if(out.candidates[i].component!=out.component){out.component=UINT32_MAX;break;}
            if(out.retained.size()==1){out.status=PlacementStatus::Unique;out.anchor=out.candidates[0].anchor;}
            else {out.status=PlacementStatus::Ambiguous;out.anchor.status=AnchorStatus::Ambiguous;}
        }
        // Convert only presentation coordinates/strands back to the caller's
        // orientation. Anchor starts and all selection decisions are unchanged.
        if(reversed)for(auto&c:out.candidates) {
            c.reverse=!c.reverse;
            for(auto&x:c.intervals)x={length-x.end,length-x.begin};
            std::reverse(c.intervals.begin(),c.intervals.end());
            for(auto&seed:c.seeds){seed.query=length-seed.query-policy_.k;seed.reverse=!seed.reverse;}
        }
        return out;
    }
};
} // namespace hic
