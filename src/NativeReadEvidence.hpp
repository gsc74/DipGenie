#pragma once
#include "NativeGraphPlacement.hpp"
#include "ReadLibraries.hpp"

namespace hic {
// Bounded parallel gzip members, serialized in physical input order. Member
// boundaries depend only on byte offsets, not worker count or record size.
// The decompressed audit is identical to the original serial gzip stream.
class NativeObservationAudit {
    static constexpr size_t chunk_bytes=1u<<20,batch_chunks=32;
    std::ofstream file;
    RowExecutor pool;
    std::string current;
    std::vector<std::string> pending;
    bool closed=false,wrote=false;

    static std::string compress(const std::string&text) {
        z_stream z{};
        if(deflateInit2(&z,1,Z_DEFLATED,15+16,8,Z_DEFAULT_STRATEGY)!=Z_OK)
            throw std::runtime_error("Cannot initialize native audit compression");
        try {
            std::string bytes(deflateBound(&z,text.size()),'\0');
            z.next_in=reinterpret_cast<Bytef*>(const_cast<char*>(text.data()));
            z.avail_in=uInt(text.size());z.next_out=reinterpret_cast<Bytef*>(bytes.data());z.avail_out=uInt(bytes.size());
            if(deflate(&z,Z_FINISH)!=Z_STREAM_END)throw std::runtime_error("Native audit compression failed");
            bytes.resize(z.total_out);deflateEnd(&z);return bytes;
        }catch(...){deflateEnd(&z);throw;}
    }
    void flush() {
        if(pending.empty())return;
        std::vector<std::string> compressed(pending.size());
        pool.for_each(pending.size(),[&](size_t i){compressed[i]=compress(pending[i]);});
        for(const auto&bytes:compressed)file.write(bytes.data(),std::streamsize(bytes.size()));
        if(!file)throw std::runtime_error("Native observation audit write failed");
        wrote=true;pending.clear();
    }
public:
    explicit NativeObservationAudit(const std::filesystem::path&path,int threads=1)
        :file(path,std::ios::binary),pool(std::min(threads,int(batch_chunks))) {
        if(!file)throw std::runtime_error("Cannot create native observation audit");
        current.reserve(chunk_bytes);pending.reserve(batch_chunks);
    }
    NativeObservationAudit(const NativeObservationAudit&)=delete;
    NativeObservationAudit&operator=(const NativeObservationAudit&)=delete;
    void write(std::string_view s) {
        if(closed)throw std::logic_error("Native observation audit is closed");
        while(!s.empty()) {
            const size_t n=std::min(s.size(),chunk_bytes-current.size());
            current.append(s.data(),n);
            s.remove_prefix(n);
            if(current.size()==chunk_bytes) {
                pending.push_back(std::move(current));current.clear();current.reserve(chunk_bytes);
                if(pending.size()==batch_chunks)flush();
            }
        }
    }
    void close() {
        if(closed)return;
        if(!current.empty()||(!wrote&&pending.empty()))pending.push_back(std::move(current));
        flush();file.close();closed=true;
        if(!file)throw std::runtime_error("Native observation audit close failed");
    }
};
struct NativePairPolicy {
    static constexpr const char* version="unique-graph-pair-exact-anchor-v1";
    // This bounds selected-reference-anchor separation, NOT insert length.
    // Unknown geometry remains relation-only. No such gate applies to Hi-C.
    uint64_t maximum_sr_anchor_separation=5000;
    // Hi-C ligation junctions: a split or partially supported Hi-C end is
    // resolved to its 5'-most placement (resolve_hic_five_prime). Short-read
    // ends are never resolved this way.
    uint32_t hic_five_prime_min_support=60,hic_five_prime_competitor_percent=80;
};
// A Hi-C read that crosses a ligation junction yields several placements that
// cover disjoint parts of the read (Split) or one placement that covers only
// part of it (InsufficientSupport). Its 5' part is adjacent to the restriction
// site of the contact, the convention of Hi-C processing with bwa mem -5SPM and
// pairtools --walks-policy 5unique. The 5'-most candidate (smallest read
// start; ties by larger support, then candidate order) is accepted only if it
// has at least min_support exactly supported bases and no other candidate
// covering at least half of its read span has competitor_percent of its
// support. Markers are then taken from that candidate's intervals only. The
// read orientation is the sequenced one, so this is a protocol rule, not an
// anchor tie-break. Returns false and leaves the end unchanged otherwise.
inline bool resolve_hic_five_prime(EndPlacement&e,const NativePairPolicy&p) {
    if((e.status!=PlacementStatus::Split&&e.status!=PlacementStatus::InsufficientSupport)||e.candidates.empty())return false;
    auto first_base=[](const PlacementCandidate&c){uint32_t b=UINT32_MAX;for(const auto&x:c.intervals)b=std::min(b,x.begin);return b;};
    size_t best=0;
    for(size_t i=1;i<e.candidates.size();++i) {
        const auto bi=first_base(e.candidates[i]),bb=first_base(e.candidates[best]);
        if(bi<bb||(bi==bb&&e.candidates[i].support>e.candidates[best].support))best=i;
    }
    const auto&c=e.candidates[best];
    if(c.intervals.empty()||c.support<p.hic_five_prime_min_support)return false;
    uint32_t lo=UINT32_MAX,hi=0;
    for(const auto&x:c.intervals){lo=std::min(lo,x.begin);hi=std::max(hi,x.end);}
    std::vector<unsigned char> covered;
    for(size_t i=0;i<e.candidates.size();++i) {
        if(i==best)continue;
        const auto&o=e.candidates[i];
        covered.assign(hi-lo,0);
        for(const auto&x:o.intervals)for(uint32_t q=std::max(x.begin,lo);q<std::min(x.end,hi);++q)covered[q-lo]=1;
        uint64_t overlap=0;for(auto v:covered)overlap+=v;
        if(overlap*2>=uint64_t(hi-lo)&&uint64_t(o.support)*100>=uint64_t(c.support)*p.hic_five_prime_competitor_percent)return false;
    }
    e.status=PlacementStatus::Unique;e.component=c.component;e.anchor=c.anchor;
    e.retained={best};e.best_support=c.support;
    return true;
}
inline const char* placement_status_name(PlacementStatus s) {
    switch(s) {
        case PlacementStatus::NoSeeds:return "no_seeds";
        case PlacementStatus::InsufficientSupport:return "insufficient";
        case PlacementStatus::Unique:return "unique";
        case PlacementStatus::Ambiguous:return "ambiguous";
        case PlacementStatus::Split:return "split";
    }
    throw std::logic_error("Invalid native placement status");
}
inline const char* anchor_status_name(AnchorStatus s) {
    switch(s) {
        case AnchorStatus::Known:return "known";
        case AnchorStatus::NoSeeds:return "no_seeds";
        case AnchorStatus::InsufficientSupport:return "insufficient";
        case AnchorStatus::Ambiguous:return "ambiguous";
        case AnchorStatus::Split:return "split";
        case AnchorStatus::Unprojected:return "unprojected";
    }
    throw std::logic_error("Invalid native anchor status");
}
inline const char* native_observation_header() {
    return "library\tsource\tmolecule\tadmission\tdistance_status\tanchor_separation\tqc\tend1\tend2\tmarkers1\tmarkers2\n";
}
inline PairAdmissionResult native_pair_evidence(const SequenceIndex&index,const NativeLocalizer&localizer,
        const FastqRead&a,const FastqRead&b,const std::string&id,MoleculeKind source,uint32_t library,
        const NativePairPolicy&policy={},bool audit=true) {
    if(index.k!=localizer.policy().k||index.w!=localizer.policy().w)
        throw std::invalid_argument("Native localization and target selection disagree on k/w");
    EndPlacement ends[2]={localizer.locate(a.sequence),localizer.locate(b.sequence)};
    bool five_prime[2]={false,false};
    if(source==MoleculeKind::HiC)
        for(int i=0;i<2;++i)five_prime[i]=resolve_hic_five_prime(ends[i],policy);
    const auto distance=selected_anchor_distance(ends[0],ends[1]);
    PairAdmissionResult result;auto&m=result.marker;
    m.id=id;m.kind=source;m.library=library;
    m.geometry=distance.status==AnchorDistanceStatus::Known?MoleculeGeometry::Known:
               distance.status==AnchorDistanceStatus::Trans?MoleculeGeometry::Trans:MoleculeGeometry::Unknown;
    if(distance.status==AnchorDistanceStatus::Known) {
        if(distance.separation>INT64_MAX)throw std::overflow_error("Native anchor separation");
        m.reference_distance=int64_t(distance.separation);
    }
    std::string admission="accepted";
    // Conservative locus admission: a component alone does not establish a
    // physical placement. Unknown reference geometry is still permitted when
    // graph placement itself is unique (e.g. an alternate sequence).
    if(ends[0].status!=PlacementStatus::Unique||ends[1].status!=PlacementStatus::Unique)admission="unresolved_end";
    else if(ends[0].component!=ends[1].component)admission="cross_component";
    else if(source==MoleculeKind::ShortRead&&distance.status==AnchorDistanceStatus::Trans)admission="sr_trans";
    else if(source==MoleculeKind::ShortRead&&distance.status==AnchorDistanceStatus::Known&&
            distance.separation>policy.maximum_sr_anchor_separation)admission="sr_anchor_separation";
    else {
        m.localization_qc=ends[0].best_support==a.sequence.size()&&ends[1].best_support==b.sequence.size()?1:2;
        m.ends[0]=native_interval_markers(index,a.sequence,ends[0].candidates.at(ends[0].retained.at(0)).intervals);
        m.ends[1]=native_interval_markers(index,b.sequence,ends[1].candidates.at(ends[1].retained.at(0)).intervals);
        if(m.ends[0].empty()||m.ends[1].empty())admission="missing_markers";
    }
    if(audit) {
        std::ostringstream row;
        row<<library<<'\t'<<(source==MoleculeKind::ShortRead?"SR":"HIC")<<'\t'<<id<<'\t'<<admission<<'\t';
        row<<(distance.status==AnchorDistanceStatus::Known?"known":distance.status==AnchorDistanceStatus::Trans?"trans":"unknown")<<'\t';
        if(distance.status==AnchorDistanceStatus::Known)row<<distance.separation;else row<<'.';
        row<<'\t'<<unsigned(m.localization_qc);
        for(int i=0;i<2;++i) {
            const auto&e=ends[i];
            row<<'\t'<<(five_prime[i]?"five_prime":placement_status_name(e.status))<<",";
            if(e.component_unique())row<<e.component;else row<<'.';
            row<<','<<e.best_support<<','<<e.competing_support<<','<<e.skipped_repetitive_seeds
               <<','<<anchor_status_name(e.anchor.status)<<',';
            if(e.anchor.status==AnchorStatus::Known)row<<e.anchor.contig<<':'<<e.anchor.start<<':'<<e.anchor.word;else row<<'.';
            // Candidate records: component/space/strand/support/[query spans]/
            // [exact seed word:physical position:query position:context].
            for(const auto&c:e.candidates) {
                row<<';'<<c.component<<'/'<<c.space<<'/'<<(c.reverse?'-':'+')<<'/'<<c.support<<'/';
                for(const auto&r:c.intervals)row<<r.begin<<'-'<<r.end<<',';
                row<<'/';for(const auto&s:c.seeds)row<<s.word<<':'<<s.position<<':'<<s.query<<':'<<s.context<<',';
            }
        }
        for(const auto&markers:m.ends){row<<'\t';for(size_t x:markers)row<<x<<',';}
        row<<'\n';result.audit=row.str();
    }
    if(admission!="accepted"){m.ends[0].clear();m.ends[1].clear();}
    return result;
}
using LibraryProgressSink=std::function<void(const ReadLibrary&,uint64_t,uint64_t,bool)>;
inline LibraryEvidence read_native_library_evidence(SequenceIndex&index,std::vector<ReadLibrary> libraries,
        const NativeLocalizer&localizer,int threads=1,const NativePairPolicy&policy={},const PairAuditSink&audit={},
        const LibraryProgressSink&progress={}) {
    validate_libraries(libraries);
    for(const auto&l:libraries) {
        if(l.name.find_first_of("\t\r\n")!=std::string::npos)throw std::invalid_argument("Native library name contains a record separator");
    }
    if(audit) {
        const auto&p=localizer.policy();std::ostringstream metadata;
        metadata<<"#native_observation_schema\t1\n#localizer\t"<<NativeLocalizer::version
            <<"\n#graph_coordinates\t"<<native_graph_context_version<<"\n#pair_policy\t"<<NativePairPolicy::version
            <<"\n#k\t"<<p.k<<"\n#w\t"<<p.w<<"\n#maximum_loci\t"<<p.maximum_loci
            <<"\n#minimum_seeds\t"<<p.minimum_seeds<<"\n#minimum_support_percent\t"<<p.minimum_support_percent
            <<"\n#competitor_percent\t"<<p.competitor_percent<<"\n#maximum_read_length\t"<<p.maximum_read_length
            <<"\n#maximum_sr_anchor_separation\t"<<policy.maximum_sr_anchor_separation
            <<"\n#hic_split\tfive-prime"
            <<"\n#hic_five_prime_min_support\t"<<policy.hic_five_prime_min_support
            <<"\n#hic_five_prime_competitor_percent\t"<<policy.hic_five_prime_competitor_percent
            <<"\n#accepted_observation_budget_units\t1024\n#qc\t1=full_query_support,2=partial_query_support\n";
        for(size_t i=0;i<libraries.size();++i)metadata<<"#library\t"<<i<<'\t'<<libraries[i].name<<'\t'<<libraries[i].kind<<'\n';
        for(size_t i=0;i<localizer.contexts().size();++i) {
            const auto&s=localizer.contexts()[i];
            if(s.name.find_first_of("\t\r\n")!=std::string::npos)throw std::invalid_argument("Native context name contains a record separator");
            metadata<<"#context\t"<<i<<'\t'<<s.name<<'\t'<<s.component<<'\t'<<s.space<<'\t'<<s.sequence.size()<<'\t';
            if(s.reference_contig!=UINT32_MAX)metadata<<s.reference_contig;else metadata<<'.';
            metadata<<'\n';
        }
        audit(metadata.str());audit(native_observation_header());
    }
    return read_library_evidence_with(index,libraries,[&](const ReadLibrary&l,uint32_t library) {
        const auto kind=l.kind=="SR"?MoleculeKind::ShortRead:MoleculeKind::HiC;
        PairAdmissionFunction admit=[&](const FastqRead&a,const FastqRead&b,const std::string&id) {
            return native_pair_evidence(index,localizer,a,b,id,kind,library,policy,bool(audit));
        };
        PairProgressSink batch_progress;
        if(progress) {
            progress(l,0,0,false);
            batch_progress=[&](uint64_t pairs,uint64_t linked){progress(l,pairs,linked,false);};
        }
        auto raw=read_evidence(index,l.r1,l.r2,{},l.kind=="SR",threads,admit,audit,batch_progress);
        // A malformed stream never reports completion. Empty libraries are
        // rejected by read_library_evidence_with, not announced as successes.
        if(progress&&raw.pairs)progress(l,raw.pairs,raw.markers.size(),true);
        return raw;
    });
}
} // namespace hic
