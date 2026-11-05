#pragma once
#include "NativeGraph.hpp"
#include "PresenceDP.hpp"
#include "ShardedWords.hpp"
#include "FlatHash.hpp"
#include <deque>
#include <fstream>
#include <set>
#include <string_view>
#include <iterator>
#include <functional>

namespace hic {
inline int base_code(char c) {
    switch(std::toupper(static_cast<unsigned char>(c))){case'A':return 0;case'C':return 1;case'G':return 2;case'T':return 3;default:return -1;}
}
template<class F> void words(std::string_view s,unsigned k,F visit) {
    if(k<1||k>31)throw std::invalid_argument("2-bit word length");
    const uint64_t mask=(uint64_t(1)<<(2*k))-1;uint64_t a=0,b=0;size_t valid=0;
    for(size_t p=0;p<s.size();++p) {
        int c=base_code(s[p]);if(c<0){a=b=0;valid=0;continue;}
        a=((a<<2)|uint64_t(c))&mask;b=(b>>2)|(uint64_t(3-c)<<(2*(k-1)));++valid;
        if(valid>=k)visit(std::min(a,b),p+1-k,b<a);
    }
}
// Exact lexicographic canonical-word selection used by v1, including its
// rightmost tie rule. 2-bit identities replace a collision-prone hash key.
template<class F> void minimizers(std::string_view s,unsigned k,unsigned w,F visit) {
    if(!w)throw std::invalid_argument("zero minimizer window");
    std::deque<std::pair<uint64_t,size_t>> q;size_t previous=SIZE_MAX,run=0;uint64_t last=UINT64_MAX;
    words(s,k,[&](uint64_t word,size_t p,bool reverse){
        (void)reverse;
        if(previous==SIZE_MAX||p!=previous+1){q.clear();run=0;last=UINT64_MAX;}
        previous=p;++run;
        while(!q.empty()&&q.back().first>=word)q.pop_back();
        q.emplace_back(word,p);
        while(!q.empty()&&q.front().second+w<=p)q.pop_front();
        if(run>=w&&q.front().first!=last){last=q.front().first;visit(last,q.front().second);}
    });
}
struct FastqRead {std::string id,sequence,quality;int header_mate=0;};
class FastqInput {
    TextInput input;
public:
    explicit FastqInput(const std::string&path):input(path){}
    bool next(FastqRead&r) {
        std::string header,plus;
        if(!input.getline(header))return false;
        if(!input.getline(r.sequence)||!input.getline(plus)||!input.getline(r.quality)||header.empty()||header[0]!='@'||plus.empty()||plus[0]!='+'||r.sequence.size()!=r.quality.size())
            throw std::runtime_error("Malformed four-line FASTQ");
        const size_t split=header.find_first_of(" \t");
        r.id=header.substr(1,split-1);r.header_mate=0;
        if(r.id.empty())throw std::runtime_error("Empty FASTQ read name");
        if(split!=std::string::npos) {
            const auto p=header.find_first_not_of(" \t",split);
            if(p!=std::string::npos&&p+1<header.size()&&header[p+1]==':'&&(header[p]=='1'||header[p]=='2'))
                r.header_mate=header[p]-'0';
        }
        return true;
    }
};
inline std::string molecule_id(const std::string&id,int mate) {
    const auto suffix="/"+std::to_string(mate);
    if(id.size()<2||id.substr(id.size()-2)!=suffix)throw std::runtime_error("Explicit /1 /2 mate labels required");
    return id.substr(0,id.size()-2);
}
inline std::string paired_molecule_id(const FastqRead&a,const FastqRead&b) {
    auto id=[](const FastqRead&r,int expected) {
        if(r.header_mate&&r.header_mate!=expected)throw std::runtime_error("FASTQ header mate contradicts input file");
        auto name=r.id;
        if(name.size()>=2&&name[name.size()-2]=='/') {
            if(name.back()!=char('0'+expected))throw std::runtime_error("FASTQ suffix mate contradicts input file");
            name.resize(name.size()-2);
        }
        if(name.empty())throw std::runtime_error("Empty physical-molecule name");
        return name;
    };
    auto one=id(a,1),two=id(b,2);
    if(one!=two)throw std::runtime_error("FASTQ mates differ");
    // Unmarked names are accepted only as matching records in explicit R1/R2
    // files. No singleton rescue, positional name guessing or pair synthesis.
    return one;
}

struct WordTarget {
    uint64_t word=0,count=0,location_start=0,location_end=0;
    // All exact occurrences in read ends, without per-end deduplication.
    // Copy-number evidence for tandem strings; `count` keeps one vote/end.
    uint64_t occurrence_count=0;
    DonorMask carriers;
    bool location_set=false,unique_locus=true;
    std::vector<std::pair<int,DonorMask>> occurrences;
    // Positive copy number -> number of panel haplotypes with that count.
    // Zero-copy haplotypes are H-popcount(carriers). Kept during indexing,
    // rather than treating a multi-copy word as a single-copy allele.
    std::vector<std::pair<uint32_t,uint32_t>> copy_histogram;
    // (path, copies) for every path carrying two or more copies, in path
    // order. Carriers absent from this list have exactly one copy.
    std::vector<std::pair<uint32_t,uint32_t>> donor_copies;
    uint32_t copies(size_t h)const {
        if(!carriers.test(h))return 0;
        auto it=std::lower_bound(donor_copies.begin(),donor_copies.end(),std::make_pair(uint32_t(h),uint32_t(0)));
        return it!=donor_copies.end()&&it->first==h?it->second:1;
    }
};
struct SequenceIndex {
    unsigned k=31,w=25;
    std::vector<WordTarget> targets;
    FlatMap<uint64_t,size_t> ids;
    static SequenceIndex build(const NativeGraph&g,const CommonSegments&parts,
                               unsigned k=31,unsigned w=25,int threads=1) {
        (void)checked_pair_domain(g.paths.size());
        if(k<1||k>31||!w||threads<1||parts.sequences.size()!=g.paths.size()||
           parts.cuts.size()!=g.paths.size())
            throw std::invalid_argument("Sequence index dimensions/options");
        SequenceIndex out;out.k=k;out.w=w;
        RowExecutor pool(threads);
        constexpr size_t discovery_chunk=1u<<20;
        std::vector<uint64_t> order;
        // Parallel window discovery with k+w-2 bases of right context.
        // Every complete minimizer window is evaluated exactly once. Repeated
        // emitted identities across chunk boundaries disappear in set-union.
        // Process one haplotype at a time: scratch does not multiply by H.
        for(const auto&s:parts.sequences) {
            const size_t context=size_t(k)+w-2;
            const size_t windows=s.size()>context?s.size()-context:0;
            const size_t chunks=windows/discovery_chunk+(windows%discovery_chunk!=0);
            std::vector<std::vector<uint64_t>> selected(chunks);
            pool.for_each(chunks,[&](size_t i) {
                const size_t begin=i*discovery_chunk;
                const size_t n=std::min(discovery_chunk,windows-begin);
                auto&v=selected[i];
                minimizers(std::string_view(s).substr(begin,n+context),k,w,
                           [&](uint64_t word,size_t){v.push_back(word);});
                sort_words(v);v.erase(std::unique(v.begin(),v.end()),v.end());
            });
            while(selected.size()>1) {
                std::vector<std::vector<uint64_t>> merged((selected.size()+1)/2);
                pool.for_each(merged.size(),[&](size_t i) {
                    auto&a=selected[2*i];
                    if(2*i+1==selected.size()){merged[i]=std::move(a);return;}
                    auto&b=selected[2*i+1];auto&dest=merged[i];
                    dest.reserve(a.size()+b.size());
                    std::set_union(a.begin(),a.end(),b.begin(),b.end(),std::back_inserter(dest));
                    std::vector<uint64_t>().swap(a);std::vector<uint64_t>().swap(b);
                });
                selected=std::move(merged);
            }
            if(selected.empty())continue;
            if(order.empty()){order=std::move(selected.front());continue;}
            std::vector<uint64_t> merged;merged.reserve(order.size()+selected.front().size());
            std::set_union(order.begin(),order.end(),selected.front().begin(),selected.front().end(),
                           std::back_inserter(merged));
            order=std::move(merged);
        }
        out.ids.reserve(order.size());out.targets.resize(order.size());
        for(size_t x=0;x<order.size();++x){out.ids.emplace(order[x],x);out.targets[x].word=order[x];}
        std::vector<uint64_t>().swap(order);
        if(out.targets.empty())return out;
        // The lookup is immutable in all worker tasks (including const find).
        const auto&ids=out.ids;
        const size_t shards=std::min(out.targets.size(),size_t(threads)*4);
        const size_t per_shard=out.targets.size()/shards+(out.targets.size()%shards!=0);
        struct Occurrence {size_t id;uint64_t start,end;int block;};
        struct Chunk {std::vector<Occurrence> data;std::vector<size_t> offsets;};
        constexpr size_t scratch_bytes=256u<<20;
        const size_t chunk_bases=std::max(size_t(1024),std::min(size_t(1u<<16),
            scratch_bytes/(2*sizeof(Occurrence)*size_t(threads)*2)));
        const size_t batch_limit=std::max(size_t(1),std::min(size_t(threads)*2,
            scratch_bytes/(2*sizeof(Occurrence)*chunk_bases)));
        std::vector<uint32_t> count(out.targets.size());
        for(size_t h=0;h<g.paths.size();++h) {
            const auto&walk=g.paths[h].walk;const auto&sequence=parts.sequences[h];
            const auto&cuts=parts.cuts[h];
            if(cuts.size()<2||cuts.front()!=0||cuts.back()!=sequence.size()||
               !std::is_sorted(cuts.begin(),cuts.end()))
                throw std::invalid_argument("Sequence index cuts");
            // One offset per graph occurrence replaces one uint64 per base.
            // Two monotone cursors recover exactly the former physical start
            // and end coordinates, including reverse-oriented small vertices.
            std::vector<size_t> offsets(walk.size()+1);
            for(size_t i=0;i<walk.size();++i) {
                const size_t length=g.labels.at(walk[i]/2).size();
                if(length>UINT32_MAX||length>SIZE_MAX-offsets[i])
                    throw std::overflow_error("Vertex coordinates");
                offsets[i+1]=offsets[i]+length;
            }
            if(offsets.back()!=sequence.size())throw std::logic_error("Panel spelling length mismatch");
            const size_t starts=sequence.size()>=k?sequence.size()-k+1:0;
            const size_t total_chunks=starts/chunk_bases+(starts%chunk_bases!=0);
            for(size_t first=0;first<total_chunks;first+=batch_limit) {
                std::vector<Chunk> batch(std::min(batch_limit,total_chunks-first));
                pool.for_each(batch.size(),[&](size_t i) {
                    const size_t begin=(first+i)*chunk_bases;
                    const size_t n=std::min(chunk_bases,starts-begin);
                    size_t left=size_t(std::upper_bound(offsets.begin(),offsets.end(),begin)-offsets.begin()-1);
                    size_t right=left;
                    size_t block=size_t(std::lower_bound(cuts.begin()+1,cuts.end(),begin+k)-cuts.begin()-1);
                    std::vector<Occurrence> found;found.reserve(n);
                    auto coordinate=[&](size_t p,size_t&cursor) {
                        while(cursor+1<walk.size()&&p>=offsets[cursor+1])++cursor;
                        return (uint64_t(walk[cursor])<<32)|uint32_t(p-offsets[cursor]);
                    };
                    auto flipped=[&](uint64_t p) {
                        const uint32_t v=uint32_t(p>>32),q=uint32_t(p);
                        return (uint64_t(v^1)<<32)|(g.labels[v/2].size()-1-q);
                    };
                    words(std::string_view(sequence).substr(begin,n+k-1),k,
                          [&](uint64_t word,size_t relative,bool reverse) {
                        auto it=ids.find(word);if(it==ids.end())return;
                        const size_t p=begin+relative;
                        while(p+k>cuts[block+1])++block;
                        const uint64_t a=coordinate(p,left),b=coordinate(p+k-1,right);
                        found.push_back({it->second,reverse?flipped(b):a,reverse?flipped(a):b,int(block)});
                    });
                    // Stable owner partition. Owners update disjoint contiguous
                    // target ranges, so no per-hit mutex or atomic is required.
                    auto&chunk=batch[i];chunk.offsets.assign(shards+1,0);
                    for(const auto&r:found)++chunk.offsets[r.id/per_shard+1];
                    for(size_t s=0;s<shards;++s)chunk.offsets[s+1]+=chunk.offsets[s];
                    auto cursor=chunk.offsets;chunk.data.resize(found.size());
                    for(const auto&r:found)chunk.data[cursor[r.id/per_shard]++]=r;
                });
                pool.for_each(shards,[&](size_t owner) {
                    // Chunk order and record order preserve the serial first
                    // locus, carrier bits, copy histogram and occurrence list.
                    for(const auto&chunk:batch)for(size_t i=chunk.offsets[owner];i<chunk.offsets[owner+1];++i) {
                        const auto&r=chunk.data[i];auto&t=out.targets[r.id];
                        t.carriers.set(h);
                        auto it=std::lower_bound(t.occurrences.begin(),t.occurrences.end(),r.block,
                                                [](const auto&a,int b){return a.first<b;});
                        if(it==t.occurrences.end()||it->first!=r.block)
                            it=t.occurrences.insert(it,{r.block,DonorMask{}});
                        it->second.set(h);
                        auto&copies=count[r.id];
                        if(copies==UINT32_MAX)throw std::overflow_error("Panel string multiplicity");
                        if(copies++)t.unique_locus=false;
                        if(!t.location_set){t.location_set=true;t.location_start=r.start;t.location_end=r.end;}
                        else if(t.location_start!=r.start||t.location_end!=r.end)t.unique_locus=false;
                    }
                });
            }
            pool.for_each(shards,[&](size_t owner) {
                const size_t last=std::min(out.targets.size(),(owner+1)*per_shard);
                for(size_t x=owner*per_shard;x<last;++x)if(count[x]) {
                    auto&hist=out.targets[x].copy_histogram;
                    auto it=std::lower_bound(hist.begin(),hist.end(),count[x],[](auto a,uint32_t n){return a.first<n;});
                    if(it==hist.end()||it->first!=count[x])hist.insert(it,{count[x],1});
                    else ++it->second;
                    if(count[x]>1)out.targets[x].donor_copies.emplace_back(uint32_t(h),count[x]);
                    count[x]=0;
                }
            });
        }
        return out;
    }
    std::vector<std::pair<size_t,size_t>> hits(const std::string&s)const {
        std::vector<std::pair<size_t,size_t>> result;
        minimizers(s,k,w,[&](uint64_t word,size_t position){auto it=ids.find(word);if(it!=ids.end())result.push_back({it->second,position});});
        return result;
    }
    void count_end(const std::vector<std::pair<size_t,size_t>>&hits) {
        std::vector<size_t> seen;for(auto[id,p]:hits){(void)p;seen.push_back(id);}
        std::sort(seen.begin(),seen.end());seen.erase(std::unique(seen.begin(),seen.end()),seen.end());
        for(size_t id:seen){if(targets[id].count==UINT64_MAX)throw std::overflow_error("read count");++targets[id].count;}
    }
};

struct ContactEnd {int mate=0;char strand='+';std::string cigar;int mapq=30;};
// Reference-coordinate distance is metadata, not a guessed graph distance.
// Negative sentinels distinguish missing coordinates from interchromosomal
// contacts. Neither may be treated as a zero-distance same-chromosome link.
constexpr int64_t unknown_contact_distance=-1;
constexpr int64_t interchromosomal_contact_distance=-2;
struct InputContact {
    ContactEnd ends[2];char type[3]="UU";
    int64_t reference_distance=unknown_contact_distance;
};
inline void check_contact_mates(const InputContact&contact) {
    const int a=contact.ends[0].mate,b=contact.ends[1].mate;
    if(!((a==1&&b==2)||(a==2&&b==1)))
        throw std::runtime_error("Hi-C phase links require one end from each FASTQ mate");
}
inline std::pair<size_t,size_t> aligned_query_range(const ContactEnd&e,size_t length) {
    size_t q=0,reference=0,begin=SIZE_MAX,end=0,number=0;
    bool digits=false,body=false,trailing_clip=false,matched=false;
    for(unsigned char c:e.cigar) {
        if(c>='0'&&c<='9') {
            const size_t digit=size_t(c-'0');
            if(number>(SIZE_MAX-digit)/10)throw std::runtime_error("Contact CIGAR length overflow");
            number=number*10+digit;digits=true;continue;
        }
        if(!digits||!number)throw std::runtime_error("Invalid contact CIGAR");
        if(c=='H')throw std::runtime_error("Hard-clipped Hi-C contact; -Y required");
        const bool alignment=c=='M'||c=='='||c=='X';
        const bool query=alignment||c=='I'||c=='S';
        const bool ref=alignment||c=='D'||c=='N';
        if(!query&&!ref&&c!='P')throw std::runtime_error("Unknown CIGAR operation");
        if(c=='S') {if(body)trailing_clip=true;}
        else {
            if(trailing_clip)throw std::runtime_error("Internal soft clipping in contact CIGAR");
            body=true;
        }
        if(ref) {
            if(number>SIZE_MAX-reference)throw std::runtime_error("Contact reference span overflow");
            reference+=number;
        }
        if(query) {
            if(number>length-q)throw std::runtime_error("Contact CIGAR exceeds FASTQ length");
            if(c!='S')begin=std::min(begin,q);
            q+=number;
            if(c!='S')end=q;
        }
        matched=matched||alignment;
        number=0;digits=false;
    }
    if(digits||q!=length||begin==SIZE_MAX||!matched)throw std::runtime_error("Contact CIGAR/read length mismatch");
    if(e.strand=='-')return{length-end,length-begin};
    if(e.strand!='+')throw std::runtime_error("Invalid contact strand");
    return{begin,end};
}
enum class MoleculeKind : unsigned char { ShortRead,HiC };
enum class MoleculeGeometry : unsigned char {Legacy,Known,Unknown,Trans};
struct MoleculeMarkers {
    std::string id;std::vector<size_t> ends[2];
    int64_t reference_distance=unknown_contact_distance;
    uint32_t library=UINT32_MAX;
    MoleculeKind kind=MoleculeKind::ShortRead;
    MoleculeGeometry geometry=MoleculeGeometry::Legacy;
    // Native QC: 1=both query ends fully supported, 2=partial exact support.
    // This is a deterministic support class, not calibrated MAPQ.
    unsigned char localization_qc=0;
};
struct PairAdmissionResult {MoleculeMarkers marker;std::string audit;};
using PairAdmissionFunction=std::function<PairAdmissionResult(const FastqRead&,const FastqRead&,const std::string&)>;
using PairAuditSink=std::function<void(std::string_view)>;
// Invoked by the coordinator after each fully consumed batch. Counts are
// per library; no read prepass or mutable access to inference data is exposed.
using PairProgressSink=std::function<void(uint64_t,uint64_t)>;
struct RawEvidence {
    uint64_t pairs=0,contact_pairs=0;
    uint64_t target_end_observations=0,minimizer_target_end_observations=0;
    std::vector<MoleculeMarkers> markers;
    std::map<uint64_t,uint64_t> count_histogram;
};
inline RawEvidence read_evidence(SequenceIndex&index,const std::string&r1,const std::string&r2,
                                std::unordered_map<std::string,InputContact> contacts,
                                bool sr_pairs=false,int threads=1,const PairAdmissionFunction&admit={},
                                const PairAuditSink&audit={},const PairProgressSink&progress={}) {
    if(threads<1)throw std::invalid_argument("Positive evidence worker count required");
    if(sr_pairs&&!contacts.empty())throw std::invalid_argument("SR pairs do not take Hi-C contacts");
    if(admit&&!contacts.empty())throw std::invalid_argument("Native evidence does not take alignment-derived contacts");
    for(const auto&entry:contacts)check_contact_mates(entry.second);
    FastqInput inputs[2]={FastqInput(r1),FastqInput(r2)};
    RawEvidence result;RowExecutor pool(threads);
    std::unordered_set<std::string> native_names;
    const auto&ids=index.ids;
    const size_t shards=std::min(size_t(256),size_t(threads)*4);
    std::vector<FlatMap<uint64_t,uint64_t>> all_read_counts(shards);
    auto owner=[&](uint64_t x) {
        x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;x^=x>>31;
        return size_t(x%shards);
    };
    // Keep the original word-based owner, but carry resolved target IDs so
    // neither count reduction needs a second word->target hash lookup.
    auto target_owner=[&](uint64_t id){return owner(index.targets[size_t(id)].word);};
    constexpr size_t batch_pairs=4096,task_pairs=64;
    struct Pair {std::string id;InputContact contact;bool phase=false;};
    struct Task {
        ShardedWords words,targets,occurrences;
        std::vector<MoleculeMarkers> markers;
        std::vector<std::string> audit;
        uint64_t target_votes=0,minimizer_votes=0;
    };
    std::vector<Task> tasks((batch_pairs+task_pairs-1)/task_pairs);
    for(;;) {
        std::vector<FastqRead> reads[2];
        pool.for_each(2,[&](size_t side) {
            reads[side].reserve(batch_pairs);
            FastqRead next;
            while(reads[side].size()<batch_pairs&&inputs[side].next(next))
                reads[side].push_back(std::move(next));
        });
        if(reads[0].size()!=reads[1].size())throw std::runtime_error("Different mate counts");
        const size_t n=reads[0].size();if(!n)break;
        if(n>UINT64_MAX-result.pairs)throw std::overflow_error("Read pair count");
        result.pairs+=n;
        std::vector<Pair> pairs(n);
        pool.for_each(n,[&](size_t i){pairs[i].id=paired_molecule_id(reads[0][i],reads[1][i]);});
        // Consume contacts in input order, preserving duplicate-ID behavior.
        // Workers never mutate or search the name table.
        for(auto&p:pairs) {
            if(admit&&!native_names.insert(p.id).second)
                throw std::invalid_argument("Duplicate physical-molecule name in native read library: "+p.id);
            auto found=contacts.find(p.id);
            p.phase=bool(admit)||sr_pairs||found!=contacts.end();
            if(!p.phase)continue;
            ++result.contact_pairs;
            if(!sr_pairs&&!admit){p.contact=std::move(found->second);contacts.erase(found);}
        }
        tasks.resize((n+task_pairs-1)/task_pairs);
        pool.for_each(tasks.size(),[&](size_t task) {
            auto&dest=tasks[task];
            dest.words.clear();dest.targets.clear();dest.occurrences.clear();
            dest.markers.clear();dest.audit.clear();dest.target_votes=dest.minimizer_votes=0;
            const size_t finish=std::min(n,(task+1)*task_pairs);
            std::vector<uint64_t> distinct,matched;
            for(size_t i=task*task_pairs;i<finish;++i) {
                std::vector<std::pair<size_t,size_t>> hits[2];
                for(size_t end=0;end<2;++end) {
                    distinct.clear();
                    minimizers(reads[end][i].sequence,index.k,index.w,[&](uint64_t word,size_t p) {
                        distinct.push_back(word);
                        // Native linkage uses positional localizer observations,
                        // not this legacy marker list. Raw counts are unchanged.
                        if(!admit){auto it=ids.find(word);if(it!=ids.end())hits[end].push_back({it->second,p});}
                    });
                    std::sort(distinct.begin(),distinct.end());
                    distinct.erase(std::unique(distinct.begin(),distinct.end()),distinct.end());
                    for(uint64_t word:distinct)dest.words.push_back(word);
                    for(uint64_t word:distinct)dest.minimizer_votes+=ids.find(word)!=ids.end();
                    // The panel's minimizers define the fixed target universe.
                    // Count every exact target observed in an end, regardless
                    // of whether that end's local window selects it again.
                    // One end contributes at most once to a target, including
                    // repeated copies. Phase markers retain the old minimizer
                    // selection and alignment-range restrictions above/below.
                    matched.clear();
                    words(reads[end][i].sequence,index.k,[&](uint64_t word,size_t,bool) {
                        auto found=ids.find(word);
                        if(found!=ids.end())matched.push_back(found->second);
                    });
                    std::sort(matched.begin(),matched.end());
                    for(uint64_t word:matched)dest.occurrences.push_back(word);
                    matched.erase(std::unique(matched.begin(),matched.end()),matched.end());
                    dest.target_votes+=matched.size();
                    for(uint64_t word:matched)dest.targets.push_back(word);
                }
                if(!pairs[i].phase)continue;
                if(admit) {
                    // Dosage has already consumed the ORIGINAL full ends.
                    // Localization changes linkage admission only.
                    auto localized=admit(reads[0][i],reads[1][i],pairs[i].id);
                    if(audit)dest.audit.push_back(std::move(localized.audit));
                    auto&m=localized.marker;
                    if(!m.ends[0].empty()&&!m.ends[1].empty())dest.markers.push_back(std::move(m));
                    continue;
                }
                MoleculeMarkers molecule;molecule.id=pairs[i].id;
                molecule.kind=sr_pairs?MoleculeKind::ShortRead:MoleculeKind::HiC;
                if(!sr_pairs)molecule.reference_distance=pairs[i].contact.reference_distance;
                for(size_t end=0;end<2;++end) {
                    const size_t side=sr_pairs?end:size_t(pairs[i].contact.ends[end].mate-1);
                    const auto range=sr_pairs?std::make_pair(size_t(0),reads[side][i].sequence.size()):
                        aligned_query_range(pairs[i].contact.ends[end],reads[side][i].sequence.size());
                    for(auto[x,p]:hits[side])
                        if(index.targets[x].unique_locus&&p>=range.first&&p+index.k<=range.second)
                            molecule.ends[end].push_back(x);
                    auto&v=molecule.ends[end];std::sort(v.begin(),v.end());
                    v.erase(std::unique(v.begin(),v.end()),v.end());
                }
                if(!molecule.ends[0].empty()&&!molecule.ends[1].empty())
                    dest.markers.push_back(std::move(molecule));
            }
            thread_local ShardedWords::Scratch scratch;
            dest.words.finish(shards,owner,scratch);
            dest.targets.finish(shards,target_owner,scratch);
            dest.occurrences.finish(shards,target_owner,scratch);
        });
        // Each identity has one owner. The all-read-minimizer histogram is
        // unchanged; panel target counts have their own exact-match reduction.
        pool.for_each(shards,[&](size_t shard) {
            auto&counts=all_read_counts[shard];
            for(const auto&task:tasks) {
                const auto&v=task.words[shard];
                for(size_t a=0;a<v.size();) {
                    size_t b=a+1;while(b<v.size()&&v[b]==v[a])++b;
                    const uint64_t add=b-a;auto&count=counts[v[a]];
                    if(add>UINT64_MAX-count)throw std::overflow_error("Calibration read count");
                    count+=add;
                    a=b;
                }
                const auto&matched=task.targets[shard];
                for(size_t a=0;a<matched.size();) {
                    size_t b=a+1;while(b<matched.size()&&matched[b]==matched[a])++b;
                    const uint64_t add=b-a;
                    auto&target_count=index.targets[size_t(matched[a])].count;
                    if(add>UINT64_MAX-target_count)throw std::overflow_error("Read count");
                    target_count+=add;a=b;
                }
                const auto&all=task.occurrences[shard];
                for(size_t a=0;a<all.size();) {
                    size_t b=a+1;while(b<all.size()&&all[b]==all[a])++b;
                    auto&n=index.targets[size_t(all[a])].occurrence_count;
                    if(uint64_t(b-a)>UINT64_MAX-n)throw std::overflow_error("Read occurrence count");
                    n+=b-a;a=b;
                }
            }
        });
        for(const auto&task:tasks) {
            if(task.target_votes>UINT64_MAX-result.target_end_observations||
               task.minimizer_votes>UINT64_MAX-result.minimizer_target_end_observations)
                throw std::overflow_error("Target observation audit");
            result.target_end_observations+=task.target_votes;
            result.minimizer_target_end_observations+=task.minimizer_votes;
        }
        for(auto&task:tasks)for(auto&marker:task.markers)
            result.markers.push_back(std::move(marker));
        // The audit sink is called serially in physical input order, after
        // parallel localization. It may be a compressed, streaming writer.
        if(audit)for(const auto&task:tasks)for(const auto&record:task.audit)audit(record);
        if(progress)progress(result.pairs,result.markers.size());
    }
    if(!contacts.empty())throw std::runtime_error("Contact IDs are not a subset of these raw FASTQs");
    std::vector<std::map<uint64_t,uint64_t>> histograms(shards);
    pool.for_each(shards,[&](size_t shard) {
        for(auto[word,count]:all_read_counts[shard]){(void)word;++histograms[shard][count];}
        FlatMap<uint64_t,uint64_t>().swap(all_read_counts[shard]);
    });
    for(const auto&hist:histograms)for(auto[count,n]:hist)result.count_histogram[count]+=n;
    return result;
}
}
