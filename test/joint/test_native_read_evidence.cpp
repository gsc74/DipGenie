#include "TestTemp.hpp"
#include "NativeReadEvidence.hpp"
#include <cassert>
#include <iostream>
#include <random>
#include <thread>
#include <tuple>
#include <unistd.h>

using namespace hic;
static std::string dna(size_t n,uint64_t seed) {
    std::mt19937_64 rng(seed);std::string s(n,'A');for(auto&c:s)c="ACGT"[rng()%4];return s;
}
static void write_fastq(const std::filesystem::path&p,const std::vector<std::string>&reads,int side,bool duplicate=false) {
    std::ofstream out(p);
    for(size_t i=0;i<reads.size();++i)out<<"@pair"<<(duplicate?0:i)<<'/'<<side<<'\n'<<reads[i]<<"\n+\n"<<std::string(reads[i].size(),'I')<<'\n';
    out.close();assert(out);
}
static SequenceIndex target_index(const std::vector<std::string>&sequences,unsigned k,unsigned w) {
    SequenceIndex index;index.k=k;index.w=w;std::set<uint64_t> selected;
    for(const auto&s:sequences)minimizers(s,k,w,[&](auto word,size_t){selected.insert(word);});
    for(auto word:selected){index.ids[word]=index.targets.size();WordTarget t;t.word=word;index.targets.push_back(t);}
    return index;
}
int main() {
    NativeLocalizer::Policy policy;policy.k=15;policy.w=7;
    const auto ref=dna(12000,1234),alt=dna(800,2345),other=dna(4000,3456);
    NativeGraph g;g.names={"a","b","alt"};g.labels={ref.substr(0,6000),ref.substr(6000),alt};
    g.paths={{"CHM13",{0,2}},{"alternate",{0,4,2}}};
    for(auto [a,b]:{std::pair<uint32_t,uint32_t>{0,2},{0,4},{4,2}}) {
        g.edges.insert(NativeGraph::edge(a,b));g.edges.insert(NativeGraph::edge(b^1,a^1));
    }
    auto contexts=native_graph_contexts(g,0);
    PlacementSequence outside;outside.name="other";outside.sequence=other;outside.component=1;outside.space=100;outside.reference_contig=1;
    contexts.push_back(std::move(outside));
    NativeLocalizer localizer(contexts,policy,4);
    auto targets=target_index({g.spell(0),g.spell(1),other},policy.k,policy.w);
    std::vector<std::string> one={ref.substr(200,150),ref.substr(200,150),alt.substr(100,150),
        ref.substr(4500,80)+other.substr(500,70),ref.substr(200,150),std::string(150,'N')};
    std::vector<std::string> two={ref.substr(1000,150),ref.substr(9000,150),alt.substr(400,150),
        ref.substr(1000,150),other.substr(200,150),ref.substr(1000,150)};
    std::string pattern_buffer=test_temp_pattern("dipgenie-native-evidence");char*pattern=pattern_buffer.data();char*tmp=mkdtemp(pattern);assert(tmp);
    const std::filesystem::path directory(tmp);
    // Independent gzip readers must see exactly the original byte stream,
    // including empty streams, member boundaries, oversized records and a
    // short final batch. Compression scheduling cannot reorder the audit.
    auto inflate=[](const std::filesystem::path&path) {
        auto f=gzopen(path.c_str(),"rb");assert(f);std::string result;char block[8192];int n;
        while((n=gzread(f,block,sizeof(block)))>0)result.append(block,size_t(n));
        assert(n==0&&gzclose(f)==Z_OK);return result;
    };
    const std::string payload=dna((33u<<20)+71,9941)+"\n";
    std::string gzip_bytes;
    for(int threads:{1,2,8}) {
        const auto path=directory/("audit-"+std::to_string(threads)+".gz");
        NativeObservationAudit writer(path,threads);
        writer.write(std::string_view(payload).substr(0,19));
        writer.write(std::string_view(payload).substr(19));writer.close();writer.close();
        assert(inflate(path)==payload);
        std::ifstream packed(path,std::ios::binary);std::string bytes((std::istreambuf_iterator<char>(packed)),{});
        if(gzip_bytes.empty())gzip_bytes=bytes;else assert(bytes==gzip_bytes);
        bool failed=false;try{writer.write("closed");}catch(const std::logic_error&){failed=true;}assert(failed);
    }
    {const auto path=directory/"empty-audit.gz";NativeObservationAudit writer(path,2);writer.close();assert(inflate(path).empty());}
    if(std::filesystem::exists("/dev/full")) {
        bool failed=false;try{NativeObservationAudit writer("/dev/full",2);writer.write("record");writer.close();}
        catch(const std::runtime_error&){failed=true;}assert(failed);
    }
    std::vector<ReadLibrary> libraries;
    for(const auto&kind:{"SR","HIC"}) {
        const auto prefix=directory/kind;
        write_fastq(prefix.string()+".1.fq",one,1);write_fastq(prefix.string()+".2.fq",two,2);
        libraries.push_back({kind,kind,prefix.string()+".1.fq",prefix.string()+".2.fq"});
    }
    auto expected=targets;
    const auto baseline=read_evidence(expected,libraries[0].r1,libraries[0].r2,{},true,1);
    assert(baseline.pairs==6);
    std::string canonical_audit;
    for(int threads:{1,2,4,8}) {
        auto index=targets;std::string audit;
        const auto coordinator=std::this_thread::get_id();
        using Event=std::tuple<std::string,uint64_t,uint64_t,bool>;
        std::vector<Event> events;
        auto native=read_native_library_evidence(index,libraries,localizer,threads,{},
            [&](std::string_view row){audit.append(row);},
            [&](const ReadLibrary&l,uint64_t pairs,uint64_t linked,bool complete) {
                assert(std::this_thread::get_id()==coordinator);
                events.emplace_back(l.name,pairs,linked,complete);
            });
        const std::vector<Event> expected_events={{"SR",0,0,false},{"SR",6,2,false},{"SR",6,2,true},
            {"HIC",0,0,false},{"HIC",6,4,false},{"HIC",6,4,true}};
        assert(events==expected_events);
        auto silent_index=targets;std::string silent_audit;
        const auto silent=read_native_library_evidence(silent_index,libraries,localizer,threads,{},
            [&](std::string_view row){silent_audit.append(row);});
        assert(audit==silent_audit&&native.combined.count_histogram==silent.combined.count_histogram);
        assert(native.combined.markers.size()==silent.combined.markers.size());
        for(size_t x=0;x<index.targets.size();++x)assert(index.targets[x].count==silent_index.targets[x].count);
        assert(native.counts.size()==2&&native.combined.pairs==12&&native.combined.markers.size()==6);
        assert((native.contact_pairs==std::vector<uint64_t>{6,6}));
        assert((native.marker_pairs==std::vector<uint64_t>{2,4}));
        assert(native.combined.target_end_observations==2*baseline.target_end_observations);
        assert(native.combined.minimizer_target_end_observations==2*baseline.minimizer_target_end_observations);
        for(auto [count,n]:baseline.count_histogram)assert(native.combined.count_histogram.at(count)==2*n);
        for(size_t x=0;x<index.targets.size();++x) {
            assert(index.targets[x].count==2*expected.targets[x].count);
            for(const auto&lib:native.counts) {
                assert(lib.counts[x]==expected.targets[x].count);
                assert(lib.occurrences[x]==expected.targets[x].occurrence_count);
            }
        }
        size_t sr=0,hic=0,unknown=0;
        for(const auto&m:native.combined.markers) {
            if(m.kind==MoleculeKind::ShortRead){++sr;assert(m.library==0&&m.id!="pair1");}
            else{++hic;assert(m.library==1);}
            if(m.geometry==MoleculeGeometry::Unknown){++unknown;assert(m.reference_distance==unknown_contact_distance&&m.id=="pair2");}
            else assert(m.geometry==MoleculeGeometry::Known&&m.reference_distance>0);
            // The split Hi-C read (pair3) is resolved to its 5' placement with partial support.
            assert(m.localization_qc==(m.kind==MoleculeKind::HiC&&m.id=="pair3"?2:1));
        }
        assert(sr==2&&hic==4&&unknown==2);
        assert(audit.find("sr_anchor_separation")!=std::string::npos&&audit.find("unprojected")!=std::string::npos);
        assert(audit.find("cross_component")!=std::string::npos&&audit.find("split")!=std::string::npos);
        assert(audit.find("five_prime,")!=std::string::npos);
        if(canonical_audit.empty())canonical_audit=audit;else assert(audit==canonical_audit);
    }
    // Exact identities, selected anchors and association survive reversal of
    // either mate and mate exchange. Hi-C has no FR or insert-size gate.
    for(size_t pair:{size_t(0),size_t(1),size_t(2)}) {
        FastqRead a,b;a.sequence=one[pair];b.sequence=two[pair];
        const auto base=native_pair_evidence(targets,localizer,a,b,"id",MoleculeKind::HiC,0);
        for(bool reverse1:{false,true})for(bool reverse2:{false,true})for(bool swap:{false,true}) {
            a.sequence=reverse1?reverse_complement(one[pair]):one[pair];
            b.sequence=reverse2?reverse_complement(two[pair]):two[pair];
            if(swap)std::swap(a,b);
            auto test=native_pair_evidence(targets,localizer,a,b,"id",MoleculeKind::HiC,0);
            assert(test.marker.reference_distance==base.marker.reference_distance&&test.marker.geometry==base.marker.geometry);
            for(int e=0;e<2;++e)assert(test.marker.ends[e]==base.marker.ends[swap?1-e:e]);
        }
    }
    // Library-scoped identity: the same physical name in the two different
    // libraries above is valid; repeating it within one library is rejected.
    auto duplicate=libraries[0];
    duplicate.r1=(directory/"duplicate.1.fq").string();duplicate.r2=(directory/"duplicate.2.fq").string();
    write_fastq(duplicate.r1,one,1,true);write_fastq(duplicate.r2,two,2,true);
    bool rejected=false;try{auto index=targets;read_native_library_evidence(index,{duplicate},localizer,2);}catch(const std::invalid_argument&){rejected=true;}assert(rejected);
    // Truncation and conflicting mate labels fail before a successful result.
    auto malformed=libraries[0];malformed.r2=(directory/"short.2.fq").string();
    write_fastq(malformed.r2,{two[0]},2);
    rejected=false;try{auto index=targets;read_native_library_evidence(index,{malformed},localizer,2);}catch(const std::runtime_error&){rejected=true;}assert(rejected);
    write_fastq(malformed.r2,two,1);
    rejected=false;try{auto index=targets;read_native_library_evidence(index,{malformed},localizer,2);}catch(const std::runtime_error&){rejected=true;}assert(rejected);
    // A multi-batch stream remains input ordered, including a short final
    // batch. Truncation after a successful batch must not announce completion.
    auto stream=libraries[0];stream.r1=(directory/"stream.1.fq").string();stream.r2=(directory/"stream.2.fq").string();
    std::vector<std::string> unplaced(4100,std::string(20,'N'));
    write_fastq(stream.r1,unplaced,1);write_fastq(stream.r2,unplaced,2);
    using Counts=std::tuple<uint64_t,uint64_t,bool>;
    std::vector<Counts> counts;
    const LibraryProgressSink observe=[&](const ReadLibrary&,uint64_t p,uint64_t linked,bool done){counts.emplace_back(p,linked,done);};
    auto stream_index=targets;
    const auto streamed=read_native_library_evidence(stream_index,{stream},localizer,4,{}, {},observe);
    const std::vector<Counts> expected_counts={{0,0,false},{4096,0,false},{4100,0,false},{4100,0,true}};
    assert(counts==expected_counts&&streamed.combined.pairs==4100&&streamed.combined.markers.empty());
    unplaced.pop_back();write_fastq(stream.r2,unplaced,2);counts.clear();
    rejected=false;try{auto index=targets;read_native_library_evidence(index,{stream},localizer,4,{}, {},observe);}
    catch(const std::runtime_error&){rejected=true;}
    assert(rejected&&counts==std::vector<Counts>({{0,0,false},{4096,0,false}}));
    // Empty gzip files have nonzero file size, but are not valid libraries.
    for(auto*p:{&stream.r1,&stream.r2}) {
        *p+=".gz";auto gz=gzopen(p->c_str(),"wb");assert(gz);assert(gzclose(gz)==Z_OK);
    }
    counts.clear();rejected=false;
    try{auto index=targets;read_native_library_evidence(index,{stream},localizer,4,{}, {},observe);}
    catch(const std::invalid_argument&){rejected=true;}
    assert(rejected&&counts==std::vector<Counts>({{0,0,false}}));
    std::filesystem::remove_all(directory);
    std::cout<<"PASS native raw SR/Hi-C: strict identity, source-specific admission, unknown geometry, unchanged exact dosage, durable audit and 1/2/4/8-thread parity\n";
}
