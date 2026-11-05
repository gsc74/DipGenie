#pragma once
// Schema-7 graph payload adapter. The old hashed positional anchors are NOT
// imported: SequenceIndex must rebuild the active lexicographic k31/w25 data.
// Component node labels already encode traversal orientation (schema contract).
#include "NativeGraph.hpp"
#include <array>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>

namespace hic::gindex {
struct Metadata {
    uint32_t id=0; uint64_t length=0,vertices=0,edges=0,paths=0,bytes=0,checksum=0;
    std::string chromosome,file;
};
inline std::vector<Metadata> components(const std::filesystem::path&root) {
    std::ifstream in(root/"manifest.dgi");std::string line;
    if(!std::getline(in,line)||line!="DIPGENIE_GRAPH_INDEX\t7")
        throw std::runtime_error("Expected schema-7 gindex");
    std::vector<Metadata> records;
    while(std::getline(in,line)) {
        std::istringstream row(line);std::string kind;row>>kind;
        // Benchmark-specific sample exclusion belongs in input preparation,
        // not in a general-purpose graph loader.
        if(kind!="component")continue;
        Metadata result;std::string checksum;
        if(!(row>>result.id>>std::quoted(result.chromosome)>>result.length>>result.vertices>>result.edges
             >>result.paths>>std::quoted(result.file)>>result.bytes>>checksum))
            throw std::runtime_error("Malformed graph-index component record");
        size_t used=0;result.checksum=std::stoull(checksum,&used,16);
        if(used!=checksum.size()||std::filesystem::path(result.file).filename()!=result.file)
            throw std::runtime_error("Invalid graph shard filename/checksum");
        for(const auto&previous:records)if(previous.id==result.id||previous.chromosome==result.chromosome)
            throw std::runtime_error("Duplicate graph component");
        records.push_back(std::move(result));
    }
    if(records.empty())throw std::runtime_error("Empty graph index");
    return records;
}
inline Metadata metadata(const std::filesystem::path&root) {
    auto records=components(root);
    if(records.size()!=1)throw std::runtime_error("Select a component from this multi-component index");
    return records[0];
}
class Reader {
    std::ifstream in;uint64_t left,hash=14695981039346656037ULL,expected;
public:
    Reader(const std::filesystem::path&p,uint64_t size,uint64_t checksum):in(p,std::ios::binary),left(size),expected(checksum) {
        if(!in||std::filesystem::file_size(p)!=size)throw std::runtime_error("Graph payload size mismatch");
    }
    void raw(void*p,size_t n) {
        if(n>left||!in.read(static_cast<char*>(p),std::streamsize(n)))throw std::runtime_error("Truncated graph payload");
        auto*q=static_cast<unsigned char*>(p);for(size_t i=0;i<n;++i){hash^=q[i];hash*=1099511628211ULL;}left-=n;
    }
    template<class T>T number(){static_assert(std::is_unsigned<T>::value);std::array<unsigned char,sizeof(T)> b{};raw(b.data(),b.size());T x=0;for(size_t i=0;i<b.size();++i)x|=T(b[i])<<(8*i);return x;}
    uint64_t count(size_t width=1){auto n=number<uint64_t>();if(n>left/width)throw std::runtime_error("Payload count exceeds remaining bytes");return n;}
    std::string string(){auto n=count();std::string s(n,'\0');raw(s.data(),s.size());return s;}
    void discard(size_t n){std::array<char,65536>b{};while(n){auto c=std::min(n,b.size());raw(b.data(),c);n-=c;}}
    template<class T>std::vector<T> flat(){auto n=count(sizeof(T));std::vector<T> v(n);const uint16_t one=1;if(*reinterpret_cast<const unsigned char*>(&one)==1)raw(v.data(),n*sizeof(T));else for(auto&x:v)x=number<T>();return v;}
    template<class T>std::vector<std::vector<T>> nested(){auto n=count(8);std::vector<std::vector<T>>v(n);for(auto&x:v)x=flat<T>();return v;}
    void finish(){if(left||hash!=expected)throw std::runtime_error("Graph payload checksum/trailing-byte mismatch");}
};
inline NativeGraph load(const std::filesystem::path&root,const Metadata&meta) {
    Reader in(root/meta.file,meta.bytes,meta.checksum);
    std::array<char,8> magic{};in.raw(magic.data(),magic.size());
    if(std::string(magic.data(),8)!="DGCIDX01"||in.number<uint32_t>()!=3||in.number<uint32_t>()!=meta.id
       ||in.number<uint64_t>()!=meta.length||in.number<uint8_t>()!=1||in.string()!=meta.chromosome)
        throw std::runtime_error("Graph payload header disagrees with manifest");
    auto adjacency=in.nested<uint32_t>();auto lengths=in.flat<uint32_t>();NativeGraph g;
    auto n=in.count(8);if(n!=meta.vertices||adjacency.size()!=n||lengths.size()!=n||n>UINT32_MAX/2)
        throw std::runtime_error("Graph vertex dimensions");
    g.labels.reserve(n);g.names.reserve(n);
    for(size_t v=0;v<n;++v){g.labels.push_back(in.string());g.names.push_back(std::to_string(v));if(g.labels.back().size()!=lengths[v])throw std::runtime_error("Node sequence length mismatch");}
    auto orientations=in.flat<uint8_t>();if(!orientations.empty()&&orientations.size()!=n)throw std::runtime_error("Node orientation dimensions");
    auto walks=in.nested<uint32_t>();auto strands=in.nested<uint8_t>();
    if(walks.size()!=meta.paths||(!strands.empty()&&strands.size()!=walks.size()))throw std::runtime_error("Walk dimensions");
    auto names=in.count(8);if(names!=walks.size())throw std::runtime_error("Path name dimensions");
    g.paths.resize(names);
    for(size_t h=0;h<names;++h){g.paths[h].name=in.string();g.paths[h].walk=std::move(walks[h]);if(!strands.empty()&&strands[h].size()!=g.paths[h].walk.size())throw std::runtime_error("Path strand dimensions");}
    if(in.count(8)!=names)throw std::runtime_error("Path sequence dimensions");
    // These optional spelling caches are redundant, not the graph input.
    for(size_t h=0;h<names;++h){auto bytes=in.count();in.discard(bytes);}
    if(in.count(32)!=names)throw std::runtime_error("Anchor catalog dimensions");
    for(size_t h=0;h<names;++h){(void)in.string();(void)in.number<uint64_t>();(void)in.number<uint64_t>();(void)in.number<uint64_t>();}
    in.finish();uint64_t edges=0;
    g.edges.reserve(size_t(meta.edges)*2);
    for(uint32_t v=0;v<n;++v){auto&row=adjacency[v];std::sort(row.begin(),row.end());edges+=row.size();for(auto w:row){if(w>=n)throw std::runtime_error("Graph edge leaves component");g.edges.insert(NativeGraph::edge(2*v,2*w));g.edges.insert(NativeGraph::edge(2*w+1,2*v+1));}}
    if(edges!=meta.edges)throw std::runtime_error("Edge count differs from manifest");
    uint64_t reference_length=0;
    for(size_t h=0;h<g.paths.size();++h){auto&p=g.paths[h];for(size_t i=0;i<p.walk.size();++i){auto v=p.walk[i];if(v>=n)throw std::runtime_error("Walk leaves component");if(i&&!std::binary_search(adjacency[p.walk[i-1]].begin(),adjacency[p.walk[i-1]].end(),v))throw std::runtime_error("Walk has no graph edge");if(!h)reference_length+=g.labels[v].size();}for(auto&v:p.walk)v*=2;}
    if(reference_length!=meta.length)throw std::runtime_error("Reference spelling length mismatch");
    return g;
}
inline NativeGraph load(const std::filesystem::path&root){return load(root,metadata(root));}
}
