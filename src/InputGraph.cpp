#include "InputGraph.hpp"
#include "PanGenieVcf.hpp"
#include <numeric>
#include <set>

namespace hic {
namespace {
namespace fs=std::filesystem;
std::vector<ComponentReference> read_reference(const fs::path&path) {
    TextInput in(path.string());std::string line;std::vector<ComponentReference> refs;
    std::set<std::string> names;
    while(in.getline(line)) {
        if(line.empty())continue;
        if(line[0]=='>') {
            auto name=line.substr(1,line.find_first_of(" \t")-1);
            if(name.empty()||!names.insert(name).second)throw std::runtime_error("Invalid reference dictionary");
            refs.push_back({name,{}});
        } else {
            if(refs.empty())throw std::runtime_error("Reference FASTA header missing");
            for(unsigned char c:line)refs.back().sequence+=char(std::toupper(c));
        }
    }
    if(refs.empty())throw std::runtime_error("Empty reference");
    for(const auto&r:refs)if(r.sequence.empty())throw std::runtime_error("Empty reference record");
    return refs;
}
uint64_t checksum(const fs::path&path) {
    std::ifstream in(path,std::ios::binary);uint64_t h=14695981039346656037ULL;
    if(!in)throw std::runtime_error("Missing graph shard");
    std::array<char,65536> buffer{};
    while(in.read(buffer.data(),buffer.size())||in.gcount())
        for(std::streamsize i=0;i<in.gcount();++i){h^=static_cast<unsigned char>(buffer[size_t(i)]);h*=1099511628211ULL;}
    if(!in.eof())throw std::runtime_error("Graph shard read failed");
    return h;
}
NativeGraph vcf_graph(const std::string&vcf,const std::string&reference) {
    auto panel=pangenie_vcf::read_panel(vcf,reference);
    std::unique_ptr<gfa_t,decltype(&gfa_destroy)> owner(panel.graph,gfa_destroy);
    const auto*g=owner.get();NativeGraph result;
    if(g->n_seg>UINT32_MAX/2)throw std::overflow_error("VCF graph vertex capacity");
    result.labels.reserve(g->n_seg);result.names.reserve(g->n_seg);
    for(uint64_t i=0;i<g->n_seg;++i) {
        result.labels.emplace_back(g->seg[i].seq);result.names.emplace_back(g->seg[i].name);
    }
    result.edges.reserve(g->n_arc);
    for(uint64_t i=0;i<g->n_arc;++i)if(!g->arc[i].del)
        result.edges.insert(NativeGraph::edge(uint32_t(g->arc[i].v),uint32_t(g->arc[i].w)));
    for(uint32_t h=0;h<g->n_walk;++h) {
        const auto&w=g->walk[h];NativeGraph::Path p;
        if(w.snid<0||uint32_t(w.snid)>=g->n_sseq)throw std::runtime_error("VCF contig identity");
        const std::string contig=g->sseq[w.snid].name;
        p.name=std::string(w.sample)=="REF"?contig:std::string(w.sample)+"#"+std::to_string(w.hap)+"#"+contig;
        for(uint64_t i=0;i<w.n_v;++i) {
            if(w.v[i]>=2*g->n_seg)throw std::runtime_error("VCF path outside graph");
            p.walk.push_back(uint32_t(w.v[i]));
        }
        result.paths.push_back(std::move(p));
    }
    return result;
}
GraphInput split(NativeGraph g,const std::set<std::string>&reference_names,const std::string&prefix) {
    const size_t V=g.labels.size();std::vector<uint32_t> parent(V);
    std::iota(parent.begin(),parent.end(),0);
    auto root=[&](uint32_t v){while(parent[v]!=v){parent[v]=parent[parent[v]];v=parent[v];}return v;};
    for(auto e:g.edges){const auto a=root(uint32_t(e>>32)/2),b=root(uint32_t(e)/2);if(a!=b)parent[b]=a;}
    std::map<uint32_t,size_t> component;std::vector<NativeGraph> parts;
    std::vector<uint32_t> local(V);std::vector<size_t> membership(V);
    for(uint32_t v=0;v<V;++v) {
        auto[it,new_component]=component.emplace(root(v),parts.size());if(new_component)parts.emplace_back();
        auto&part=parts[it->second];membership[v]=it->second;local[v]=uint32_t(part.labels.size());
        part.labels.push_back(std::move(g.labels[v]));part.names.push_back(std::move(g.names[v]));
    }
    for(auto e:g.edges) {
        uint32_t a=uint32_t(e>>32),b=uint32_t(e);
        if(membership[a/2]!=membership[b/2])throw std::logic_error("Graph component split");
        parts[membership[a/2]].edges.insert(NativeGraph::edge(2*local[a/2]+(a&1),2*local[b/2]+(b&1)));
    }
    for(auto&p:g.paths) {
        if(p.walk.empty())throw std::runtime_error("Empty panel path");
        const auto c=membership[p.walk[0]/2];
        for(auto&v:p.walk){if(membership[v/2]!=c)throw std::runtime_error("Path spans disconnected components");v=2*local[v/2]+(v&1);}
        parts[c].paths.push_back(std::move(p));
    }
    auto storage=std::make_shared<std::vector<NativeGraph>>();GraphInput input;
    for(auto&part:parts) {
        if(part.paths.empty())continue; // No panel walks: not an inference component.
        std::vector<size_t> hits;
        for(size_t h=0;h<part.paths.size();++h) {
            const auto&name=part.paths[h].name;
            if(reference_names.count(name)||(!prefix.empty()&&(name==prefix||name.rfind(prefix+"#",0)==0)))hits.push_back(h);
        }
        if(hits.size()!=1)throw std::runtime_error("Each component requires exactly one complete --REF path");
        std::swap(part.paths[0],part.paths[hits[0]]);
        input.references.push_back({part.paths[0].name,part.spell(0)});storage->push_back(std::move(part));
    }
    if(storage->empty())throw std::runtime_error("No graph components with panel walks");
    input.load=[storage](size_t c){return std::move(storage->at(c));};
    return input;
}
}

GraphInput graph_input(const std::map<std::string,std::string>&args,int threads) {
    if(threads<1)throw std::invalid_argument("Positive graph-input thread count required");
    if(args.count("--gfa")) {
        if(!args.count("--REF"))throw std::invalid_argument("GFA input requires --REF path/sample prefix");
        return split(NativeGraph::load(args.at("--gfa"),threads),{},args.at("--REF"));
    }
    if(args.count("--vcf")) {
        if(!args.count("--reference"))throw std::invalid_argument("VCF input requires --reference FASTA");
        std::set<std::string> names;for(const auto&r:read_reference(args.at("--reference")))names.insert(r.name);
        return split(vcf_graph(args.at("--vcf"),args.at("--reference")),names,{});
    }
    const fs::path root=args.at("--gindex");std::ifstream manifest(root/"manifest.dgi");std::string line;
    if(!std::getline(manifest,line))throw std::runtime_error("Missing graph manifest");
    GraphInput input;input.references=read_reference(root/"reference.fa");
    if(line=="DIPGENIE_GRAPH_INDEX\t7") {
        const auto records=gindex::components(root);
        if(records.size()!=input.references.size())throw std::runtime_error("Indexed reference/component count mismatch");
        // Legacy dictionary order is not assumed to match component order.
        std::vector<ComponentReference> ordered;
        for(const auto&meta:records) {
            auto p=std::find_if(input.references.begin(),input.references.end(),[&](const auto&r){return r.name==meta.chromosome;});
            if(p==input.references.end()||p->sequence.size()!=meta.length)throw std::runtime_error("Legacy reference/component mismatch");
            ordered.push_back(std::move(*p));
        }
        input.references=std::move(ordered);
        input.load=[root,records](size_t c){return gindex::load(root,records.at(c));};
    } else if(line=="DIPGENIE_LINKAGE_INDEX\t1") {
        struct Shard{std::string file,reference;uint64_t bytes=0,hash=0;};std::vector<Shard> shards;
        while(std::getline(manifest,line)) {
            if(line.empty())continue;
            std::istringstream row(line);std::string kind;Shard shard;
            if(!(row>>kind>>std::quoted(shard.file)>>std::quoted(shard.reference)>>shard.bytes>>std::hex>>shard.hash)||kind!="component"
                ||shard.file.empty()||fs::path(shard.file).filename()!=shard.file)
                throw std::runtime_error("Invalid linkage graph manifest");
            if(shards.size()>=input.references.size()||shard.reference!=input.references[shards.size()].name)
                throw std::runtime_error("Reference/component order mismatch");
            shards.push_back(std::move(shard));
        }
        if(shards.size()!=input.references.size())throw std::runtime_error("Reference/component count mismatch");
        input.load=[root,shards,threads](size_t c) {
            const auto&s=shards.at(c);const auto path=root/s.file;
            if(fs::file_size(path)!=s.bytes||checksum(path)!=s.hash)throw std::runtime_error("Graph shard checksum mismatch");
            auto g=NativeGraph::load(path.string(),threads);
            if(g.paths.at(0).name!=s.reference)throw std::runtime_error("Indexed REF path mismatch");
            return g;
        };
    } else throw std::runtime_error("Unsupported graph index schema");
    // A modified reference must never silently route reads against sequence
    // different from the graph spelling. Validate again when loading a shard.
    std::vector<uint64_t> expected;
    for(const auto&r:input.references) {
        uint64_t h=14695981039346656037ULL;
        for(unsigned char c:r.sequence){h^=c;h*=1099511628211ULL;}
        expected.push_back(h);
    }
    auto load=input.load;
    input.load=[load,expected](size_t c) {
        auto graph=load(c);uint64_t h=14695981039346656037ULL;
        for(unsigned char b:graph.spell(0)){h^=b;h*=1099511628211ULL;}
        if(h!=expected.at(c))throw std::runtime_error("Reference/graph spelling checksum mismatch");
        return graph;
    };
    return input;
}

void write_graph_index(const GraphInput&input,const fs::path&dir) {
    if(!fs::create_directory(dir))throw std::runtime_error("Index directory exists");
    std::ofstream refs(dir/"reference.fa"),manifest(dir/"manifest.pending");
    manifest<<"DIPGENIE_LINKAGE_INDEX\t1\n";
    for(size_t c=0;c<input.references.size();++c) {
        const auto&r=input.references[c];auto g=input.load(c);
        if(g.spell(0)!=r.sequence)throw std::runtime_error("Reference/index spelling mismatch");
        refs<<'>'<<r.name<<'\n';for(size_t p=0;p<r.sequence.size();p+=80)refs<<r.sequence.substr(p,80)<<'\n';
        const std::string filename="component-"+std::to_string(c)+".gfa.gz";const auto path=dir/filename;
        gzFile out=gzopen(path.c_str(),"wb1");if(!out)throw std::runtime_error("Cannot create graph shard");
        auto emit=[&](const std::string&s){size_t p=0;while(p<s.size()){const auto n=std::min(size_t(1<<20),s.size()-p);if(gzwrite(out,s.data()+p,unsigned(n))!=int(n))throw std::runtime_error("Graph shard write failed");p+=n;}};
        try {
            emit("H\tVN:Z:1.0\n");
            for(size_t v=0;v<g.labels.size();++v)emit("S\t"+std::to_string(v)+"\t"+g.labels[v]+"\n");
            std::vector<uint64_t> edges(g.edges.begin(),g.edges.end());std::sort(edges.begin(),edges.end());
            for(auto e:edges) {
                uint32_t a=uint32_t(e>>32),b=uint32_t(e);if(e>NativeGraph::edge(b^1,a^1))continue;
                emit("L\t"+std::to_string(a/2)+"\t"+(a&1?"-":"+")+"\t"+std::to_string(b/2)+"\t"+(b&1?"-":"+")+"\t0M\n");
            }
            for(const auto&p:g.paths) {
                emit("P\t"+p.name+"\t");size_t i=0;
                for(auto v:p.walk)emit(std::string(i++?",":"")+std::to_string(v/2)+(v&1?"-":"+"));
                emit("\t*\n");
            }
        }catch(...){gzclose(out);throw;}
        if(gzclose(out)!=Z_OK)throw std::runtime_error("Graph shard close failed");
        manifest<<"component "<<std::quoted(filename)<<' '<<std::quoted(r.name)<<' '<<fs::file_size(path)<<' '<<std::hex<<checksum(path)<<std::dec<<'\n';
    }
    refs.close();manifest.close();if(!refs||!manifest)throw std::runtime_error("Graph index write failed");
    fs::rename(dir/"manifest.pending",dir/"manifest.dgi");
}
}
