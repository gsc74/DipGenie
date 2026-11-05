#include "TestTemp.hpp"
#include "NativeGraph.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <unistd.h>

namespace {
using hic::NativeGraph;
void require(bool yes,const std::string&why){if(!yes)throw std::runtime_error(why);}

// Frozen pre-optimization parser, kept as a differential test oracle only.
// Preserve node insertion order, every oriented step, reverse edges and
// rejection behavior, not only the spelled strings.
NativeGraph reference_load(const std::string&file) {
    NativeGraph g;std::unordered_map<std::string,uint32_t> ids;
    auto node=[&](const std::string&name) {
        auto[it,inserted]=ids.emplace(name,uint32_t(ids.size()));
        if(inserted){if(ids.size()>UINT32_MAX/2)throw std::runtime_error("Too many oriented vertices");g.names.push_back(name);g.labels.emplace_back();}
        return it->second;
    };
    hic::TextInput in(file);std::string line;
    while(in.getline(line)) {
        if(line.empty()||line[0]=='#')continue;
        auto f=hic::fields(line);
        if(f[0]=="S") {
            if(f.size()<3||f[2]=="*")throw std::runtime_error("GFA requires explicit segment sequences");
            auto n=node(f[1]);if(!g.labels[n].empty())throw std::runtime_error("Duplicate segment");
            g.labels[n]=f[2];for(char&c:g.labels[n])c=std::toupper(static_cast<unsigned char>(c));
        } else if(f[0]=="L") {
            if(f.size()<6||f[5]!="0M"||(f[2]!="+"&&f[2]!="-")||(f[4]!="+"&&f[4]!="-"))
                throw std::runtime_error("GFA links must use explicit zero overlap and strand");
            uint32_t a=2*node(f[1])+(f[2]=="-"),b=2*node(f[3])+(f[4]=="-");
            g.edges.insert(NativeGraph::edge(a,b));g.edges.insert(NativeGraph::edge(b^1,a^1));
        } else if(f[0]=="P"||f[0]=="W") {
            NativeGraph::Path p;
            if(f[0]=="P") {
                if(f.size()<3)throw std::runtime_error("Truncated P record");
                p.name=f[1];
                for(auto&s:hic::fields(f[2],',')) {
                    if(s.size()<2||(s.back()!='+'&&s.back()!='-'))throw std::runtime_error("Invalid oriented path step");
                    p.walk.push_back(2*node(s.substr(0,s.size()-1))+(s.back()=='-'));
                }
            } else {
                if(f.size()<7)throw std::runtime_error("Truncated W record");
                p.name=f[1]+"#"+f[2]+"#"+f[3];size_t i=0;
                while(i<f[6].size()) {
                    char strand=f[6][i++];if(strand!='>'&&strand!='<')throw std::runtime_error("Invalid W orientation");
                    size_t end=f[6].find_first_of("><",i);if(end==std::string::npos)end=f[6].size();
                    if(end==i)throw std::runtime_error("Empty W segment");
                    p.walk.push_back(2*node(f[6].substr(i,end-i))+(strand=='<'));i=end;
                }
            }
            if(p.walk.empty())throw std::runtime_error("Empty panel path");
            g.paths.push_back(std::move(p));
        }
    }
    if(g.paths.empty())throw std::runtime_error("No panel paths");
    for(auto&s:g.labels)if(s.empty())throw std::runtime_error("Undefined graph segment");
    for(auto&p:g.paths)for(size_t i=1;i<p.walk.size();++i)
        if(!g.edges.count(NativeGraph::edge(p.walk[i-1],p.walk[i])))throw std::runtime_error("Panel step missing from GFA links: "+p.name);
    return g;
}
void same(const NativeGraph&a,const NativeGraph&b) {
    require(a.names==b.names,"node names/order changed");
    require(a.labels==b.labels,"node sequences changed");
    require(a.edges==b.edges,"oriented edges changed");
    require(a.paths.size()==b.paths.size(),"path count changed");
    for(size_t h=0;h<a.paths.size();++h)
        require(a.paths[h].name==b.paths[h].name&&a.paths[h].walk==b.paths[h].walk,"oriented path/order changed");
}
struct Temp {
    std::filesystem::path dir;
    Temp(){std::string name_buffer=test_temp_pattern("dipgenie-graph-loader");char*name=name_buffer.data();auto p=mkdtemp(name);if(!p)throw std::runtime_error("mkdtemp");dir=p;}
    ~Temp(){std::error_code error;std::filesystem::remove_all(dir,error);}
};
void write(const std::filesystem::path&p,const std::string&s,bool zipped) {
    if(zipped) {
        gzFile f=gzopen(p.c_str(),"wb1");require(f!=nullptr,"gzip open");
        const int n=gzwrite(f,s.data(),unsigned(s.size()));const int rc=gzclose(f);
        require(n==int(s.size())&&rc==Z_OK,"gzip write");
    } else {std::ofstream f(p,std::ios::binary);f<<s;f.close();require(bool(f),"fixture write");}
}
std::string error(const std::string&file,int threads,bool reference) {
    try {if(reference)(void)reference_load(file);else(void)NativeGraph::load(file,threads);}
    catch(const std::exception&e){return e.what();}
    return {};
}
void check(const std::filesystem::path&p,const std::string&s,bool zipped,int threads=4) {
    write(p,s,zipped);
    const auto expected=error(p.string(),1,true);
    if(!expected.empty()) {require(error(p.string(),threads,false)==expected,"different rejection: "+expected);return;}
    const auto reference=reference_load(p.string());
    same(reference,NativeGraph::load(p.string(),1));same(reference,NativeGraph::load(p.string(),threads));
}
}

int main(int argc,char**argv)try {
    if(argc==4) {
        // Optional real-shard diagnostic: compare complete in-memory graphs,
        // not just a sample or a hash. Loader-only, not an inference run.
        const std::string file=argv[1],order=argv[3];const int threads=std::stoi(argv[2]);
        require(order=="new-first"||order=="old-first","invalid benchmark order");
        NativeGraph old,newer;double old_s=0,new_s=0;
        auto old_load=[&]{const auto start=std::chrono::steady_clock::now();old=reference_load(file);old_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
        auto new_load=[&]{const auto start=std::chrono::steady_clock::now();newer=NativeGraph::load(file,threads);new_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
        if(order=="new-first"){new_load();old_load();}else{old_load();new_load();}
        same(old,newer);
        std::cout<<"PASS complete graph equality\nold_seconds\t"<<old_s<<"\nnew_seconds\t"<<new_s
                 <<"\nthreads\t"<<threads<<"\nvertices\t"<<old.labels.size()<<"\npaths\t"<<old.paths.size()<<'\n';
        return 0;
    }
    require(argc==1,"Usage: test_graph_loader [graph.gfa[.gz] threads old-first|new-first]");
    Temp tmp;const auto p=tmp.dir/"graph.gfa.gz";
    std::mt19937_64 rng(9143);
    const std::vector<std::string> awkward={"01","1","0","4294967296","-2","+3","xlong-name","0000000000000","17x"};
    for(unsigned trial=0;trial<400;++trial) {
        const unsigned n=2+rng()%8,H=1+rng()%9;std::vector<std::string> names;
        for(unsigned i=0;i<n;++i)names.push_back(trial%3?std::to_string(i):awkward[i]);
        std::string labels,links,paths;
        for(unsigned i=0;i<n;++i)labels+="S\t"+names[i]+"\tacGTn\tXX:Z:ignored\n";
        // Small fully connected oriented graph, including repeat cycles.
        for(unsigned a=0;a<2*n;++a)for(unsigned b=0;b<2*n;++b)
            links+="L\t"+names[a/2]+"\t"+(a&1?"-":"+")+"\t"+names[b/2]+"\t"+(b&1?"-":"+")+"\t0M\n";
        for(unsigned h=0;h<H;++h) {
            const unsigned length=1+rng()%40;std::string steps;
            for(unsigned i=0;i<length;++i) {
                const unsigned v=rng()%(2*n);
                if(h&1)steps+=(v&1?"<":">")+names[v/2];
                else steps+=(i?",":"")+names[v/2]+(v&1?"-":"+");
            }
            paths+=(h&1)?"W\tsample\t"+std::to_string(h)+"\tchr1\t0\t*\t"+steps+"\n":
                           "P\thap"+std::to_string(h)+"\t"+steps+"\t*\n";
        }
        std::string text="H\tVN:Z:1.0\n# ignored\n\n";
        text+=(trial%3==0?paths+links+labels:trial%3==1?links+labels+paths:labels+links+paths);
        if(trial%2==0)text.pop_back();
        if(trial%5==0)for(size_t i=0;(i=text.find('\n',i))!=std::string::npos;i+=2)text.insert(i,1,'\r');
        check(p,text,trial%2);
    }
    const std::string base="S\t0\tA\nS\t1\tC\nL\t0\t+\t1\t+\t0M\n";
    for(const auto&bad:std::vector<std::string>{
        "", "S\t0\tA\n", "P\th\t0+\n", "S\t0\t*\nP\th\t0+\n",
        base+"S\t0\tA\nP\th\t0+\n", base+"P\th\t0+,\n",base+"P\th\t0\n",
        base+"P\th\t0+,1-\n",base+"W\ts\t1\tc\t0\t*\t>0<\n",
        base+"W\ts\t1\tc\t0\t*\t\n",base+"W\ts\t1\tc\t0\t*\t0\n",
        base+"L\t0\t?\t1\t+\t0M\nP\th\t0+\n",base+"L\t0\t+\t1\t+\t1M\nP\th\t0+\n",
        base+"P\tfirst-bad\t1+,0+\nP\tsecond-bad\t0-,1+\n"})
        for(bool gzip:{false,true})check(p,bad,gzip,8);
    std::string long_path="S\t0\tA\nL\t0\t+\t0\t+\t0M\n";
    for(unsigned h=0;h<8;++h) {
        long_path+="P\th"+std::to_string(h)+"\t0+";
        for(unsigned i=1;i<400000;++i)long_path+=",0+";
        long_path+="\t*\n";
    }
    check(p,long_path,true,64); // each P record exceeds TextInput's 1 MiB buffer
    require(!error(p.string(),0,false).empty(),"zero workers accepted");
    std::cout<<"PASS graph parser: numeric/arbitrary names, forward references, P/W, cycles, gzip/plain, CRLF, long lines, invalid records and thread parity\n";
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
