#pragma once
#include "RowExecutor.hpp"
#include <algorithm>
#include <charconv>
#include <cctype>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <zlib.h>

namespace hic {
// zlib transparently reads uncompressed GFA as well as gzip streams. Some
// existing PGGB files have a .gz suffix but contain plain text.
class TextInput {
    gzFile file=nullptr;
    std::vector<char> buffer;
    size_t begin=0,end=0;
public:
    explicit TextInput(const std::string& path):buffer(1<<20) {
        file=gzopen(path.c_str(),"rb");if(!file)throw std::runtime_error("Cannot open "+path);
    }
    TextInput(const TextInput&)=delete;
    ~TextInput(){if(file)gzclose(file);}
    bool getline(std::string& line) {
        line.clear();
        for(;;) {
            if(begin==end) {
                int n=gzread(file,buffer.data(),unsigned(buffer.size()));
                if(n<0){int error=0;const char*why=gzerror(file,&error);throw std::runtime_error(why);}
                begin=0;end=size_t(n);if(!n)return !line.empty();
            }
            auto p=std::find(buffer.begin()+begin,buffer.begin()+end,'\n');
            size_t stop=size_t(p-buffer.begin());line.append(buffer.data()+begin,stop-begin);begin=stop;
            if(begin<end){++begin;if(!line.empty()&&line.back()=='\r')line.pop_back();return true;}
        }
    }
};
inline std::vector<std::string> fields(const std::string& line,char separator='\t') {
    std::vector<std::string> result;size_t a=0;
    for(;;){size_t b=line.find(separator,a);result.push_back(line.substr(a,b-a));if(b==std::string::npos)break;a=b+1;}
    return result;
}
inline char complement(char c) {
    switch(std::toupper(static_cast<unsigned char>(c))){case 'A':return 'T';case 'C':return 'G';case 'G':return 'C';case 'T':return 'A';default:return 'N';}
}
inline std::string reverse_complement(std::string s){std::reverse(s.begin(),s.end());for(char&c:s)c=complement(c);return s;}

struct NativeGraph {
    struct Path {std::string name;std::vector<uint32_t> walk;};
    std::vector<std::string> names,labels;
    std::vector<Path> paths;
    std::unordered_set<uint64_t> edges;
    static uint64_t edge(uint32_t a,uint32_t b){return(uint64_t(a)<<32)|b;}
    std::string spell_vertex(uint32_t v) const {
        const auto&s=labels.at(v/2);return(v&1)?reverse_complement(s):s;
    }
    std::string spell(size_t h) const {
        size_t length=0;
        for(uint32_t v:paths.at(h).walk) {
            const auto n=labels.at(v/2).size();
            if(n>SIZE_MAX-length)throw std::overflow_error("Panel spelling length");
            length+=n;
        }
        std::string result;result.reserve(length);
        for(uint32_t v:paths[h].walk) {
            const auto&s=labels[v/2];
            if(!(v&1))result+=s;
            else for(auto i=s.rbegin();i!=s.rend();++i)result+=complement(*i);
        }
        return result;
    }
    static NativeGraph load(const std::string& file,int threads=1) {
        if(threads<1)throw std::invalid_argument("Positive graph-loader thread count required");
        NativeGraph g;std::unordered_map<std::string,uint32_t> ids;
        auto node=[&](std::string_view name) {
            // Saved shards use canonical decimal IDs in insertion order.
            // Check the spelling as well: arbitrary numeric names, leading
            // zeros and forward references must retain the general GFA rules.
            if(!name.empty()&&name.size()<=10) {
                uint32_t n=0;const auto parsed=std::from_chars(name.data(),name.data()+name.size(),n);
                if(parsed.ec==std::errc{}&&parsed.ptr==name.data()+name.size()&&n<g.names.size()
                   &&std::string_view(g.names[n])==name)return n;
            }
            auto [it,inserted]=ids.emplace(std::string(name),uint32_t(ids.size()));
            if(inserted){if(ids.size()>UINT32_MAX/2)throw std::runtime_error("Too many oriented vertices");g.names.emplace_back(name);g.labels.emplace_back();}
            return it->second;
        };
        TextInput in(file);std::string line;std::vector<std::string_view> f;f.reserve(8);
        while(in.getline(line)) {
            if(line.empty()||line[0]=='#')continue;
            // Views are consumed before the next read; never retain a view
            // into the reusable line buffer in graph storage.
            f.clear();size_t a=0;
            for(;;) {
                const size_t b=line.find('\t',a);
                f.emplace_back(line.data()+a,(b==std::string::npos?line.size():b)-a);
                if(b==std::string::npos)break;
                a=b+1;
            }
            if(f[0]=="S") {
                if(f.size()<3||f[2]=="*")throw std::runtime_error("GFA requires explicit segment sequences");
                auto n=node(f[1]);if(!g.labels[n].empty())throw std::runtime_error("Duplicate segment");
                g.labels[n]=f[2];for(char&c:g.labels[n])c=std::toupper(static_cast<unsigned char>(c));
            } else if(f[0]=="L") {
                if(f.size()<6||f[5]!="0M"||(f[2]!="+"&&f[2]!="-")||(f[4]!="+"&&f[4]!="-"))
                    throw std::runtime_error("GFA links must use explicit zero overlap and strand");
                uint32_t a=2*node(f[1])+(f[2]=="-"),b=2*node(f[3])+(f[4]=="-");
                g.edges.insert(edge(a,b));g.edges.insert(edge(b^1,a^1));
            } else if(f[0]=="P"||f[0]=="W") {
                Path p;
                if(f[0]=="P") {
                    if(f.size()<3)throw std::runtime_error("Truncated P record");
                    p.name=f[1];
                    const auto steps=f[2];p.walk.reserve(size_t(std::count(steps.begin(),steps.end(),','))+1);
                    size_t first=0;
                    for(;;) {
                        const size_t last=steps.find(',',first);
                        const auto s=steps.substr(first,last==std::string_view::npos?last:last-first);
                        if(s.size()<2||(s.back()!='+'&&s.back()!='-'))throw std::runtime_error("Invalid oriented path step");
                        p.walk.push_back(2*node(s.substr(0,s.size()-1))+(s.back()=='-'));
                        if(last==std::string_view::npos)break;
                        first=last+1;
                    }
                } else {
                    if(f.size()<7)throw std::runtime_error("Truncated W record");
                    p.name=std::string(f[1])+"#"+std::string(f[2])+"#"+std::string(f[3]);
                    size_t i=0;
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
        // All graph mutation has finished. Each worker validates one walk
        // against the same immutable edge set; report errors in input order.
        std::vector<unsigned char> invalid(g.paths.size());const auto&edges=g.edges;
        RowExecutor pool(int(std::min(size_t(threads),g.paths.size())));
        pool.for_each(g.paths.size(),[&](size_t h) {
            const auto&w=g.paths[h].walk;
            for(size_t i=1;i<w.size();++i)if(!edges.count(edge(w[i-1],w[i]))) {invalid[h]=1;break;}
        });
        for(size_t h=0;h<g.paths.size();++h)if(invalid[h])
            throw std::runtime_error("Panel step missing from GFA links: "+g.paths[h].name);
        return g;
    }
};

struct CommonSegments {
    std::vector<uint32_t> anchors;
    std::vector<std::vector<size_t>> cuts; // base coordinates on each panel string
    std::vector<std::string> sequences;
    size_t blocks()const{return cuts.at(0).size()-1;}
    // Every junction lies after a shared k-1 base context. Literal k-mers
    // crossing that junction therefore do not depend on the previous donor.
    static CommonSegments build(const NativeGraph&g,size_t k,int threads=1) {
        if(k<2||threads<1)throw std::invalid_argument("segment options");
        const size_t H=g.paths.size(),V=g.labels.size();
        if(!H||!V||V>SIZE_MAX/2)throw std::invalid_argument("Empty/oversized graph");
        RowExecutor pool(threads);
        const size_t packed=(2*V+63)/64;
        std::vector<std::vector<uint64_t>> present(H,std::vector<uint64_t>(packed));
        // Bounded per-haplotype counters, rather than worker-count copies of
        // chromosome-sized coordinate arrays.
        const size_t batch=std::max(size_t(1),std::min(size_t(threads),size_t(256u<<20)/V));
        for(size_t first=0;first<H;first+=batch)
            pool.for_each(std::min(batch,H-first),[&](size_t i) {
                const size_t h=first+i;std::vector<unsigned char> count(V);
                const auto&walk=g.paths[h].walk;auto&row=present[h];
                for(uint32_t v:walk){if(v/2>=V)throw std::invalid_argument("Path vertex");if(count[v/2]<2)++count[v/2];}
                for(uint32_t v:walk)if(count[v/2]==1&&g.labels[v/2].size()>=k-1)
                    row[v/64]|=uint64_t(1)<<(v%64);
            });
        std::vector<uint64_t> common(packed,UINT64_MAX);
        pool.for_each(packed,[&](size_t i){for(const auto&row:present)common[i]&=row[i];});
        std::vector<std::vector<uint64_t>>().swap(present);
        std::vector<uint32_t> candidates;
        std::vector<size_t> candidate_id(2*V,SIZE_MAX);
        for(uint32_t v:g.paths[0].walk)if(common[v/64]>>(v%64)&1) {
            candidate_id[v]=candidates.size();candidates.push_back(v);
        }
        std::vector<uint64_t>().swap(common);
        std::vector<std::vector<size_t>> positions(H,std::vector<size_t>(candidates.size(),SIZE_MAX));
        pool.for_each(H,[&](size_t h) {
            for(size_t i=0;i<g.paths[h].walk.size();++i) {
                const auto x=candidate_id[g.paths[h].walk[i]];
                if(x!=SIZE_MAX)positions[h][x]=i;
            }
        });
        std::vector<size_t>().swap(candidate_id);
        CommonSegments result;std::vector<size_t> previous(H,0);bool first=true;
        for(size_t x=0;x<candidates.size();++x) {
            bool valid=true;size_t terminal=0;
            for(size_t h=0;h<H;++h) {
                if(positions[h][x]==SIZE_MAX)throw std::logic_error("Missing common anchor");
                valid=valid&&(first||positions[h][x]>previous[h]);
                terminal+=positions[h][x]+1==g.paths[h].walk.size();
            }
            // A cut k-1 bases into an anchor of exactly k-1 bases ends a panel
            // haplotype on which the anchor is last. Such a boundary is not
            // strictly interior on every haplotype unless it ends all of them.
            if(terminal&&terminal<H&&g.labels[candidates[x]/2].size()==k-1)valid=false;
            if(valid){result.anchors.push_back(candidates[x]);for(size_t h=0;h<H;++h)previous[h]=positions[h][x];first=false;}
        }
        std::vector<std::vector<size_t>>().swap(positions);
        result.cuts.resize(H);result.sequences.resize(H);
        std::vector<uint64_t> selected(packed);
        for(uint32_t v:result.anchors)selected[v/64]|=uint64_t(1)<<(v%64);
        pool.for_each(H,[&](size_t h) {
            auto&cuts=result.cuts[h];cuts.reserve(result.anchors.size()+2);cuts.push_back(0);size_t offset=0;
            auto&sequence=result.sequences[h];sequence=g.spell(h);
            for(uint32_t v:g.paths[h].walk) {
                if(selected[v/64]>>(v%64)&1)cuts.push_back(offset+k-1);
                offset+=g.labels[v/2].size();
            }
            cuts.push_back(offset);
        });
        for(auto&cuts:result.cuts)if(cuts.size()!=result.anchors.size()+2||!std::is_sorted(cuts.begin(),cuts.end()))
            throw std::runtime_error("Inconsistent common anchors");
        return result;
    }
    std::string mosaic(const std::vector<size_t>& donors)const {
        if(donors.size()!=blocks())throw std::invalid_argument("Incomplete mosaic");
        std::string result;
        for(size_t b=0;b<donors.size();++b){size_t h=donors[b];result+=sequences.at(h).substr(cuts.at(h)[b],cuts[h][b+1]-cuts[h][b]);}
        return result;
    }
};
}
