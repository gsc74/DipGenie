#include "VariantOutput.hpp"
#include "bgzf.h"
#include "hts.h"
#include <algorithm>
#include <array>
#include <cctype>
#include <climits>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unistd.h>
#include <zlib.h>

namespace dipcall_vcf::output {
namespace {
namespace fs=std::filesystem;
void need(bool ok,const std::string&message) {if(!ok)throw std::runtime_error("VCF output: "+message);}
bool compressed(const std::string&s){return s.size()>=3&&s.compare(s.size()-3,3,".gz")==0;}
std::vector<std::string> split(const std::string&s,char delimiter) {
    std::vector<std::string> v;size_t p=0;
    do {const auto e=s.find(delimiter,p);v.push_back(s.substr(p,e-p));if(e==std::string::npos)break;p=e+1;}while(true);
    return v;
}
uint64_t integer(const std::string&s) {
    need(!s.empty(),"empty integer");uint64_t n=0;
    for(char c:s){need(c>='0'&&c<='9'&&n<=(UINT64_MAX-(c-'0'))/10,"invalid/overflowing integer");n=n*10+(c-'0');}
    return n;
}
struct Input {
    gzFile f;
    explicit Input(const std::string&p):f(gzopen(p.c_str(),"rb")){need(f,"cannot open "+p);}
    ~Input(){gzclose(f);}
    bool line(std::string&s) {
        s.clear();char b[65536];
        while(gzgets(f,b,sizeof b)) {s+=b;if(s.back()=='\n'){s.pop_back();if(!s.empty()&&s.back()=='\r')s.pop_back();return true;}}
        int e=0;gzerror(f,&e);need(e==Z_OK||e==Z_STREAM_END,"truncated compressed input");return !s.empty();
    }
};
struct Header {
    std::vector<std::string> lines,names;
    std::map<std::string,std::string> definitions;
    std::map<std::string,size_t> ids;
    std::string samples;
    void add(const std::string&s) {
        need(samples.empty(),"header after sample line");
        if(s.rfind("#CHROM\t",0)==0){auto c=split(s,'\t');need(c.size()==10&&!c[9].empty(),"expected one sample");samples=s;return;}
        need(s.rfind("##",0)==0,"invalid VCF header");
        auto id=s.find("=<ID=");
        if(id!=std::string::npos) {
            auto end=s.find_first_of(",>",id+5);need(end!=std::string::npos,"malformed header definition");
            auto name=s.substr(id+5,end-id-5);need(!name.empty(),"empty header ID");
            auto key=s.substr(0,id)+":"+name;
            const auto [it,added]=definitions.emplace(key,s);need(added||it->second==s,"conflicting header definition "+key);
            if(!added)return;
            if(s.rfind("##contig=",0)==0){ids.emplace(name,names.size());names.push_back(name);}
        } else if(std::find(lines.begin(),lines.end(),s)!=lines.end())return;
        lines.push_back(s);
    }
    void merge(const Header&h) {
        need(samples==h.samples,"component sample/header mismatch");
        auto saved=std::move(samples);samples.clear();for(const auto&l:h.lines)add(l);samples=std::move(saved);
    }
    void valid()const{need(!samples.empty()&&!names.empty(),"missing sample or contig dictionary");}
};
struct Record {
    std::array<std::string,10> c;
    uint64_t pos=0;
    size_t rid=0;
    Record()=default;
    Record(const std::string&line,const Header&h) {
        auto v=split(line,'\t');need(v.size()==10,"expected ten VCF columns");
        std::move(v.begin(),v.end(),c.begin());pos=integer(c[1]);need(pos&&pos<=INT_MAX,"position outside supported reference range");
        const auto it=h.ids.find(c[0]);need(it!=h.ids.end(),"undeclared contig "+c[0]);rid=it->second;
        need(!c[3].empty()&&!c[4].empty(),"empty allele");
    }
    std::string text()const{std::string s;for(size_t i=0;i<c.size();++i){if(i)s+='\t';s+=c[i];}return s+'\n';}
};
bool less(const Record&a,const Record&b){return std::tie(a.rid,a.pos,a.c[3],a.c[4])<std::tie(b.rid,b.pos,b.c[3],b.c[4]);}
Header read_header(Input&in,std::string&first) {
    Header h;
    while(in.line(first)) {if(first.empty())continue;if(first[0]!='#')break;h.add(first);first.clear();}
    h.valid();return h;
}
std::map<std::string,std::string> read_reference(const std::string&file) {
    Input in(file);std::map<std::string,std::string> result;std::string l,name;
    while(in.line(l)) {
        if(l.empty())continue;
        if(l[0]=='>'){name=l.substr(1,l.find_first_of(" \t")-1);need(!name.empty()&&result.emplace(name,"").second,"duplicate/empty reference record");}
        else {need(!name.empty(),"FASTA sequence before header");for(unsigned char c:l){need(std::strchr("ACGTNRYWSKMBDHVacgtnrywskmbdhv",c),"invalid reference base");result[name]+=char(std::toupper(c));}}
    }
    need(!result.empty(),"empty reference");return result;
}
char reference_base(char c){return c=='A'||c=='C'||c=='G'||c=='T'?c:'N';}
// Left-normalize a literal biallelic edit: repeatedly remove a common suffix,
// extending with the preceding reference base when an allele becomes empty,
// then trim a common prefix while retaining one anchoring base. No atomization.
void left_normalize(Record&r,const std::string&reference) {
    auto&a=r.c[3];auto&b=r.c[4];size_t p=r.pos-1;
    need(p<=reference.size()&&a.size()<=reference.size()-p&&reference.compare(p,a.size(),a)==0,
         "REF mismatch at "+r.c[0]+":"+r.c[1]);
    // A spanning-deletion allele denotes sequence removed by an upstream
    // call, not a literal edit at this position. Preserve its REF and POS,
    // as bcftools norm does; still check REF against the reference above.
    if(b=="*")return;
    need(a!=b,"identical REF/ALT");
    while(a.back()==b.back()) {
        if((a.size()==1||b.size()==1)&&p==0)break;
        a.pop_back();b.pop_back();
        if(a.empty()||b.empty()) {
            // Pad in chunks so a long homopolymer does not cause quadratic
            // insertions for every shifted base of a long indel.
            const size_t n=std::min<size_t>(p,100);std::string pad=reference.substr(p-n,n);
            for(char&c:pad)c=reference_base(c);
            a.insert(0,pad);b.insert(0,pad);p-=n;
        }
    }
    size_t left=0;while(left+1<a.size()&&left+1<b.size()&&a[left]==b[left])++left;
    a.erase(0,left);b.erase(0,left);r.pos=p+left+1;r.c[1]=std::to_string(r.pos);
}
// Write records and CSI together using the already bundled HTS primitives.
// Construct the index explicitly: the old generic tbx helper assumes >=1
// record and dereferences null for an otherwise valid header-only VCF.
void write(const Header&h,const std::vector<Record>&rows,const std::string&dest) {
    need(!fs::exists(dest)&&!fs::exists(dest+".csi")&&!fs::exists(dest+".tbi"),"refusing to overwrite output/index");
    h.valid();const bool zip=compressed(dest);
    std::unique_ptr<BGZF,decltype(&bgzf_close)> out(nullptr,bgzf_close);
    std::ofstream plain;
    if(zip)out.reset(bgzf_open(dest.c_str(),"w"));else plain.open(dest);
    need(zip?bool(out):bool(plain),"cannot create "+dest);
    auto put=[&](const std::string&s){if(zip)need(bgzf_write(out.get(),s.data(),s.size())==ssize_t(s.size()),"BGZF write failed");else{plain<<s;need(bool(plain),"VCF write failed");}};
    for(const auto&s:h.lines)put(s+"\n");
    put(h.samples+"\n");
    std::unique_ptr<hts_idx_t,decltype(&hts_idx_destroy)> idx(nullptr,hts_idx_destroy);
    if(zip){idx.reset(hts_idx_init(int(h.names.size()),HTS_FMT_CSI,bgzf_tell(out.get()),14,6));need(bool(idx),"CSI allocation failed");}
    size_t previous_rid=0;uint64_t previous_pos=0;bool first=true;
    for(const auto&r:rows) {
        need(first||r.rid>previous_rid||(r.rid==previous_rid&&r.pos>=previous_pos),"records not sorted");
        need(r.rid<h.names.size()&&h.names[r.rid]==r.c[0],"record/header contig mismatch");
        need(r.c[3].size()<=size_t(INT_MAX)-r.pos+1,"variant exceeds supported coordinate range");
        put(r.text());
        if(zip)need(hts_idx_push(idx.get(),r.rid,r.pos-1,r.pos-1+r.c[3].size(),bgzf_tell(out.get()),1)==0,"CSI record rejected");
        previous_rid=r.rid;previous_pos=r.pos;first=false;
    }
    if(zip) {
        hts_idx_finish(idx.get(),bgzf_tell(out.get()));
        std::vector<uint8_t> meta;
        auto word=[&](uint32_t v){for(unsigned i=0;i<4;++i)meta.push_back(uint8_t(v>>(i*8)));};
        uint64_t name_bytes=0;for(const auto&s:h.names)name_bytes+=s.size()+1;
        need(name_bytes<=INT_MAX-28,"CSI names overflow");
        for(uint32_t n:{uint32_t(2),uint32_t(1),uint32_t(2),uint32_t(0),uint32_t('#'),uint32_t(0),uint32_t(name_bytes)})word(n);
        for(const auto&s:h.names){meta.insert(meta.end(),s.begin(),s.end());meta.push_back(0);}
        hts_idx_set_meta(idx.get(),meta.size(),meta.data(),1);
        need(bgzf_close(out.release())==0,"BGZF close failed");
        // Work files live in a private directory in production. Check that the
        // void upstream index writer produced a readable index before publish.
        hts_idx_save(idx.get(),dest.c_str(),HTS_FMT_CSI);
        std::unique_ptr<hts_idx_t,decltype(&hts_idx_destroy)> checked(hts_idx_load_direct((dest+".csi").c_str(),HTS_FMT_CSI),hts_idx_destroy);
        need(bool(checked),"CSI write failed");
    } else {plain.close();need(bool(plain),"VCF close failed");}
}
} // namespace

void normalize(const std::string&reference,const std::string&input,const std::string&dest) {
    const auto refs=read_reference(reference);Input in(input);std::string line;const auto h=read_header(in,line);
    std::vector<Record> rows;
    do {
        if(line.empty())continue;
        Record base(line,h);const auto filters=split(base.c[6],';');
        if(std::find(filters.begin(),filters.end(),".")==filters.end()&&std::find(filters.begin(),filters.end(),"PASS")==filters.end())continue;
        need(base.c[7]=="."&&base.c[8]=="GT:AD","internal caller contract requires INFO '.' and FORMAT GT:AD");
        for(auto j:{3,4})for(char&c:base.c[j])c=char(std::toupper(static_cast<unsigned char>(c)));
        const auto alts=split(base.c[4],',');auto values=split(base.c[9],':');need(values.size()==2,"malformed GT:AD");
        auto gt=split(values[0],'|'),ad=split(values[1],',');need(gt.size()==2&&ad.size()==alts.size()+1,"invalid diploid genotype or AD cardinality");
        for(const auto&v:gt)if(v!=".")need(integer(v)<=alts.size(),"GT allele outside ALT list");
        for(const auto&v:ad)integer(v);
        auto ref=refs.find(base.c[0]);need(ref!=refs.end(),"contig absent from reference");
        for(size_t a=0;a<alts.size();++a) {
            Record r=base;r.c[4]=alts[a];
            for(auto j:{3,4}) {
                if(j==4&&r.c[j]=="*")continue;
                for(char c:r.c[j])need(std::strchr("ACGTN",c),"non-literal or ambiguous caller allele at "+r.c[0]+":"+r.c[1]);
            }
            left_normalize(r,ref->second);
            std::array<std::string,2> g;
            for(size_t k=0;k<2;++k)g[k]=gt[k]=="."?".":(integer(gt[k])==a+1?"1":"0");
            r.c[9]=g[0]+"|"+g[1]+":"+ad[0]+","+ad[a+1];rows.push_back(std::move(r));
        }
    }while(in.line(line));
    std::stable_sort(rows.begin(),rows.end(),less);
    // First identical allele survives, even if its genotype is later filtered,
    // matching norm -d exact followed by view's called-nonreference predicate.
    rows.erase(std::unique(rows.begin(),rows.end(),[](const auto&a,const auto&b){return !less(a,b)&&!less(b,a);}),rows.end());
    rows.erase(std::remove_if(rows.begin(),rows.end(),[](const auto&r){const auto gt=r.c[9].substr(0,r.c[9].find(':'));return gt.find('.')!=std::string::npos||gt=="0|0";}),rows.end());
    write(h,rows,dest);
}
void publish(const std::string&input,const std::string&dest) {
    Input in(input);std::string line;const auto h=read_header(in,line);std::vector<Record> rows;
    do {if(!line.empty())rows.emplace_back(line,h);}while(in.line(line));write(h,rows,dest);
}
void combine(const std::vector<std::string>&inputs,const std::string&dest) {
    need(!inputs.empty(),"no component calls");Header header;std::vector<Record> rows;
    for(const auto&file:inputs) {
        Input in(file);std::string line;auto h=read_header(in,line);
        if(header.samples.empty())header=h;else header.merge(h);
        do {if(!line.empty())rows.emplace_back(line,header);}while(in.line(line));
    }
    std::stable_sort(rows.begin(),rows.end(),less);write(header,rows,dest);
}
} // namespace dipcall_vcf::output
