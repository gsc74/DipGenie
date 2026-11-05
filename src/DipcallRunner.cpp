#include "DipcallRunner.hpp"
#include "AssemblyCalling.hpp"
#include "VariantOutput.hpp"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <unistd.h>
#include <zlib.h>

namespace dipcall_vcf {
namespace {
namespace fs = std::filesystem;
void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error("dipcall: " + message);
}
bool suffix(const std::string& s, const std::string& tail) {
    return s.size() >= tail.size() && s.compare(s.size()-tail.size(), tail.size(), tail)==0;
}
bool sample_ok(const std::string& s) {
    return !s.empty() && std::all_of(s.begin(),s.end(),[](unsigned char c) {
        return (c>='a'&&c<='z')||(c>='A'&&c<='Z')||(c>='0'&&c<='9')||c=='_'||c=='-';
    });
}
void save(const fs::path& file,const std::string& text) {
    std::ofstream out(file);out<<text;out.close();require(bool(out),"cannot write "+file.string());
}

struct Gzip {
    gzFile file;
    explicit Gzip(const fs::path& p):file(gzopen(p.c_str(),"rb")){require(file,"cannot read "+p.string());}
    ~Gzip(){if(file)gzclose(file);}
    bool line(std::string& out) {
        out.clear();char buffer[65536];
        while(gzgets(file,buffer,sizeof buffer)) {
            out+=buffer;
            if(!out.empty()&&out.back()=='\n'){out.pop_back();if(!out.empty()&&out.back()=='\r')out.pop_back();return true;}
        }
        int error=Z_OK;gzerror(file,&error);require(error==Z_OK||error==Z_STREAM_END,"truncated compressed input");
        return !out.empty();
    }
};

// Each invocation calls one component. Keep the reference name/coordinates;
// dipcall -a also accepts numeric or non-chromosome reference names.
std::string reference(const fs::path& source,const fs::path& dest,std::string& internal) {
    Gzip in(source);std::ofstream out(dest);std::string line,name;uint64_t bases=0;
    while(in.line(line)) {
        if(line.empty())continue;
        if(line[0]=='>') {
            require(name.empty(),"expected one reference record per graph component");
            name=line.substr(1,line.find_first_of(" \t")-1);
            require(!name.empty()&&name.find_first_of("\t\r\n,<>")==std::string::npos,"invalid reference name");
            internal=name;
            out<<'>'<<internal<<'\n';
        } else {
            require(!name.empty(),"reference sequence precedes its header");
            for(unsigned char c:line)require(c&&std::strchr("ACGTNRYWSKMBDHVacgtnrywskmbdhv",c),"invalid reference base");
            bases+=line.size();out<<line<<'\n';
        }
    }
    out.close();require(bool(out)&&bases>0,"empty or unwritable reference");return name;
}
void haplotype(const fs::path& file,const std::string& sequence,unsigned lane) {
    require(!sequence.empty(),"empty inferred haplotype");
    for(unsigned char c:sequence)require(c&&std::strchr("ACGTNRYWSKMBDHVacgtnrywskmbdhv",c),"invalid haplotype base");
    std::ofstream out(file);out<<">hap"<<lane<<'\n';
    for(size_t i=0;i<sequence.size();i+=80)out.write(sequence.data()+i,std::min(size_t(80),sequence.size()-i)).put('\n');
    out.close();require(bool(out),"cannot write haplotype input");
}
std::vector<std::string> fields(const std::string& value,char delimiter) {
    std::vector<std::string> result;size_t begin=0;
    while(true) {
        const auto end=value.find(delimiter,begin);
        result.push_back(value.substr(begin,end-begin));
        if(end==std::string::npos)return result;
        begin=end+1;
    }
}
void validate_sample(const fs::path& vcf,const std::string& sample,bool phased=false,bool allow_unphased=false) {
    Gzip in(vcf);std::string line;bool header=false;
    while(in.line(line)) {
        if(line.rfind("#CHROM\t",0)==0) {
            require(!header,"duplicate VCF sample header");
            require(std::count(line.begin(),line.end(),'\t')==9&&suffix(line,"\t"+sample),"unexpected VCF sample or ploidy columns");
            header=true;
        }else if(!line.empty()&&line[0]!='#') {
            require(header,"VCF records precede sample header");
            const auto columns=fields(line,'\t');
            require(columns.size()==10,"unexpected VCF record column count");
            if(!phased)continue;
            const auto format=fields(columns[8],':'), values=fields(columns[9],':');
            const auto gt=std::find(format.begin(),format.end(),"GT");
            require(gt!=format.end(),"VCF record lacks GT");
            const auto index=size_t(gt-format.begin());
            require(index<values.size(),"VCF record lacks sample genotype");
            auto alleles=fields(values[index],'|');
            if(alleles.size()==1&&allow_unphased)alleles=fields(values[index],'/');
            require(alleles.size()==2,"expected phased diploid GT");
            const auto maximum=fields(columns[4],',').size();
            for(const auto& allele:alleles) {
                require(!allele.empty(),"empty genotype allele");size_t number=0;
                for(unsigned char c:allele) {
                    require(c>='0'&&c<='9',"noncalled or unphased genotype");
                    require(number<=maximum/10,"genotype allele out of range");
                    number=10*number+(c-'0');
                    require(number<=maximum,"genotype allele out of range");
                }
            }
        }
    }
    require(header,"caller did not produce a VCF sample header");
}
} // namespace

size_t assign_phase_sets(const std::string& source,const std::string& dest,
                         const std::vector<std::pair<uint64_t,uint64_t>>& intervals) {
    require(std::is_sorted(intervals.begin(),intervals.end()),"phase-set intervals must be sorted");
    for(size_t i=1;i<intervals.size();++i)require(intervals[i-1].second<=intervals[i].first,"phase-set intervals overlap");
    Gzip in(source);std::ofstream out(dest);std::string line;size_t separate=0;
    while(in.line(line)) {
        if(line.rfind("#CHROM\t",0)==0)out<<"##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">\n";
        if(line.empty()||line[0]=='#'){out<<line<<'\n';continue;}
        auto columns=fields(line,'\t');require(columns.size()==10,"invalid VCF record");
        auto format=fields(columns[8],':'),values=fields(columns[9],':');
        require(!format.empty()&&format[0]=="GT"&&std::find(format.begin(),format.end(),"PS")==format.end(),"unexpected VCF FORMAT");
        const auto alleles=fields(values[0],'|');
        if(alleles.size()==2&&alleles[0]!=alleles[1]) {
            const uint64_t pos=std::stoull(columns[1])-1;
            auto it=std::upper_bound(intervals.begin(),intervals.end(),std::make_pair(pos,UINT64_MAX));
            uint64_t ps=1;
            if(it!=intervals.begin()&&pos<std::prev(it)->second){ps=std::prev(it)->first+2;++separate;}
            columns[8]+=":PS";columns[9]+=":"+std::to_string(ps);
        }
        for(size_t i=0;i<columns.size();++i)out<<(i?"\t":"")<<columns[i];
        out<<'\n';
    }
    out.close();require(bool(out),"cannot write phase-set VCF");
    return separate;
}

size_t assign_block_phase_sets(const std::string& source,const std::string& dest,const std::vector<Config::Block>& blocks) {
    for(size_t i=0;i<blocks.size();++i) {
        require(blocks[i].begin<blocks[i].end&&blocks[i].ps>0&&blocks[i].ps<1000000000ULL,"invalid phase block");
        require(!i||blocks[i-1].end<=blocks[i].begin,"phase blocks must be sorted and disjoint");
    }
    Gzip in(source);std::ofstream out(dest);std::string line;size_t outside=0,index=0;
    while(in.line(line)) {
        if(line.rfind("#CHROM\t",0)==0)out<<"##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">\n";
        if(line.empty()||line[0]=='#'){out<<line<<'\n';continue;}
        auto columns=fields(line,'\t');require(columns.size()==10,"invalid VCF record");
        auto format=fields(columns[8],':'),values=fields(columns[9],':');
        require(!format.empty()&&format[0]=="GT"&&std::find(format.begin(),format.end(),"PS")==format.end(),"unexpected VCF FORMAT");
        const auto alleles=fields(values[0],'|');
        if(alleles.size()==2&&alleles[0]!=alleles[1]) {
            const uint64_t pos=std::stoull(columns[1])-1;
            auto it=std::upper_bound(blocks.begin(),blocks.end(),pos,[](uint64_t p,const Config::Block&b){return p<b.begin;});
            uint64_t ps;
            if(it!=blocks.begin()&&pos<std::prev(it)->end)ps=std::prev(it)->ps;
            else {ps=2000000000ULL+index;++outside;}
            columns[8]+=":PS";columns[9]+=":"+std::to_string(ps);
        }
        ++index;
        for(size_t i=0;i<columns.size();++i)out<<(i?"\t":"")<<columns[i];
        out<<'\n';
    }
    out.close();require(bool(out),"cannot write block phase-set VCF");
    return outside;
}

void local_haplotypes(const std::string& reference,const std::vector<SiteRecord>& records,size_t i,size_t flank,
                      std::array<std::string,2>& haps,std::array<std::pair<size_t,size_t>,2>& span) {
    const auto& centre=records.at(i);
    const uint64_t lo=centre.pos>flank?centre.pos-flank:0,hi=std::min<uint64_t>(reference.size(),centre.pos+centre.ref.size()+flank);
    size_t first=i;while(first>0&&records[first-1].pos+64>=lo)--first;
    for(int h=0;h<2;++h) {
        std::string out;uint64_t cur=lo;span[h]={0,0};
        for(size_t j=first;j<records.size()&&records[j].pos<hi;++j) {
            const auto& r=records[j];
            if(r.pos<cur||r.pos<lo)continue;
            const std::string& a=r.allele[h];
            if(a==r.ref&&j!=i)continue;
            out+=reference.substr(cur,r.pos-cur);
            if(j==i)span[h]={out.size(),out.size()+a.size()};
            out+=a;cur=r.pos+r.ref.size();
        }
        if(cur<hi)out+=reference.substr(cur,hi-cur);
        haps[h]=std::move(out);
    }
}

size_t assign_site_phase_sets(const std::string& source,const std::string& dest,const std::string& reference,
    const std::function<bool(const std::array<std::string,2>&,const std::array<std::pair<size_t,size_t>,2>&)>& supported,
    bool any_block) {
    Gzip in(source);std::string line;std::vector<std::string> header;std::vector<std::vector<std::string>> body;
    while(in.line(line)){if(line.empty())continue;if(line[0]=='#')header.push_back(line);else body.push_back(fields(line,'\t'));}
    std::vector<SiteRecord> records;records.reserve(body.size());
    for(const auto& c:body) {
        require(c.size()==10,"invalid VCF record");
        auto format=fields(c[8],':'),values=fields(c[9],':');
        require(!format.empty()&&format[0]=="GT","GT must be first");
        auto gt=fields(values[0],'|');if(gt.size()!=2)gt=fields(values[0],'/');
        require(gt.size()==2,"diploid GT required");
        auto alts=fields(c[4],',');SiteRecord r;r.pos=std::stoull(c[1])-1;r.ref=c[3];
        for(int h=0;h<2;++h){const auto k=gt[h]=="."?0:std::stoul(gt[h]);r.allele[h]=k?alts.at(k-1):r.ref;}
        const auto ps=std::find(format.begin(),format.end(),"PS");
        const std::string id=ps!=format.end()?values[size_t(ps-format.begin())]:".";
        r.linked=r.allele[0]!=r.allele[1]&&(any_block?(id!="."&&std::stoull(id)<1000000000ULL):id=="1");
        records.push_back(std::move(r));
    }
    for(size_t i=1;i<records.size();++i)require(records[i-1].pos<=records[i].pos,"VCF must be sorted");
    std::ofstream out(dest);size_t changed=0;
    for(const auto& h:header)out<<h<<'\n';
    for(size_t i=0;i<body.size();++i) {
        auto c=body[i];
        if(records[i].linked) {
            std::array<std::string,2> haps;std::array<std::pair<size_t,size_t>,2> span;
            local_haplotypes(reference,records,i,48,haps,span);
            if(!supported(haps,span)) {
                auto format=fields(c[8],':'),values=fields(c[9],':');
                values[size_t(std::find(format.begin(),format.end(),"PS")-format.begin())]=std::to_string(1000000000ULL+i);
                std::string joined;for(size_t k=0;k<values.size();++k)joined+=(k?":":"")+values[k];
                c[9]=joined;++changed;
            }
        }
        for(size_t k=0;k<c.size();++k)out<<(k?"\t":"")<<c[k];
        out<<'\n';
    }
    out.close();require(bool(out),"cannot write site phase-set VCF");
    return changed;
}

std::string run(const Config& cfg,const std::vector<std::string>& haps) {
    require(haps.size()==2,"exactly two haplotypes are required");
    require(sample_ok(cfg.sample),"sample must contain only letters, digits, '_' or '-'");
    require(cfg.threads>0,"thread count must be positive");
    require(suffix(cfg.output_vcf,".vcf")||suffix(cfg.output_vcf,".vcf.gz"),"output must end in .vcf or .vcf.gz");
    auto output=fs::absolute(cfg.output_vcf);fs::create_directories(output.parent_path());
    require(!fs::exists(output)&&!fs::exists(output.string()+".csi")&&!fs::exists(output.string()+".tbi"),
            "refusing to overwrite an existing VCF or index");
    require(output!=fs::absolute(cfg.reference_fasta),"output cannot overwrite reference");
    auto pattern=(output.parent_path()/".dipgenie-dipcall-XXXXXX").string();
    std::vector<char> buffer(pattern.begin(),pattern.end());buffer.push_back(0);
    char* created=mkdtemp(buffer.data());require(created,"cannot create private caller directory");
    const fs::path work(created);
    fprintf(stderr,"[M::dipcall] work=%s threads=%u sample=%s\n",work.c_str(),cfg.threads,cfg.sample.c_str());
    std::string internal_reference;
    reference(cfg.reference_fasta,work/"reference.fa",internal_reference);
    haplotype(work/"hap1.fa",haps[0],1);haplotype(work/"hap2.fa",haps[1],2);
    const auto prefix=(work/"calls").string();
    const auto stats=embedded::align((work/"reference.fa").string(),
        {(work/"hap1.fa").string(),(work/"hap2.fa").string()},prefix,cfg.threads,cfg.cpu_ids);
    embedded::pileup((work/"reference.fa").string(),{prefix+".hap1.bam",prefix+".hap2.bam"},prefix+".pair.vcf");
    embedded::pair_vcf(prefix+".pair.vcf",prefix+".dip.vcf",cfg.sample);
    embedded::callable_bed({prefix+".hap1.paf",prefix+".hap2.paf"},prefix);
    save(work/"backend.tsv","backend\tversion\nminimap2\t2.20-r1061\nhtsbox\t3ca3e06\ndipcall\t0.3-all-contigs-diploid\n");
    fprintf(stderr,"[M::dipcall] retained alignments=%llu,%llu; calling complete\n",
        (unsigned long long)stats[0].retained_alignments,(unsigned long long)stats[1].retained_alignments);
    require(fs::is_regular_file(work/"calls.dip.bed"),"missing callable BED");
    validate_sample(work/"calls.dip.vcf",cfg.sample);
    output::normalize((work/"reference.fa").string(),(work/"calls.dip.vcf").string(),(work/"sorted.vcf").string());
    const bool compressed=suffix(cfg.output_vcf,".gz");
    const std::string result=compressed?"final.vcf.gz":"final.vcf";
    if(cfg.separate.empty()&&cfg.blocks.empty())output::publish((work/"sorted.vcf").string(),(work/result).string());
    else {
        if(!cfg.blocks.empty()) {
            const auto outside=assign_block_phase_sets((work/"sorted.vcf").string(),(work/"phase-sets.vcf").string(),cfg.blocks);
            fprintf(stderr,"[M::dipcall] read-linked phase blocks=%zu; heterozygous calls outside blocks=%zu\n",cfg.blocks.size(),outside);
        } else {
            const auto separate=assign_phase_sets((work/"sorted.vcf").string(),(work/"phase-sets.vcf").string(),cfg.separate);
            fprintf(stderr,"[M::dipcall] heterozygous calls in separate phase sets=%zu\n",separate);
        }
        std::string last="phase-sets.vcf";
        if(cfg.site_supported) {
            std::string sequence;{Gzip ref(work/"reference.fa");std::string l;while(ref.line(l))if(!l.empty()&&l[0]!='>')sequence+=l;}
            for(char&ch:sequence)ch=char(std::toupper(static_cast<unsigned char>(ch)));
            const auto unsupported=assign_site_phase_sets((work/last).string(),(work/"site-sets.vcf").string(),sequence,cfg.site_supported,!cfg.blocks.empty());
            fprintf(stderr,"[M::dipcall] linked-block heterozygous calls without phase support (own phase set)=%zu\n",unsupported);
            last="site-sets.vcf";
        }
        output::publish((work/last).string(),(work/result).string());
    }
    validate_sample(work/result,cfg.sample,true);
    // Work is a private sibling of the destination: rename is atomic. Keep
    // logs, raw dipcall calls, callable BED and haplotype order on success.
    fs::rename(work/result,output);
    std::error_code ec;
    fs::remove(output.string()+".tbi",ec);fs::remove(output.string()+".csi",ec);
    if(compressed)fs::rename(work/(result+".csi"),output.string()+".csi");
    return output.string();
}
} // namespace dipcall_vcf
