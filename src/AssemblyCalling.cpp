// C++ adapters for the pinned minimap2/HTSbox libraries and dipcall 0.3.
// The vcfpair/samflt policies derive from lh3/dipcall commit
// 6bd5d7724699491f215aeb5fb628490ebf2cc3ae (MIT, 2019- DFCI).
// See third_party/DIPCALL-LICENSE.txt; upstream algorithms keep their licenses.
#include "AssemblyCalling.hpp"
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>
#include <sched.h>
#include "mmpriv.h"
#include "sam.h"
#include "../third_party/htsbox/faidx.h"
extern "C" int dg_embedded_pileup(int,char**,FILE*);

namespace dipcall_vcf::embedded {
namespace {
namespace fs=std::filesystem;
void need(bool b,const std::string&m){if(!b)throw std::runtime_error("assembly caller: "+m);}
using Bam=std::unique_ptr<bam1_t,decltype(&bam_destroy1)>;
using Header=std::unique_ptr<bam_hdr_t,decltype(&bam_hdr_destroy)>;
struct Text {kstring_t v{0,0,nullptr};~Text(){std::free(v.s);}};
struct Regions {
    int n=0;mm_reg1_t *v=nullptr;
    ~Regions(){for(int i=0;i<n;++i)std::free(v[i].p);std::free(v);}
};
struct Queries {
    int n=0;mm_bseq1_t *v=nullptr;
    ~Queries(){for(int i=0;i<n;++i){std::free(v[i].name);std::free(v[i].seq);std::free(v[i].qual);std::free(v[i].comment);}std::free(v);}
};
struct Bgzf {
    BGZF*p;
    explicit Bgzf(const std::string&path,const char*mode):p(bgzf_open(path.c_str(),mode)){need(p,"cannot open "+path);}
    ~Bgzf(){if(p)bgzf_close(p);}
    void close(){auto*q=p;p=nullptr;need(bgzf_close(q)==0,"BAM close failed");}
};
bool sam_filter(const bam1_t*b) {
    if((b->core.flag&BAM_FSECONDARY)||b->core.qual<5)return false;
    uint64_t length=0;const auto*cigar=bam_get_cigar(b);
    for(uint32_t i=0;i<b->core.n_cigar;++i) {
        const auto op=bam_cigar_op(cigar[i]);
        if(op==BAM_CMATCH||op==BAM_CINS||op==BAM_CDEL)length+=bam_cigar_oplen(cigar[i]);
    }
    return length>=50000;
}
std::vector<std::string> split(const std::string&s,char delim) {
    std::vector<std::string> out;size_t pos=0;
    do{auto end=s.find(delim,pos);out.push_back(s.substr(pos,end-pos));
       if(end==std::string::npos)break;
       pos=end+1;}while(true);
    return out;
}
uint64_t number(const std::string&s) {
    need(!s.empty(),"empty integer");uint64_t n=0;
    for(unsigned char c:s){need(c>='0'&&c<='9',"invalid integer");need(n<=(UINT64_MAX-(c-'0'))/10,"integer overflow");n=n*10+c-'0';}
    return n;
}
void absent(const std::string&path){need(!fs::exists(path),"refusing to overwrite "+path);}
struct Affinity {
    cpu_set_t saved{};bool active=false;
    explicit Affinity(const std::vector<uint32_t>&cpus) {
        if(cpus.empty())return;
        cpu_set_t desired;CPU_ZERO(&desired);
        for(auto cpu:cpus){need(cpu<CPU_SETSIZE,"CPU id outside affinity range");CPU_SET(cpu,&desired);}
        need(sched_getaffinity(0,sizeof(saved),&saved)==0,"cannot read CPU affinity");
        need(sched_setaffinity(0,sizeof(desired),&desired)==0,"cannot set CPU affinity");active=true;
    }
    ~Affinity(){if(active)sched_setaffinity(0,sizeof(saved),&saved);}
};
}

std::array<AlignmentStats,2> align(const std::string&ref,
    const std::array<std::string,2>&fasta,const std::string&prefix,unsigned threads,const std::vector<uint32_t>&cpus) {
    need(threads>0&&threads<=INT_MAX,"invalid thread count");
    Affinity main_affinity(cpus);
    need(fs::is_regular_file(ref),"reference not found");
    for(unsigned h=0;h<2;++h) {
        need(fs::is_regular_file(fasta[h]),"haplotype FASTA not found");
        for(const auto*ext:{".bam",".paf"})absent(prefix+".hap"+std::to_string(h+1)+ext);
    }
    mm_idxopt_t io{};mm_mapopt_t mo{};mm_set_opt(nullptr,&io,&mo);
    need(mm_set_opt("asm5",&io,&mo)==0,"asm5 unavailable");
    mo.flag|=MM_F_CIGAR|MM_F_OUT_SAM|MM_F_OUT_CS|MM_F_OUT_CG|MM_F_PAF_NO_HIT;
    std::unique_ptr<mm_idx_reader_t,decltype(&mm_idx_reader_close)> reader(mm_idx_reader_open(ref.c_str(),&io,nullptr),mm_idx_reader_close);
    need(bool(reader),"cannot open assembly reference");
    std::unique_ptr<mm_idx_t,decltype(&mm_idx_destroy)> index(mm_idx_reader_read(reader.get(),int(threads)),mm_idx_destroy);
    need(index&&index->n_seq,"empty assembly reference");
    // The production caller supplies one component, never a split WGS index.
    // Reject multipart input rather than silently compute wrong MAPQ/calls.
    std::unique_ptr<mm_idx_t,decltype(&mm_idx_destroy)> extra(mm_idx_reader_read(reader.get(),1),mm_idx_destroy);
    need(!extra,"multipart reference index is unsupported for one-component calling");
    mm_mapopt_update(&mo,index.get());
    std::ostringstream header_text;header_text<<"@HD\tVN:1.6\tSO:coordinate\n";
    for(uint32_t i=0;i<index->n_seq;++i)header_text<<"@SQ\tSN:"<<index->seq[i].name<<"\tLN:"<<index->seq[i].len<<'\n';
    const std::string text=header_text.str();
    std::array<AlignmentStats,2> stats{};std::exception_ptr error;
    // Mapping the two complete haplotypes is independent. Splitting them into
    // artificial windows would change dipcall alignments and is not performed.
    #pragma omp parallel for num_threads(std::min(threads,2u)) schedule(static)
    for(int h=0;h<2;++h)try {
        Affinity worker_affinity(cpus);
        Header header(sam_hdr_parse(int(text.size()),text.c_str()),bam_hdr_destroy);
        need(bool(header),"invalid BAM dictionary");
        std::unique_ptr<mm_bseq_file_t,decltype(&mm_bseq_close)> input(mm_bseq_open(fasta[h].c_str()),mm_bseq_close);
        need(bool(input),"cannot read haplotype");
        std::unique_ptr<mm_tbuf_t,decltype(&mm_tbuf_destroy)> buffer(mm_tbuf_init(),mm_tbuf_destroy);
        need(bool(buffer),"mapping buffer allocation failed");
        std::vector<Bam> records;std::set<std::string> names;
        std::ofstream paf(prefix+".hap"+std::to_string(h+1)+".paf");need(bool(paf),"cannot write PAF");
        while(!mm_bseq_eof(input.get())) {
            Queries queries;queries.v=mm_bseq_read(input.get(),mo.mini_batch_size,0,&queries.n);
            for(int q=0;q<queries.n;++q) {
                auto&t=queries.v[q];need(t.l_seq>0&&names.insert(t.name).second,"empty/duplicate assembly record");
                ++stats[h].queries;Regions regs;
                regs.v=mm_map(index.get(),t.l_seq,t.seq,&regs.n,buffer.get(),&mo,t.name);
                stats[h].alignments+=regs.n;Text line;
                if(!regs.n){mm_write_paf(&line.v,index.get(),&t,nullptr,nullptr,mo.flag);paf<<line.v.s<<'\n';}
                const mm_reg1_t *all=regs.v;
                for(int j=0;j<regs.n;++j) {
                    mm_write_paf(&line.v,index.get(),&t,&regs.v[j],mm_tbuf_get_km(buffer.get()),mo.flag);
                    paf<<line.v.s<<'\n';
                    mm_write_sam2(&line.v,index.get(),&t,0,j,1,&regs.n,&all,mm_tbuf_get_km(buffer.get()),mo.flag);
                    Bam b(bam_init1(),bam_destroy1);need(bool(b),"BAM allocation failed");
                    need(sam_parse1(&line.v,header.get(),b.get())>=0,"cannot parse generated SAM");
                    if(sam_filter(b.get())){records.push_back(std::move(b));++stats[h].retained_alignments;}
                }
            }
        }
        need(stats[h].queries>0,"empty assembly FASTA");
        paf.close();need(bool(paf),"PAF write failed");
        // Match samtools coordinate sort's position/strand order; stable ties
        // retain query and alignment order without any thread-dependent rank.
        std::stable_sort(records.begin(),records.end(),[](const Bam&a,const Bam&b){
            return std::make_tuple(a->core.tid,a->core.pos,bool(a->core.flag&BAM_FREVERSE))
                 < std::make_tuple(b->core.tid,b->core.pos,bool(b->core.flag&BAM_FREVERSE));});
        Bgzf out(prefix+".hap"+std::to_string(h+1)+".bam","w");
        need(bam_hdr_write(out.p,header.get())>=0,"BAM header write failed");
        for(const auto&b:records)need(bam_write1(out.p,b.get())>=0,"BAM write failed");
        out.close();
    }catch(...){
        #pragma omp critical(dipgenie_assembly_error)
        {if(!error)error=std::current_exception();}
    }
    if(error)std::rethrow_exception(error);
    return stats;
}

void pileup(const std::string&reference,const std::array<std::string,2>&bam,const std::string&output) {
    absent(output);need(fs::is_regular_file(reference),"missing reference");
    std::vector<std::pair<std::string,uint32_t>> dictionary;
    for(size_t h=0;h<bam.size();++h) {
        Bgzf in(bam[h],"r");Header header(bam_hdr_read(in.p),bam_hdr_destroy);need(bool(header),"invalid BAM header");
        std::vector<std::pair<std::string,uint32_t>> current;
        for(int i=0;i<header->n_targets;++i)current.emplace_back(header->target_name[i],header->target_len[i]);
        if(!h)dictionary=current;else need(current==dictionary,"BAM dictionaries differ");
    }
    need(!dictionary.empty(),"empty BAM dictionary");
    static std::mutex upstream_mutex;std::lock_guard<std::mutex> lock(upstream_mutex);
    // Index only the caller's private reference, never a shared benchmark file.
    need(fs::is_regular_file(reference+".fai")||fai_build(reference.c_str())==0,"FASTA indexing failed");
    std::unique_ptr<faidx_t,decltype(&fai_destroy)> fai(fai_load(reference.c_str()),fai_destroy);
    need(bool(fai),"invalid FASTA index");
    for(const auto&[name,length]:dictionary)need(faidx_seq_len(fai.get(),name.c_str())==int64_t(length),"BAM/reference dictionary mismatch");
    std::vector<std::string> args={"pileup","-q5","-evcf",reference,bam[0],bam[1]};std::vector<char*> argv;
    for(auto&a:args)argv.push_back(a.data());
    argv.push_back(nullptr);
    FILE*out=std::fopen(output.c_str(),"wx");need(out,"cannot create pair VCF");
    const int status=dg_embedded_pileup(int(args.size()),argv.data(),out);
    const int closed=std::fclose(out);need(status==0&&closed==0,"embedded pileup failed");
}

void pair_vcf(const std::string&input,const std::string&output,const std::string&sample) {
    need(!sample.empty()&&sample.find_first_of("\t\r\n ")==std::string::npos,"invalid sample");absent(output);
    std::ifstream in(input);need(bool(in),"cannot read pair VCF");std::ofstream out(output);bool header=false;
    std::string line;
    while(std::getline(in,line)) {
        if(!line.empty()&&line.back()=='\r')line.pop_back();
        if(line.rfind("##source=",0)==0||line.rfind("##reference=",0)==0)continue;
        if(line.rfind("#CHROM\t",0)==0) {
            need(!header,"duplicate VCF header");auto c=split(line,'\t');need(c.size()==11,"pair VCF requires two samples");
            for(const auto*id:{"HET1","HET2","GAP1","GAP2"}) {
                const bool het=id[0]=='H';out<<"##FILTER=<ID="<<id<<",Description=\""<<(het?"Heterozygous":"Uncalled")
                    <<" in the "<<(id[3]=='1'?"first":"second")<<" haplotype\">\n";
            }
            c.resize(10);c[9]=sample;for(size_t i=0;i<c.size();++i)out<<(i?"\t":"")<<c[i];out<<'\n';header=true;continue;
        }
        if(!line.empty()&&line[0]=='#'){out<<line<<'\n';continue;}
        need(header&&!line.empty(),"VCF data/header order");auto c=split(line,'\t');need(c.size()==11&&c[8]=="GT:AD","expected pileup GT:AD");
        if(c[4].find('N')!=std::string::npos)continue;
        std::array<std::string,2> ht;std::vector<uint64_t> depths;std::string filters;
        const size_t n_alleles=split(c[4],',').size()+1;
        for(size_t h=0;h<2;++h) {
            auto sample_fields=split(c[9+h],':');need(sample_fields.size()==2,"invalid pileup sample");
            auto gt=split(sample_fields[0],'/'),ad=split(sample_fields[1],',');
            need(gt.size()==2&&ad.size()==n_alleles,"invalid pileup genotype/depths");
            if(!h)depths.resize(ad.size(),0);
            for(size_t j=0;j<ad.size();++j){const auto d=number(ad[j]);need(depths[j]<=UINT64_MAX-d,"AD overflow");depths[j]+=d;}
            for(const auto&a:gt)if(a!=".")need(number(a)<n_alleles,"GT outside allele list");
            std::string filter;
            if(gt[0]==".")filter="GAP";else if(gt[0]!=gt[1])filter="HET";
            ht[h]=filter.empty()?gt[0]:".";
            if(!filter.empty()){if(!filters.empty())filters+=';';filters+=filter+std::to_string(h+1);}
        }
        c.resize(10);c[5]="30";c[6]=filters.empty()?".":filters;c[9]=ht[0]+"|"+ht[1]+":";
        for(size_t j=0;j<depths.size();++j)c[9]+=(j?",":"")+std::to_string(depths[j]);
        for(size_t i=0;i<c.size();++i)out<<(i?"\t":"")<<c[i];
        out<<'\n';
    }
    need(in.eof()&&header,"incomplete pair VCF");out.close();need(bool(out),"diploid VCF write failed");
}

void callable_bed(const std::array<std::string,2>&paf,const std::string&prefix) {
    // Region-only adaptation of minimap2 paftools.js call (pinned 2.20), not
    // its variant caller. The VCF is always derived from HTSbox pileup above.
    struct Interval {std::string contig;uint64_t start,end;};
    struct Alignment {std::string contig,line;uint64_t start,end,length;};
    auto write=[](const std::string&file,const std::vector<Interval>&regions){
        absent(file);std::ofstream out(file);
        for(const auto&r:regions)out<<r.contig<<'\t'<<r.start<<'\t'<<r.end<<'\n';
        out.close();need(bool(out),"BED write failed");};
    std::array<std::vector<Interval>,2> regions;
    for(size_t h=0;h<2;++h) {
        std::ifstream in(paf[h]);need(bool(in),"cannot read PAF");std::string line;std::vector<Alignment> alignments;
        while(std::getline(in,line)) {
            auto c=split(line,'\t');need(c.size()>=12,"incomplete PAF record");
            if(c[5]=="*")continue;
            const auto length=number(c[10]),mapq=number(c[11]);
            if(length<10000||mapq<5)continue;
            bool s1=false,s2=false,secondary=false;
            for(size_t i=12;i<c.size();++i){s1|=c[i].rfind("s1:i:",0)==0;s2|=c[i].rfind("s2:i:",0)==0;
                secondary|=c[i]=="tp:A:S"||c[i]=="tp:A:i";}
            if((s1&&!s2)||secondary)continue;
            const auto start=number(c[7]),end=number(c[8]);
            need(start<end&&end<=number(c[6]),"invalid PAF target interval");
            alignments.push_back({c[5],line,start,end,length});
        }
        need(in.eof(),"PAF read failed");
        std::sort(alignments.begin(),alignments.end(),[](const Alignment&a,const Alignment&b){
            return std::tie(a.contig,a.start,a.line)<std::tie(b.contig,b.start,b.line);});
        std::string contig;uint64_t start=0,end=0;bool counted=false;
        auto flush=[&](uint64_t stop){if(counted&&stop>start)regions[h].push_back({contig,start,stop});};
        for(const auto&a:alignments) {
            if(a.contig!=contig||a.start>=end) {
                flush(end);contig=a.contig;start=a.start;end=a.end;counted=a.length>=50000;
            }else if(a.end>end) {
                flush(a.start);start=end;end=a.end;counted=a.length>=50000;
            }else if(a.end>start){flush(a.start);start=a.end;}
        }
        flush(end);write(prefix+".hap"+std::to_string(h+1)+".bed",regions[h]);
    }
    std::vector<Interval> result;size_t i=0,j=0;
    while(i<regions[0].size()&&j<regions[1].size()) {
        const auto&a=regions[0][i];const auto&b=regions[1][j];
        if(a.contig<b.contig){++i;continue;}if(b.contig<a.contig){++j;continue;}
        const auto start=std::max(a.start,b.start),end=std::min(a.end,b.end);
        if(start<end) {
            if(!result.empty()&&result.back().contig==a.contig&&result.back().end>=start)result.back().end=std::max(result.back().end,end);
            else result.push_back({a.contig,start,end});
        }
        if(a.end<=b.end)++i;else ++j;
    }
    write(prefix+".dip.bed",result);
}
}
