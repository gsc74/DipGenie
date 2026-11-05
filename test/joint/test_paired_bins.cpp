#include "TestTemp.hpp"
#include "PairedReadBinning.hpp"
#include <cassert>
#include <random>
#include <iostream>
#include <unistd.h>
using namespace hic;
static std::string contents(const std::filesystem::path&p) {
    std::ifstream in(p);return std::string(std::istreambuf_iterator<char>(in),{});
}
int main() {
    std::string pattern_buffer=test_temp_pattern("dipgenie-paired-bins");char*pattern=pattern_buffer.data();const char*root=mkdtemp(pattern);
    if(!root)throw std::runtime_error("Temporary directory");
    const std::filesystem::path dir(root);
    std::mt19937 rng(1151);std::vector<ComponentReference> refs={{"chr1",""},{"chr6",""}};
    for(auto&r:refs){r.sequence.resize(6000);for(char&c:r.sequence)c="ACGT"[rng()%4];}
    const auto&a=refs[0].sequence;const auto&b=refs[1].sequence;
    std::vector<std::pair<std::string,std::string>> pairs={
        {a.substr(100,250),reverse_complement(a.substr(1500,250))},
        {b.substr(150,250),b.substr(2000,250)},
        {a.substr(300,250),b.substr(300,250)},
        {a.substr(500,250),std::string(250,'N')},
        {std::string(250,'N'),std::string(250,'N')},
        {a.substr(2500,250)+b.substr(1000,250),a.substr(4000,250)}
    };
    std::ofstream fq[2]={std::ofstream(dir/"R1.fq"),std::ofstream(dir/"R2.fq")};
    for(size_t i=0;i<pairs.size();++i) {
        const std::string*reads[2]={&pairs[i].first,&pairs[i].second};
        for(size_t end=0;end<2;++end)fq[end]<<"@r"<<i<<'/'<<end+1<<'\n'<<*reads[end]<<"\n+\n"<<std::string(reads[end]->size(),'I')<<'\n';
    }
    for(auto&f:fq)f.close();
    ReadLibrary library{"hic","HIC",(dir/"R1.fq").string(),(dir/"R2.fq").string()};
    for(int threads:{1,2,4,8}) {
        PairedReadBinner binner(refs,threads);
        const auto out=dir/("native"+std::to_string(threads));
        auto bins=binner.partition(library,out);
        // The chimeric end is explicitly quarantined. Both mates are preserved byte-for-byte.
        assert(bins.input_pairs==6&&bins.cross_component_pairs==1&&bins.unassigned_pairs==3);
        assert(bins.summaries[0].pairs==1&&bins.summaries[1].pairs==1);
        FastqInput r1(bins.libraries[0].r1),r2(bins.libraries[0].r2);FastqRead x,y;
        assert(r1.next(x)&&r2.next(y)&&paired_molecule_id(x,y)=="r0");
        assert(x.sequence==pairs[0].first&&y.sequence==pairs[0].second&&!r1.next(x)&&!r2.next(y));
        if(threads>1)for(auto bin:{"component-0","component-1","cross-component","unassigned"})
            for(auto file:{"R1.fastq","R2.fastq"})assert(contents(out/bin/file)==contents(dir/"native1"/bin/file));
    }
    library.kind="SR";
    PairedReadBinner sr(refs,4);auto s=sr.partition(library,dir/"sr");
    uint64_t total=s.cross_component_pairs+s.unassigned_pairs;
    for(const auto&bin:s.summaries)total+=bin.pairs;
    assert(s.input_pairs==6&&total==6);
    bool rejected=false;try{(void)sr.partition({"x","ONT",library.r1,library.r2},dir/"lr");}
    catch(const std::invalid_argument&){rejected=true;}assert(rejected);
    // Strict identity: duplicates are invalid within a library, but the same
    // names in distinct libraries are independent observations.
    PairedReadBinner native(refs,2);
    auto independent=library;independent.name="second-library";
    assert(native.partition(independent,dir/"native-second").input_pairs==6);
    std::ofstream duplicates[2]={std::ofstream(dir/"duplicate1.fq"),std::ofstream(dir/"duplicate2.fq")};
    for(size_t side=0;side<2;++side)for(int rep=0;rep<2;++rep)
        duplicates[side]<<"@same/"<<side+1<<'\n'<<a.substr(20+side*350,150)<<"\n+\n"<<std::string(150,'I')<<'\n';
    for(auto&f:duplicates)f.close();
    rejected=false;
    try{native.partition({"dup","HIC",(dir/"duplicate1.fq").string(),(dir/"duplicate2.fq").string()},dir/"duplicate-bin");}
    catch(const std::runtime_error&){rejected=true;}
    assert(rejected);
    // More bins than the descriptor cache, revisited in a second pass.
    std::vector<ComponentReference> many;
    for(size_t i=0;i<70;++i) {
        std::string sequence(1000,'A');for(char&c:sequence)c="ACGT"[rng()%4];
        many.push_back({"contig"+std::to_string(i),std::move(sequence)});
    }
    std::ofstream large[2]={std::ofstream(dir/"many1.fq"),std::ofstream(dir/"many2.fq")};
    for(int pass=0;pass<2;++pass)for(size_t c=0;c<many.size();++c)for(size_t side=0;side<2;++side) {
        const auto sequence=many[c].sequence.substr(100+side*300,250);
        large[side]<<"@p"<<pass<<"c"<<c<<'/'<<side+1<<'\n'<<sequence<<"\n+\n"<<std::string(250,'I')<<'\n';
    }
    for(auto&f:large)f.close();
    PairedReadBinner cached(many,4);
    const auto many_bins=cached.partition({"many","SR",(dir/"many1.fq").string(),(dir/"many2.fq").string()},dir/"many");
    assert(many_bins.input_pairs==140&&many_bins.unassigned_pairs==0&&many_bins.cross_component_pairs==0);
    for(size_t c=0;c<many.size();++c) {
        assert(many_bins.summaries[c].pairs==2);
        FastqInput reads(many_bins.libraries[c].r1);FastqRead read;size_t n=0;
        while(reads.next(read))++n;
        assert(n==2);
    }
    std::filesystem::remove_all(dir);
    std::cout<<"PASS paired component routing: native placement, reverse strands, cross-component and missing-end quarantine, quality preservation, no duplicated molecules, thread parity, LR rejection, descriptor cache\n";
}
