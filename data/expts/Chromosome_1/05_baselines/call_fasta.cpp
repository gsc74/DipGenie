#include "DipcallRunner.hpp"
#include "SequenceEvidence.hpp"
#include <iostream>
int main(int argc,char**argv)try {
    if(argc!=4)throw std::runtime_error("two-haplotypes.fa reference.fa calls.vcf.gz");
    hic::TextInput in(argv[1]);std::string line;std::vector<std::string> haps;
    while(in.getline(line)) {
        if(line.empty())continue;
        if(line[0]=='>')haps.emplace_back();
        else {if(haps.empty())throw std::runtime_error("Missing FASTA header");haps.back()+=line;}
    }
    if(haps.size()!=2||haps[0].empty()||haps[1].empty())throw std::runtime_error("Exactly two nonempty haplotypes required");
    dipcall_vcf::Config cfg;cfg.reference_fasta=argv[2];cfg.output_vcf=argv[3];cfg.sample="HG002";cfg.threads=64;
    std::cout<<dipcall_vcf::run(cfg,haps)<<'\n';
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
