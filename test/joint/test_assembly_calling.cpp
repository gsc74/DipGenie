#include "AssemblyCalling.hpp"
#include "VariantOutput.hpp"
#include <iostream>
#include <stdexcept>
int main(int argc,char**argv)try {
    using namespace dipcall_vcf::embedded;
    if(argc==5&&std::string(argv[1])=="normalize")dipcall_vcf::output::normalize(argv[2],argv[3],argv[4]);
    else if(argc==4&&std::string(argv[1])=="publish")dipcall_vcf::output::publish(argv[2],argv[3]);
    else if(argc>=4&&std::string(argv[1])=="combine")dipcall_vcf::output::combine(std::vector<std::string>(argv+3,argv+argc),argv[2]);
    else if(argc==5&&std::string(argv[1])=="pair")pair_vcf(argv[2],argv[3],argv[4]);
    else if(argc==5&&std::string(argv[1])=="bed")callable_bed({argv[2],argv[3]},argv[4]);
    else if(argc==6&&std::string(argv[1])=="pileup")pileup(argv[2],{argv[3],argv[4]},argv[5]);
    else if(argc==7&&std::string(argv[1])=="align") {
        const auto stats=align(argv[2],{argv[3],argv[4]},argv[5],std::stoul(argv[6]));
        for(size_t h=0;h<2;++h)std::cout<<h+1<<'\t'<<stats[h].queries<<'\t'<<stats[h].alignments<<'\t'<<stats[h].retained_alignments<<'\n';
    } else throw std::invalid_argument("align REF HAP1 HAP2 PREFIX THREADS | pileup REF BAM1 BAM2 OUTPUT | pair INPUT OUTPUT SAMPLE");
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
