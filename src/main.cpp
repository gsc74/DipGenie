#include "Inference.hpp"
#include "InputGraph.hpp"
#include "DipcallRunner.hpp"
#include "VariantOutput.hpp"
#include <chrono>
#include <iostream>
#include <set>

namespace fs=std::filesystem;
int main(int argc,char**argv)try {
    std::map<std::string,std::string> args;
    const std::set<std::string> flags={"--sr","--hic"};
    const std::set<std::string> values={"--gfa","--gindex","--vcf","--reference","--REF","--out","--write-index",
        "--r1","--r2","--h1","--h2","--sample","--output-vcf","-t","-R"};
    for(int i=1;i<argc;++i) {
        std::string key=argv[i],attached;
        if(key.size()>2&&(key.rfind("-t",0)==0||key.rfind("-R",0)==0)) {
            attached=key.substr(2);key.resize(2);
        }
        if(key=="--help"||key=="-h") {
            std::cout<<"DipGenie: joint diploid inference from paired SR and/or Hi-C reads\n"
                "  --gfa graph.gfa --REF reference_path   or --vcf panel.vcf.gz --reference ref.fa\n"
                "  --write-index panel.gindex            build a reusable graph/reference index\n"
                "  --gindex panel.gindex --out run -t 8 -R 18\n"
                "  --sr -r1 reads.R1.fq.gz -r2 reads.R2.fq.gz\n"
                "  --hic -h1 hic.R1.fq.gz -h2 hic.R2.fq.gz\n"
                "  Paired reads are localized on the graph; no external read mapper is required.\n"
                "  Combine --sr and --hic for both libraries. GFA supports cycles.\n"
                "  --output-vcf calls.vcf.gz --sample HG002 also writes phased variant calls.\n"
                "  -R is a shared recombination budget per component; default 18.\n";
            return 0;
        }
        if(key=="--version"){std::cout<<"DipGenie joint-linkage 2.0\n";return 0;}
        if(key=="-r1"||key=="-r2"||key=="-h1"||key=="-h2")key="-"+key;
        if(!flags.count(key)&&!values.count(key))throw std::invalid_argument("Unknown option: "+key);
        std::string value="1";
        if(!flags.count(key)) {
            if(!attached.empty())value=attached;
            else {if(++i==argc)throw std::invalid_argument("Missing value: "+key);value=argv[i];}
        }
        if(!args.emplace(key,value).second)throw std::invalid_argument("Duplicate option: "+key);
    }
    if(args.count("--gfa")+args.count("--vcf")+args.count("--gindex")!=1)
        throw std::invalid_argument("Choose exactly one of --gfa, --vcf or --gindex");
    auto positive=[&](const std::string&key,int fallback,int minimum) {
        if(!args.count(key))return fallback;
        size_t end=0;int value=std::stoi(args.at(key),&end);
        if(end!=args.at(key).size()||value<minimum)throw std::invalid_argument("Invalid value: "+key);
        return value;
    };
    const int threads=positive("-t",1,1);(void)positive("-R",18,0);
    auto read_graph=[&] {
        const auto start=std::chrono::steady_clock::now();
        std::cerr<<"[M::DipGenie] reading graph input\n";
        auto input=hic::graph_input(args,threads);
        std::cerr<<"[M::DipGenie] graph input ready: components="<<input.references.size()
                 <<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
        return input;
    };
    if(args.count("--write-index")) {
        if(args.count("--sr")||args.count("--hic")||args.count("--out")||args.count("--output-vcf"))
            throw std::invalid_argument("Index construction does not perform inference");
        auto input=read_graph();hic::write_graph_index(input,args.at("--write-index"));
        std::cerr<<"[M::DipGenie] index written: components="<<input.references.size()<<"; explicit full reference stored\n";
        return 0;
    }
    if(!args.count("--out"))throw std::invalid_argument("Inference requires --out");
    const auto out=fs::absolute(args.at("--out"));args["--out"]=out.string();
    if(fs::exists(out))throw std::runtime_error("Output exists; choose a fresh output directory");
    std::vector<hic::ReadLibrary> libraries;
    auto library=[&](const char*flag,const char*a,const char*b,const char*name,const char*kind) {
        if(args.count(flag)) {
            if(!args.count(a)||!args.count(b))throw std::invalid_argument(std::string(flag)+" requires both mates");
            libraries.push_back({name,kind,args.at(a),args.at(b)});
        } else if(args.count(a)||args.count(b))throw std::invalid_argument(std::string(a)+"/"+b+" require "+flag);
    };
    library("--sr","--r1","--r2","sr","SR");library("--hic","--h1","--h2","hic","HIC");
    hic::validate_libraries(libraries);
    if(args.count("--output-vcf")) {
        args["--output-vcf"]=fs::absolute(args.at("--output-vcf")).string();
        const auto&name=args.at("--output-vcf");
        const auto suffix=[&](const std::string&s){return name.size()>=s.size()&&name.compare(name.size()-s.size(),s.size(),s)==0;};
        if(!suffix(".vcf")&&!suffix(".vcf.gz"))throw std::invalid_argument("VCF output must end in .vcf or .vcf.gz");
        if(fs::exists(name)||fs::exists(name+".csi")||fs::exists(name+".tbi"))throw std::runtime_error("VCF output/index exists");
    }
    auto input=read_graph();
    auto load_component=[&](size_t c) {
        const auto start=std::chrono::steady_clock::now();
        std::cerr<<"[M::DipGenie] loading component="<<input.references.at(c).name<<'\n';
        auto graph=input.load(c);
        std::cerr<<"[M::DipGenie] component ready: vertices="<<graph.labels.size()<<" paths="<<graph.paths.size()
                 <<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
        return graph;
    };
    if(input.references.size()==1) {
        auto g=load_component(0);
        if(g.spell(0)!=input.references[0].sequence)throw std::runtime_error("Indexed reference differs from graph");
        // Use the stored reference coordinate name for FASTA-to-VCF calling.
        g.paths[0].name=input.references[0].name;args["--REF"]=g.paths[0].name;
        input.references.clear();return hic::infer_component(args,std::move(libraries),std::move(g));
    }
    fs::create_directories(out);fs::create_directory(out/"read-bins");fs::create_directory(out/"components");
    std::vector<std::vector<hic::ReadLibrary>> binned(input.references.size());
    {
        hic::PairedReadBinner binner(input.references,threads);
        std::ofstream summary(out/"read-bins.tsv");summary<<"library\tcomponent\tpairs\tbases\n";
        for(const auto&lib:libraries) {
            const auto bins=binner.partition(lib,out/"read-bins"/lib.name);
            for(size_t c=0;c<binned.size();++c) {
                const auto&s=bins.summaries[c];summary<<lib.name<<'\t'<<input.references[c].name<<'\t'<<s.pairs<<'\t'<<s.bases<<'\n';
                if(s.pairs)binned[c].push_back(bins.libraries[c]);
            }
            std::cerr<<"[M::DipGenie] binned library="<<lib.name<<" input_pairs="<<bins.input_pairs
                     <<" cross_component="<<bins.cross_component_pairs<<" unassigned="<<bins.unassigned_pairs<<'\n';
        }
        summary.close();if(!summary)throw std::runtime_error("Bin summary write failed");
    }
    // Release full reference strings after assignment, before graph/DP storage.
    for(auto&r:input.references){r.sequence.clear();r.sequence.shrink_to_fit();}
    std::ofstream combined(out/"all-components.fa"),status(out/"components.tsv");std::vector<std::string> vcfs;
    size_t completed=0;
    status<<"component\tstatus\toutput\n";
    for(size_t c=0;c<binned.size();++c) {
        const auto&name=input.references[c].name;
        if(binned[c].empty()){status<<name<<"\tNO_READ_PAIRS\tNULL\n";continue;}
        const auto directory=out/"components"/std::to_string(c);auto options=args;
        options["--out"]=directory.string();options["--REF"]=name;
        if(args.count("--output-vcf"))options["--output-vcf"]=(out/"components"/(std::to_string(c)+".vcf.gz")).string();
        auto graph=load_component(c);graph.paths.at(0).name=name;
        hic::infer_component(options,binned[c],std::move(graph));
        ++completed;
        hic::TextInput fasta((directory/"all-components.fa").string());std::string line;
        while(fasta.getline(line)){if(!line.empty()&&line[0]=='>')line+="#"+name;combined<<line<<'\n';}
        if(args.count("--output-vcf"))vcfs.push_back(options.at("--output-vcf"));
        status<<name<<"\tCOMPLETE\t"<<directory.string()<<'\n';
    }
    combined.close();status.close();if(!combined||!status)throw std::runtime_error("Combined output write failed");
    if(!completed)throw std::runtime_error("No inferred components; approximation_ratio=NULL");
    if(args.count("--output-vcf")) {
        if(vcfs.empty())throw std::runtime_error("No inferred components; approximation_ratio=NULL; no VCF generated");
        dipcall_vcf::output::combine(vcfs,args.at("--output-vcf"));
    }
    return 0;
}catch(const std::exception&e){std::cerr<<"[E::DipGenie] "<<e.what()<<'\n';return 1;}
