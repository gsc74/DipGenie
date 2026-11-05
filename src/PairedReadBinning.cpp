#include "PairedReadBinning.hpp"
#include "NativeLocalizer.hpp"
#include <climits>
#include <array>
#include <cstdlib>
#include <tuple>

namespace hic {
namespace {
struct BinFiles {
    std::ofstream mates[2];
    std::filesystem::path directory;
    explicit BinFiles(const std::filesystem::path&dir,bool resume=false):directory(dir) {
        if(!resume&&!std::filesystem::create_directory(dir))throw std::runtime_error("Read bin exists");
        const auto mode=std::ios::out|(resume?std::ios::app:std::ios::trunc);
        mates[0].open(dir/"R1.fastq",mode);mates[1].open(dir/"R2.fastq",mode);
        if(!mates[0]||!mates[1])throw std::runtime_error("Cannot create paired read bin");
    }
    void put(const std::string&id,const FastqRead&a,const FastqRead&b) {
        const FastqRead*reads[2]={&a,&b};
        for(size_t side=0;side<2;++side)
            mates[side]<<'@'<<id<<'/'<<side+1<<'\n'<<reads[side]->sequence<<"\n+\n"<<reads[side]->quality<<'\n';
        if(!mates[0]||!mates[1])throw std::runtime_error("Paired read bin write failed");
    }
    void finish(){for(auto&f:mates){f.close();if(!f)throw std::runtime_error("FASTQ bin close failed");}}
};
}
class PairedReadBinner::Impl {
public:
    std::vector<std::string> names;
    int threads;
    std::unique_ptr<NativeLocalizer> native;
    Impl(const std::vector<ComponentReference>&refs,int n):threads(n) {
        if(n<1||refs.empty()||refs.size()>INT_MAX)throw std::invalid_argument("Reference binner dimensions");
        std::set<std::string> unique;
        for(const auto&r:refs) {
            if(r.name.empty()||r.sequence.empty()||r.sequence.size()>INT_MAX||!unique.insert(r.name).second)
                throw std::invalid_argument("Empty, oversized or duplicate component reference");
            names.push_back(r.name);
        }
        std::vector<PlacementSequence> contexts;
        for(size_t c=0;c<refs.size();++c) {
            PlacementSequence s;s.name=refs[c].name;s.sequence=refs[c].sequence;
            s.component=s.space=s.reference_contig=uint32_t(c);contexts.push_back(std::move(s));
        }
        native=std::make_unique<NativeLocalizer>(std::move(contexts),NativeLocalizer::Policy{},threads);
    }
};
PairedReadBinner::PairedReadBinner(const std::vector<ComponentReference>&refs,int threads):impl_(std::make_unique<Impl>(refs,threads)){}
PairedReadBinner::~PairedReadBinner()=default;
PairedBins PairedReadBinner::partition(const ReadLibrary&library,const std::filesystem::path&dir) {
    if(library.kind!="SR"&&library.kind!="HIC")throw std::invalid_argument("Paired binning supports SR and Hi-C only");
    if(!std::filesystem::create_directory(dir))throw std::runtime_error("Paired bin directory exists");
    std::ofstream provenance(dir/"source.tsv");
    provenance<<"library\tkind\tr1\tr2\n"<<library.name<<'\t'<<library.kind<<'\t'
              <<library.r1<<'\t'<<library.r2<<'\n';
    provenance.close();if(!provenance)throw std::runtime_error("Binning provenance write failed");
    PairedBins result;const size_t components=impl_->names.size();
    result.summaries.resize(components);result.libraries.resize(components);
    result.cross_component=dir/"cross-component";result.unassigned=dir/"unassigned";
    std::vector<std::unique_ptr<BinFiles>> outputs(components+2);
    std::vector<unsigned char> created(components+2);
    std::vector<uint64_t> touched(components+2);
    size_t open_bins=0;
    for(size_t i=0;i<components;++i) {
        auto dest=dir/("component-"+std::to_string(i));
        result.libraries[i]={library.name,library.kind,(dest/"R1.fastq").string(),(dest/"R2.fastq").string()};
    }
    std::ofstream assignments(dir/"pair-components.tsv");
    assignments<<"molecule\tend1_component\tend2_component\tbin\n";
    FastqInput inputs[2]={FastqInput(library.r1),FastqInput(library.r2)};RowExecutor pool(impl_->threads);
    // Names identify physical observations within a library.
    // Repeated names are rejected, even if identical; equal placements with
    // different names are NOT silently labelled PCR duplicates.
    std::unordered_set<std::string> native_names;
    for(;;) {
        std::vector<FastqRead> reads[2];
        pool.for_each(2,[&](size_t side){FastqRead read;reads[side].reserve(4096);
            while(reads[side].size()<4096&&inputs[side].next(read))reads[side].push_back(std::move(read));});
        const size_t n=reads[0].size();if(n!=reads[1].size())throw std::runtime_error("Binning mate count mismatch");if(!n)break;
        std::vector<std::string> ids(n);std::vector<std::array<int,2>> component(n);
        for(size_t i=0;i<n;++i) {
            ids[i]=paired_molecule_id(reads[0][i],reads[1][i]);
            if(!native_names.insert(ids[i]).second)throw std::runtime_error("Duplicate molecule name in library "+library.name+": "+ids[i]);
        }
        pool.for_each((n+63)/64,[&](size_t task) {
            for(size_t i=task*64;i<std::min(n,(task+1)*64);++i)for(size_t side=0;side<2;++side) {
                const auto placement=impl_->native->locate(reads[side][i].sequence);
                component[i][side]=placement.component_unique()?int(placement.component):-1;
            }
        });
        for(size_t i=0;i<n;++i) {
            ++result.input_pairs;const int a=component[i][0],b=component[i][1];size_t bin;
            if(a<0||b<0){bin=components+1;++result.unassigned_pairs;}
            else if(a!=b){bin=components;++result.cross_component_pairs;}
            else {bin=size_t(a);auto&s=result.summaries[bin];++s.pairs;s.bases+=reads[0][i].sequence.size()+reads[1][i].sequence.size();}
            const auto destination=bin==components?result.cross_component:bin==components+1?result.unassigned:
                                    dir/("component-"+std::to_string(bin));
            if(!outputs[bin]) {
                // Bound descriptor use even for graphs with thousands of components.
                if(open_bins==64) {
                    size_t victim=SIZE_MAX;
                    for(size_t j=0;j<outputs.size();++j)if(outputs[j]&&(victim==SIZE_MAX||touched[j]<touched[victim]))victim=j;
                    outputs.at(victim)->finish();outputs[victim].reset();--open_bins;
                }
                outputs[bin]=std::make_unique<BinFiles>(destination,created[bin]!=0);
                created[bin]=1;++open_bins;
            }
            touched[bin]=result.input_pairs;
            outputs[bin]->put(ids[i],reads[0][i],reads[1][i]);
            assignments<<ids[i]<<'\t'<<(a<0?"*":impl_->names[size_t(a)])<<'\t'<<(b<0?"*":impl_->names[size_t(b)])<<'\t'<<destination.filename().string()<<'\n';
        }
    }
    for(auto&output:outputs)if(output)output->finish();
    assignments.close();if(!assignments)throw std::runtime_error("Pair-assignment audit write failed");
    uint64_t accounted=result.cross_component_pairs+result.unassigned_pairs;
    for(const auto&s:result.summaries)accounted+=s.pairs;
    if(accounted!=result.input_pairs)throw std::logic_error("Paired binning loses or duplicates molecules");
    return result;
}
}
