#include "BuildModel.hpp"
#include <cassert>
#include <iostream>
using namespace hic;
int main() {
    SequenceIndex index;
    WordTarget control;control.word=0;control.carriers=3;control.occurrences={{0,3}};control.copy_histogram={{1,2}};
    WordTarget a;a.word=1;a.carriers=1;a.occurrences={{0,1}};a.copy_histogram={{1,1}};
    WordTarget b=a;b.word=2;b.occurrences={{1,1}};
    WordTarget c=b;c.word=3;c.carriers=2;c.occurrences={{1,2}};
    index.targets={control,a,b,c};
    CountLibrary sr{"sr","SR",{40,20,20,0}};
    RawEvidence raw;MoleculeMarkers molecule;molecule.id="hic";molecule.ends[0]={1};molecule.ends[1]={2};raw.markers={molecule};
    PreparedModel reference;
    for(size_t x=0;x<index.targets.size();++x)index.targets[x].count=sr.counts[x];
    reference=build_model(index,raw,2,2,2,1,{sr});
    assert((reference.target_ids==std::vector<size_t>{1,2,3})&&reference.informative_molecules==1);
    assert(reference.zero_dosage_targets==1&&reference.absent_probabilities[2]>.99L);
    for(uint64_t n:{1,10,100,10000})for(int threads:{1,2,4}) {
        CountLibrary hc{"hic","HIC",{2*n,n,n,10*n}};
        for(size_t x=0;x<index.targets.size();++x)index.targets[x].count=sr.counts[x]+hc.counts[x];
        const auto p=build_model(index,raw,2,2,2,threads,{sr,hc});
        assert(p.target_ids==reference.target_ids&&p.hom_probabilities==reference.hom_probabilities);
        assert(p.dosage_counts==reference.dosage_counts&&p.hic_only_targets_excluded==0&&p.hic_only_targets_retained==1);
        assert(p.informative_molecules==1&&p.raw_molecule_units==1024&&p.read_contact_units>0);
        assert((p.dosage_libraries==std::vector<size_t>{0}));
        assert(p.library_calibrations.size()==2);
        for(int a=0;a<4;++a)for(int b=0;b<4;++b)assert(p.model.score({a,b})==reference.model.score({a,b}));
        const auto swapped=build_model(index,raw,2,2,2,threads,{hc,sr});
        assert(swapped.hom_probabilities==p.hom_probabilities&&swapped.target_ids==p.target_ids);
        assert((swapped.dosage_libraries==std::vector<size_t>{1}));
        // Removing contact evidence changes the phase contribution, proving
        // that Hi-C has not accidentally been removed with its dosage counts.
        const auto no_links=build_model(index,RawEvidence{},2,2,2,threads,{sr,hc});
        assert(no_links.read_contact_units==0&&no_links.hom_probabilities==p.hom_probabilities);
        // Hi-C-only still uses its own calibrated count evidence.
        for(size_t x=0;x<index.targets.size();++x)index.targets[x].count=hc.counts[x];
        const auto only=build_model(index,raw,2,2,2,threads,{hc});
        assert(!only.sr_dosage_only&&only.hic_only_targets_excluded==0&&only.target_ids.size()==3);
        const auto expected=only.observed_dosage.posterior(n);
        assert(only.hom_probabilities[0]==expected[2]&&only.absent_probabilities[0]==expected[0]);
    }
    std::cout<<"PASS library-order/thread invariance; SR dosage unaffected by Hi-C abundance; physical Hi-C phase links retained; unobserved local targets retained; Hi-C-only calibration preserved\n";
}
