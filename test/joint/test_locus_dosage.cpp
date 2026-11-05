#include "BuildModel.hpp"
#include <cassert>
#include <iostream>

int main() {
    using namespace hic;
    for(uint64_t count:{1,20,40,400})for(int threads:{1,2,4}) {
        SequenceIndex index;
        WordTarget control;control.word=0;control.count=40;control.carriers=3;
        control.occurrences={{0,3}};control.copy_histogram={{1,2}};
        WordTarget local;local.word=1;local.count=count;local.carriers=1;
        local.occurrences={{0,1}};local.copy_histogram={{1,1}};
        WordTarget cross;cross.word=2;cross.count=count;cross.carriers=3;
        cross.unique_locus=false;cross.occurrences={{0,1},{1,2}};
        cross.copy_histogram={{1,2}};
        index.targets={control,local,cross};RawEvidence raw;
        const auto p=build_model(index,raw,2,2,2,threads);
        assert((p.target_ids==std::vector<size_t>{1,2}));
        assert(p.presence_only_targets==1);
        const auto&t=p.model.targets[1];
        assert(t.reward_units==512&&t.other_reward_units==512);
        assert(p.hom_probabilities[1]==0.5L);
        for(unsigned mask:{1u,2u,3u})assert(PresenceModel::genotype_reward(t,mask)==512);
        assert(PresenceModel::genotype_reward(t,0)==0);
        const auto expected=p.observed_dosage.posterior(count);
        assert(p.hom_probabilities[0]==expected[2]&&p.absent_probabilities[0]==expected[0]);
        // All 16 two-block diploid paths preserve the exact generic score
        // and its universal upper bound, including two copies on one path.
        for(int a=0;a<4;++a)for(int b=0;b<4;++b) {
            std::vector<int> path{a,b};const auto score=p.model.score(path);
            assert(score<=p.model.universal_bound());
            assert(score==p.model.score({(a%2)*2+a/2,(b%2)*2+b/2}));
        }
        // Intact pair: 1 copy on each. Recombinant pair: 2 copies on hap1,
        // none on hap2. Neither gets artificial extra dosage reward.
        const auto presence=[&](const std::vector<int>&path){unsigned char v=0;for(const auto&[b,mask]:t.occurrences)v|=p.model.signature(mask,path[b]);return v;};
        assert(presence({1,1})==3&&presence({1,2})==1);
        assert(PresenceModel::genotype_reward(t,presence({1,1}))==PresenceModel::genotype_reward(t,presence({1,2})));
    }
    std::cout<<"PASS cross-segment dosage abstention; local three-state confidence; all observations retained; copy counterexample; full score/bound and global haplotype-exchange parity\n";
}
