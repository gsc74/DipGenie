#pragma once
#include "SequenceEvidence.hpp"
#include "ControlCalibration.hpp"
#include "DosageConfidence.hpp"
#include "CompleteDosage.hpp"
#include "CopyLevels.hpp"
#include "PanelPhase.hpp"

namespace hic {
struct CountLibrary {
    std::string name,kind;
    std::vector<uint64_t> counts;
    std::vector<uint64_t> occurrences; // optional: exact occurrences, no per-end dedup
};
// Shared by construction and provenance diagnostics: allocate exactly one
// fixed-point vote per informative physical pair, after target remapping.
// het==nullptr: equal shares a_fxy over the molecule's string pairs (v1).
// het!=nullptr: shares proportional to t_x*t_y, the pair's heterozygosity
// confidence (fragments consist of heterozygous alleles, as in HapCUT2), so
// pairs of strings that cannot be phase-informative do not dilute the vote.
// Either way sum_xy a_fxy = 1; a molecule without any heterozygous-capable
// pair is omitted.
template<class F> inline bool molecule_pair_units(const MoleculeMarkers&fragment,
        const std::vector<size_t>&remap,F emit,const std::vector<Score>*het=nullptr) {
    std::vector<std::pair<size_t,size_t>> pairs;
    for(size_t a:fragment.ends[0])for(size_t b:fragment.ends[1]) {
        size_t x=remap.at(a),y=remap.at(b);
        if(x==SIZE_MAX||y==SIZE_MAX||x==y)continue;
        if(x>y)std::swap(x,y);
        pairs.emplace_back(x,y);
    }
    std::sort(pairs.begin(),pairs.end());pairs.erase(std::unique(pairs.begin(),pairs.end()),pairs.end());
    if(pairs.empty())return false;
    uint64_t hash=1469598103934665603ULL;
    for(unsigned char c:fragment.id){hash^=c;hash*=1099511628211ULL;}
    const size_t offset=hash%pairs.size();
    constexpr Score units=1024;
    Score assigned=0;
    if(het) {
        std::vector<Score> p(pairs.size());Score total=0;
        for(size_t i=0;i<pairs.size();++i){p[i]=score_product(het->at(pairs[i].first),het->at(pairs[i].second));total=score_sum(total,p[i]);}
        if(!total)return false;
        std::vector<Score> share(pairs.size());Score floor_sum=0;
        for(size_t i=0;i<pairs.size();++i){share[i]=Score((__int128(units)*p[i])/total);floor_sum+=share[i];}
        // Deterministic leftover units, in molecule-hash order, to positive pairs.
        for(size_t i=0;floor_sum<units;i=(i+1)%pairs.size()) {
            const size_t j=(offset+i)%pairs.size();
            if(p[j]){++share[j];++floor_sum;}
        }
        for(size_t i=0;i<pairs.size();++i)if(share[i]){emit(pairs[i],share[i]);assigned+=share[i];}
    } else {
        const Score base=units/Score(pairs.size()),remainder=units%Score(pairs.size());
        for(size_t i=0;i<pairs.size();++i) {
            const Score weight=base+(i<size_t(remainder));if(!weight)continue;
            emit(pairs[(offset+i)%pairs.size()],weight);assigned+=weight;
        }
    }
    if(assigned!=units)throw std::logic_error("Physical-molecule weight not conserved");
    return true;
}
struct PreparedModel {
    PresenceModel model;
    hic_calibration::Calibration calibration;
    uint64_t calibration_controls=0;
    std::vector<size_t> target_ids;
    // Exact linkage eligibility from target construction: absent-dose words
    // excluded, and copy-expanded words map to their FIRST threshold only.
    std::vector<size_t> linkage_remap;
    size_t hom=0,het=0,informative_molecules=0;
    size_t presence_only_targets=0;
    Score contact_units=0;
    Score read_contact_units=0;
    Score raw_molecule_units=0;
    PanelSummary panel;
    std::vector<hic_calibration::Calibration> library_calibrations;
    std::vector<size_t> dosage_libraries;
    std::vector<uint64_t> dosage_counts;
    bool sr_dosage_only=false;
    size_t hic_only_targets_excluded=0;
    static constexpr Score confidence_scale=1024;
    Score hom_reward_units=0,het_reward_units=0;
    Score absent_reward_units=0;
    hic_calibration::CompleteDosage observed_dosage;
    size_t zero_dosage_targets=0,hic_only_targets_retained=0;
    size_t error_aware_targets=0;
    std::vector<long double> hom_probabilities;
    std::vector<long double> absent_probabilities;
    // 0 = ordinary target; k>0 = ">= k copies" level of a single-segment
    // multi-copy string (CopyLevels.hpp). Several levels share one string.
    std::vector<uint32_t> copy_level;
    // Optional stage-2 (orientation-only) molecule factors; see Inference.
    std::vector<PresenceContact> phase_contacts;
    size_t phase_informative_molecules=0;Score phase_raw_molecule_units=0;
    std::vector<uint64_t> occurrence_counts; // dosage libraries; copy levels only
    size_t copy_strings=0,copy_level_targets=0;
};
// Highest copy level for a single-segment string with variable panel copy
// number; 0 when the string is scored as an ordinary target.
inline uint32_t copy_levels(const WordTarget&t) {
    if(t.unique_locus||t.occurrences.size()!=1||t.copy_histogram.empty())return 0;
    const uint32_t top=t.copy_histogram.back().first;
    return top>=2&&top<=64?top:0;
}
inline PreparedModel build_model(const SequenceIndex&index,const RawEvidence&raw,int H,int blocks,int budget,int threads,
                                const std::vector<CountLibrary>&libraries={},Score phase_weight=1,bool het_vote=false,bool phase_stage=false,
                                bool panel_prior=false,long double log_odds_cap=0,bool het_markers=false,bool panel_relations=true,
                                bool prepare_read_contacts=true) {
    if(phase_weight<1)throw std::invalid_argument("Phase weight");
    PreparedModel result;auto&m=result.model;m.H=H;m.blocks=blocks;m.budget=budget;
    RowExecutor pool(threads);
    std::vector<size_t> remap(index.targets.size(),SIZE_MAX);
    (void)checked_pair_domain(size_t(H));
    const DonorMask all=DonorMask::full(size_t(H));
    // Controlled library-role model: SR supplies dosage when available;
    // all retained SR and Hi-C molecule links still supply phase evidence.
    // Hi-C-only inference retains its original count likelihood. No depth
    // threshold, truth label, or per-run objective flag is used.
    for(const auto&lib:libraries) {
        if(lib.counts.size()!=index.targets.size())throw std::invalid_argument("Library count dimensions");
        if(lib.kind!="SR"&&lib.kind!="HIC")throw std::invalid_argument("Library kind");
        result.sr_dosage_only|=lib.kind=="SR";
    }
    for(size_t l=0;l<libraries.size();++l)
        if(!result.sr_dosage_only||libraries[l].kind=="SR")result.dosage_libraries.push_back(l);
    const auto dosage_count=[&](size_t x) {
        if(libraries.empty())return index.targets[x].count;
        uint64_t count=0;
        for(size_t l:result.dosage_libraries) {
            if(libraries[l].counts[x]>UINT64_MAX-count)throw std::overflow_error("Dosage count");
            count+=libraries[l].counts[x];
        }
        return count;
    };
    hic_calibration::Controls controls;uint64_t zero_controls=0;
    for(size_t x=0;x<index.targets.size();++x)
        if(index.targets[x].unique_locus&&index.targets[x].carriers==all){controls.add(dosage_count(x));zero_controls+=dosage_count(x)==0;}
    result.calibration=hic_calibration::Calibration::fit(controls);
    result.calibration_controls=controls.n;
    std::map<uint64_t,uint64_t> dosage_histogram;
    for(size_t x=0;x<index.targets.size();++x) {
        const auto&t=index.targets[x];const auto n=dosage_count(x);
        if(t.unique_locus&&t.carriers!=all&&t.occurrences.size()==1)++dosage_histogram[n];
    }
    result.observed_dosage=hic_calibration::CompleteDosage::fit(result.calibration,dosage_histogram,zero_controls,controls.n);
    for(const auto&lib:libraries) {
        if(lib.counts.size()!=index.targets.size())throw std::invalid_argument("Library count dimensions");
        hic_calibration::Controls c;
        for(size_t x=0;x<index.targets.size();++x)
            if(index.targets[x].unique_locus&&index.targets[x].carriers==all)c.add(lib.counts[x]);
        result.library_calibrations.push_back(hic_calibration::Calibration::fit(c));
    }
    std::vector<hic_calibration::Calibration> calibrations;
    if(libraries.empty())calibrations.push_back(result.calibration);
    else for(size_t l:result.dosage_libraries)calibrations.push_back(result.library_calibrations[l]);
    for(size_t x=0;x<index.targets.size();++x) {
        auto&t=index.targets[x];
        if(!libraries.empty()) {
            uint64_t total=0;
            for(const auto&lib:libraries) {
                if(lib.counts[x]>UINT64_MAX-total)throw std::overflow_error("Combined library counts");
                total+=lib.counts[x];
            }
            if(total!=t.count)throw std::logic_error("Combined count mismatch");
        }
        // Shared, single-locus conserved strings are calibration controls,
        // not discriminative coverage targets, as in v1's anchor filtering.
        if(t.unique_locus&&t.carriers==all)continue;
        if(!dosage_count(x)){
            // A unique single-segment locus has an identifiable zero-copy
            // alternative. Repeated/cross-segment zeros remain unclassified.
            if(!t.unique_locus||t.occurrences.size()!=1){if(t.count)++result.hic_only_targets_excluded;continue;}
            ++result.zero_dosage_targets;result.hic_only_targets_retained+=t.count>0;
        }
        // Absent observations enter dosage, not the molecule-marker map.
        // Otherwise a Hi-C-only zero-SR word would dilute the fixed physical
        // fragment vote, changing phase evidence in this dosage-only test.
        if(dosage_count(x))remap[x]=result.target_ids.size();
        const uint32_t levels=dosage_count(x)?copy_levels(t):0;
        result.target_ids.push_back(x);result.copy_level.push_back(levels?1:0);
        if(levels)++result.copy_strings;
        for(uint32_t k=2;k<=levels;++k){result.target_ids.push_back(x);result.copy_level.push_back(k);++result.copy_level_targets;}
    }
    m.targets.resize(result.target_ids.size());
    result.hom_probabilities.resize(result.target_ids.size());
    result.absent_probabilities.resize(result.target_ids.size());
    result.dosage_counts.resize(result.target_ids.size());
    result.occurrence_counts.assign(result.target_ids.size(),0);
    constexpr size_t target_chunk=256;
    pool.for_each((result.target_ids.size()+target_chunk-1)/target_chunk,[&](size_t chunk) {
      std::vector<uint64_t> counts(calibrations.size());
      const size_t last=std::min(result.target_ids.size(),(chunk+1)*target_chunk);
      for(size_t target=chunk*target_chunk;target<last;++target) {
        const size_t x=result.target_ids[target];const auto&t=index.targets[x];
        if(libraries.empty())counts[0]=t.count;
        else for(size_t l=0;l<result.dosage_libraries.size();++l)
            counts[l]=libraries[result.dosage_libraries[l]].counts[x];
        result.dosage_counts[target]=dosage_count(x);
        auto hist=t.copy_histogram;
        // Unique-locus synthetic fixtures can derive their exact histogram
        // from carrier bits. Never guess multiplicity for a repeated string.
        if(hist.empty()&&t.unique_locus)hist={{1,uint32_t(t.carriers.count())}};
        uint64_t carriers=0;for(auto[k,n]:hist){(void)k;carriers+=n;}
        if(carriers!=t.carriers.count())throw std::logic_error("Panel copy/carrier mismatch");
        // A string found in several common-anchor segments is not one
        // allelic locus. Recombination can place copies from distinct donors
        // on the same inferred haplotype, even when every intact donor has
        // only one copy. Its intact-panel copy histogram therefore cannot
        // identify haplotype dosage for all paths searched by this model.
        // Retain its sequence-presence support with equal one-/two-haplotype
        // rewards. This is dosage abstention, not a new posterior estimate,
        // a depth threshold, or removal of the observed string.
        const auto confidence=t.occurrences.size()>1?
            hic_calibration::DosageConfidence{0.5L,0.5L,0.L}:
            hic_calibration::dosage_confidence(hist,H,calibrations,counts);
        // Expected reward for local targets; dosage-neutral support otherwise.
        // All three integer dosage weights sum exactly to the confidence scale.
        std::array<long double,3> probabilities{0.L,confidence.het_probability,confidence.hom_probability};
        auto occurrences=t.occurrences;
        if(t.unique_locus&&t.occurrences.size()==1) {
            if(panel_prior) {
                // Frequency-aware prior: Hardy-Weinberg dosage of the string's
                // panel frequency f (carriers among the H panel paths, with a
                // half-count pseudocount); no truth or inferred path enters.
                const long double f=(static_cast<long double>(t.carriers.count())+.5L)/(H+1.L);
                probabilities=result.observed_dosage.posterior(result.dosage_counts[target],{(1-f)*(1-f),2*f*(1-f),f*f});
            } else probabilities=result.observed_dosage.posterior(result.dosage_counts[target]);
        }
        if(const uint32_t level=result.copy_level[target]) {
            std::map<uint32_t,uint64_t> classes;
            const uint64_t zero=uint64_t(H)-t.carriers.count();
            if(zero)classes[0]=zero;
            for(auto[c,n]:t.copy_histogram)classes[c]=n;
            uint64_t copy_count=0;
            if(libraries.empty())copy_count=t.occurrence_count;
            else for(size_t l:result.dosage_libraries) {
                const auto&v=libraries[l].occurrences;
                if(v.size()!=index.targets.size())throw std::invalid_argument("Copy-level targets require occurrence counts");
                copy_count+=v[x];
            }
            result.occurrence_counts[target]=copy_count;
            const auto levels=hic_calibration::copy_level_posterior(result.observed_dosage,copy_count,classes);
            probabilities=levels.at(level-1);
            DonorMask mask;
            for(size_t h=0;h<size_t(H);++h)if(t.copies(h)>=level)mask.set(h);
            if(!mask)throw std::logic_error("Empty copy level");
            occurrences={{t.occurrences.front().first,mask}};
        }
        const auto units=log_odds_cap>0?hic_calibration::log_odds_units(probabilities,log_odds_cap,PreparedModel::confidence_scale):
                                        hic_calibration::dosage_units(probabilities);
        const Score h=units[2],tweight=units[1],absent=units[0];
        if(h<0||h>PreparedModel::confidence_scale)
            throw std::overflow_error("Genotype confidence units");
        const bool hom=h>=tweight; // reporting only; does not filter phase factors
        m.targets[target]={hom,occurrences,hom?h:tweight,hom?tweight:h,absent};
        result.hom_probabilities[target]=probabilities[2];
        result.absent_probabilities[target]=probabilities[0];
      }
    });
    for(const auto&t:m.targets) {
        if(t.occurrences.size()>1)++result.presence_only_targets;
        const Score h=t.hom?t.reward_units:t.other_reward_units;
        const Score tweight=t.hom?t.other_reward_units:t.reward_units;
        if(result.hom_reward_units>INT64_MAX-h||result.het_reward_units>INT64_MAX-tweight)
            throw std::overflow_error("Genotype confidence units");
        result.hom_reward_units+=h;result.het_reward_units+=tweight;
        result.absent_reward_units+=t.absent_reward_units;
        if(t.hom)++result.hom;else ++result.het;
    }
    result.error_aware_targets=result.observed_dosage.observations;
    // An all-conserved panel can legitimately have no discriminative target.
    // Its objective is identically zero; the DP returns a valid panel pair
    // with a zero/zero optimality certificate, not a missing-path failure.
    using PairKey=std::pair<size_t,size_t>;
    struct PairHash {
        size_t operator()(const PairKey&p)const {
            uint64_t x=uint64_t(p.first)*0x9e3779b97f4a7c15ULL+uint64_t(p.second);
            x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;x^=x>>31;
            return size_t(x);
        }
    };
    struct Contribution {PairKey pair;Score weight;};
    struct ContactTask {std::vector<std::vector<Contribution>> bins;size_t molecules=0;};
    constexpr Score units=1024;
    constexpr size_t molecules_per_task=64;
    const size_t owners=std::min(size_t(256),size_t(threads)*4);
    const size_t batch_size=molecules_per_task*std::min(size_t(256),size_t(threads)*4);
    std::vector<Score> het_units(m.targets.size());
    for(size_t x=0;x<m.targets.size();++x){const auto&t=m.targets[x];het_units[x]=t.hom?t.other_reward_units:t.reward_units;}
    // het_markers: a molecule end keeps a string as a linkage marker only if its
    // read-count dosage posterior favours one copy (t_x is the largest of the
    // three dosage weights). Strings that both haplotypes carry are then not
    // markers, however the pair places them; data only, no truth or pair.
    std::vector<unsigned char> marker(m.targets.size(),1);
    if(het_markers)for(size_t x=0;x<m.targets.size();++x){
        const auto&t=m.targets[x];const Score h=t.hom?t.reward_units:t.other_reward_units;
        marker[x]=het_units[x]>0&&het_units[x]>=h&&het_units[x]>=t.absent_reward_units;
    }
    using Contact=typename std::decay_t<decltype(m.contacts)>::value_type;
    struct Pooled {std::vector<Contact> contacts;size_t molecules=0;Score raw=0;};
    // One pooling pass: vote allocation (nullptr = equal shares) and an optional
    // heterozygous-marker filter applied to the fragment ends before pairing.
    auto pool_contacts=[&](const std::vector<Score>*vote) {
    Pooled res;
    std::vector<std::unordered_map<PairKey,Score,PairHash>> pooled(owners);
    for(size_t first=0;first<raw.markers.size();first+=batch_size) {
        const size_t finish=std::min(raw.markers.size(),first+batch_size);
        std::vector<ContactTask> tasks((finish-first+molecules_per_task-1)/molecules_per_task);
        pool.for_each(tasks.size(),[&](size_t task) {
            auto&out=tasks[task];out.bins.resize(owners);
            const size_t end=std::min(finish,first+(task+1)*molecules_per_task);
            for(size_t n=first+task*molecules_per_task;n<end;++n) {
                const MoleculeMarkers*fragment=&raw.markers[n];MoleculeMarkers filtered;
                if(het_markers) {
                    filtered.id=fragment->id;
                    for(int e=0;e<2;++e)for(size_t x:fragment->ends[e]){
                        if(het_markers&&(remap.at(x)==SIZE_MAX||!marker[remap.at(x)]))continue;
                        filtered.ends[e].push_back(x);
                    }
                    fragment=&filtered;
                }
                if(molecule_pair_units(*fragment,remap,[&](const PairKey&pair,Score weight){
                    out.bins[PairHash{}(pair)%owners].push_back({pair,weight});
                },vote))++out.molecules;
            }
            for(auto&bin:out.bins)
                std::sort(bin.begin(),bin.end(),[](const auto&a,const auto&b){return a.pair<b.pair;});
        });
        pool.for_each(owners,[&](size_t owner) {
            auto&dest=pooled[owner];
            for(const auto&task:tasks) {
                const auto&bin=task.bins[owner];
                for(size_t a=0;a<bin.size();) {
                    size_t b=a+1;Score weight=bin[a].weight;
                    while(b<bin.size()&&bin[b].pair==bin[a].pair)weight+=bin[b++].weight;
                    dest[bin[a].pair]+=weight;a=b;
                }
            }
        });
        for(const auto&task:tasks) {
            res.molecules+=task.molecules;
            res.raw+=Score(task.molecules)*units;
        }
    }
    // Confidence is still applied only after all molecules have been pooled.
    // Exact integer reduction and a final key sort preserve the former score
    // and factor ordering independently of worker scheduling.
    std::vector<std::vector<Contact>> contacts(owners);
    pool.for_each(owners,[&](size_t owner) {
        auto&out=contacts[owner];out.reserve(pooled[owner].size());
        for(const auto&[pair,weight]:pooled[owner]) {
            const Score weighted=Score((__int128(weight)*het_units[pair.first]*het_units[pair.second])/(units*units));
            if(weighted)out.push_back({pair.first,pair.second,weighted});
        }
        decltype(pooled)::value_type{}.swap(pooled[owner]);
        std::sort(out.begin(),out.end(),[](const auto&a,const auto&b){return std::make_pair(a.x,a.y)<std::make_pair(b.x,b.y);});
    });
    while(contacts.size()>1) {
        std::vector<std::vector<Contact>> merged((contacts.size()+1)/2);
        pool.for_each(merged.size(),[&](size_t i) {
            auto&a=contacts[2*i];
            if(2*i+1==contacts.size()){merged[i]=std::move(a);return;}
            auto&b=contacts[2*i+1];auto&dest=merged[i];dest.reserve(a.size()+b.size());
            std::merge(a.begin(),a.end(),b.begin(),b.end(),std::back_inserter(dest),
                [](const auto&x,const auto&y){return std::make_pair(x.x,x.y)<std::make_pair(y.x,y.y);});
            std::vector<Contact>().swap(a);std::vector<Contact>().swap(b);
        });
        contacts=std::move(merged);
    }
    res.contacts=std::move(contacts.front());
    return res;
    };
    if(prepare_read_contacts) {
        auto stage=pool_contacts(het_vote?&het_units:nullptr);
        m.contacts=std::move(stage.contacts);result.informative_molecules=stage.molecules;result.raw_molecule_units=stage.raw;
    }
    if(phase_stage&&prepare_read_contacts) {
        // Stage-2 linked-observation factors (orientation only): heterozygosity-
        // proportional votes over genome-unique markers.
        auto stage=pool_contacts(&het_units);
        result.phase_contacts=std::move(stage.contacts);
        result.phase_informative_molecules=stage.molecules;result.phase_raw_molecule_units=stage.raw;
    }
    for(const auto&contact:m.contacts)result.contact_units+=contact.units;
    if(result.contact_units>result.raw_molecule_units)throw std::logic_error("Confidence creates molecule mass");
    result.read_contact_units=result.contact_units;
    if(panel_relations)result.panel=add_panel_phase(m,index,result.target_ids,&result.dosage_counts);
    result.contact_units+=result.panel.units;
    // Diagnostic ablation only: phase_weight scales total phase mass relative
    // to total heterozygous genotype mass (1 = the documented objective).
    m.weights=normalized_weights(result.contact_units,result.het_reward_units,phase_weight);m.check();
    result.linkage_remap=std::move(remap);
    return result;
}
}
