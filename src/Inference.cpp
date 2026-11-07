#include "ReadLibraries.hpp"
#include "DipcallRunner.hpp"
#include "GindexGraph.hpp"
#include "PanelSelection.hpp"
#include "ConstantPresence.hpp"
#include "GlobalPhaseMove.hpp"
#include "FreeOrientation.hpp"
#include "ProgressLog.hpp"
#include "NativeReadEvidence.hpp"
#include "NativeFactors.hpp"
#include <chrono>
#include <filesystem>
#include <iostream>
#include <iomanip>

#include "Inference.hpp"

int hic::infer_component(std::map<std::string,std::string> args, std::vector<ReadLibrary> libraries, NativeGraph graph) {
    const std::string sample=args.count("--sample")?args.at("--sample"):"SAMPLE";
    if(sample.empty()||sample.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-")!=std::string::npos)
        throw std::invalid_argument("Invalid sample name");
    hic::validate_libraries(libraries);
    auto integer=[&](const std::string&key,int fallback) {
        if(!args.count(key))return fallback;
        size_t end=0;int n=std::stoi(args.at(key),&end);
        if(end!=args.at(key).size())throw std::invalid_argument("Invalid integer: "+key);
        return n;
    };
    const int threads=integer("-t",1),budget=integer("-R",18);
    // Hi-C ends that cross a ligation junction keep their 5'-most unique
    // placement (NativeReadEvidence.hpp).
    if(threads<1||budget<0)throw std::invalid_argument("Invalid threads or recombination budget");
    if(args.count("--output-vcf")) {
        if(!args.count("--REF")&&!args.count("--gindex"))throw std::invalid_argument("--output-vcf requires --REF panel path name/prefix for GFA");
        if(std::filesystem::exists(args["--output-vcf"]))throw std::runtime_error("Refusing to overwrite existing VCF");
    }
    std::filesystem::path out=args.at("--out");
    if(std::filesystem::exists(out))throw std::runtime_error("Output directory already exists; preserve it and choose a fresh run directory");
    std::filesystem::create_directories(out);
    auto start=std::chrono::steady_clock::now();
    std::ofstream stages(out/"stages.tsv");
    if(!stages)throw std::runtime_error("Cannot create stage ledger");
    stages<<std::setprecision(17)<<"stage\tstart_seconds\tend_seconds\telapsed_seconds\n";
    std::string active_stage="reference_and_libraries";double stage_start=0;
    auto next_stage=[&](const std::string&next) {
        const double now=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        stages<<active_stage<<'\t'<<stage_start<<'\t'<<now<<'\t'<<now-stage_start<<'\n';
        stages.flush();if(!stages)throw std::runtime_error("Stage ledger write failed");
        active_stage=next;stage_start=now;
    };
    std::ofstream detailed_log(out/"inference.log");
    if(!detailed_log)throw std::runtime_error("Cannot create inference log");
    hic::ProgressLog progress(std::cerr,detailed_log);
    auto log=[&](const std::string&s){progress.detail(s);};
    progress.stage("start: threads="+std::to_string(threads)+" R="+std::to_string(budget)+" libraries="+std::to_string(libraries.size()));
    std::vector<std::string> sequences(2);
    std::vector<std::pair<uint64_t,uint64_t>> unlinked_intervals;
    const size_t read_linked_blocks=0;
    // Per-variant phase support (stage 2 output rule): canonical word -> summed
    // calibrated orientation log-likelihood ratio (nats) of informative factors.
    std::unordered_map<uint64_t,double> word_support;
    double nats_per_unit=0,site_threshold=0;std::string site_policy="phred30";
    unsigned linkage_factor_denominator=1;
    size_t supported_targets=0;
    size_t unlinked_segments=0,copy_number_segments=0,nondiploid_windows=0;
    std::vector<unsigned char> nondiploid_window;   // per 25 kb reference window
    constexpr uint64_t copy_window=25000;
    { // Release graph, evidence, DP storage and audits before assembly calling.
    log("graph input; objective=three-state-observed-dosage+continuous-paired-and-panel-phase; lambda=1; raw-input-only=yes");
    if(args.count("--gindex")&&!args.count("--REF"))args["--REF"]=graph.paths.at(0).name;
    log("graph loaded; constructing lexicographic k31/w25 evidence");
    size_t reference_id=SIZE_MAX;
    if(args.count("--REF")) {
        std::vector<size_t> found;
        for(size_t h=0;h<graph.paths.size();++h)if(graph.paths[h].name==args["--REF"]){found={h};break;}
        if(found.empty())for(size_t h=0;h<graph.paths.size();++h)
            if(graph.paths[h].name.rfind(args["--REF"]+"#",0)==0)found.push_back(h);
        if(found.size()!=1)throw std::runtime_error("--REF must resolve to exactly one complete panel path");
        reference_id=found[0];
        std::ofstream ref(out/"reference.fa");ref<<'>'<<graph.paths[reference_id].name<<'\n'<<graph.spell(reference_id)<<'\n';
        ref.close();if(!ref)throw std::runtime_error("Reference FASTA write failed");
    }
    if(reference_id==SIZE_MAX)throw std::invalid_argument("Inference requires a designated --REF path");
    {std::ofstream f(out/"resolved-libraries.tsv");f<<"library\tkind\tr1\tr2\n";
     for(const auto&l:libraries)f<<l.name<<'\t'<<l.kind<<'\t'<<l.r1<<'\t'<<l.r2<<'\n';
     f.close();if(!f)throw std::runtime_error("Cannot write resolved read libraries");}
    // Keep reference first before whole-sample panel selection.
    if(reference_id!=SIZE_MAX&&reference_id!=0)std::swap(graph.paths[0],graph.paths[reference_id]);
    next_stage("panel_selection_and_sequence_index");
    progress.stage("selecting panel and indexing sequence evidence");
    hic::select_panel(graph,libraries,out,threads);
    auto parts=hic::CommonSegments::build(graph,31,threads);
    log("graph parsed: panel="+std::to_string(graph.paths.size())+" joint segments="+std::to_string(parts.blocks()));
    auto index=hic::SequenceIndex::build(graph,parts,31,25,threads);
    log("indexed canonical strings="+std::to_string(index.targets.size()));
    progress.stage("index ready: paths="+std::to_string(graph.paths.size())+" segments="+std::to_string(parts.blocks())+" targets="+std::to_string(index.targets.size()));
    progress.stage("counting reads and collecting paired evidence");
    hic::LibraryEvidence evidence;
    next_stage("native_localizer_index");
    progress.stage("building native graph-context localizer");
    hic::NativeLocalizer::Policy policy;policy.k=index.k;policy.w=index.w;
    hic::NativeLocalizer localizer(hic::native_graph_contexts(graph,0),policy,threads);
    const auto&is=localizer.stats();
    {
        std::ofstream f(out/"native-index.tsv");
        f<<"context_bases\t"<<is.context_bases<<"\nselected_words\t"<<is.selected_words
         <<"\nretained_occurrences\t"<<is.retained_occurrences<<"\nphysical_loci\t"<<is.physical_loci
         <<"\ndiscarded_repetitive_words\t"<<is.discarded_words
         <<"\npeak_ingest_hit_scratch_bytes\t"<<is.peak_ingest_scratch_bytes
         <<"\ningest_batches\t"<<is.ingest_batches<<'\n';
        f.close();if(!f)throw std::runtime_error("Native index ledger write failed");
    }
    log("native index: contexts="+std::to_string(localizer.contexts().size())+
        " words="+std::to_string(is.selected_words)+" retained_occurrences="+std::to_string(is.retained_occurrences)+
        " peak_ingest_hit_scratch_bytes="+std::to_string(is.peak_ingest_scratch_bytes));
    next_stage("native_localization_counts_and_audit");
    progress.stage("localizing raw paired reads; no external read mapper");
    hic::NativeObservationAudit audit(out/"native-observations.tsv.gz",threads);
    const hic::NativePairPolicy pair_policy;
    evidence=hic::read_native_library_evidence(index,libraries,localizer,threads,pair_policy,
                 [&](std::string_view row){audit.write(row);},
                 [&](const hic::ReadLibrary&library,uint64_t pairs,uint64_t linked,bool complete) {
                     progress.pairs("reads "+library.kind+"/"+library.name,pairs,linked,complete);
                 });
    audit.close();
    next_stage("model_and_factor_preparation");
    auto&raw=evidence.combined;
    log("raw pairs="+std::to_string(raw.pairs)+" linkage candidates="+std::to_string(raw.contact_pairs)+
        " admitted marker pairs="+std::to_string(raw.markers.size()));
    progress.stage("reads complete: pairs="+std::to_string(raw.pairs)+"; calibrating dosage and phase factors");
    log("target counts=all exact per-end matches; exact_observations="+std::to_string(raw.target_end_observations)+
        " minimizer_selected_observations="+std::to_string(raw.minimizer_target_end_observations)+" phase_markers=unchanged_minimizers");
    {
        std::ofstream f(out/"calibration-histogram.tsv");f<<"read_end_count\tper_library_read_minimizer_types\n";
        for(auto[count,n]:raw.count_histogram)f<<count<<'\t'<<n<<'\n';
        f.close();if(!f)throw std::runtime_error("Calibration histogram write failed");
    }
    constexpr int phase_weight=1;
    // Two-stage phase evidence. Stage 1 (DP, donor pairs and genotypes): every
    // panel-unique marker with equal molecule shares, as in v1. Stage 2
    // (orientation only, after the donor pairs are fixed): genome-unique
    // markers with heterozygosity-proportional shares and calibrated
    // distance confidence. Orientation changes never alter a donor pair.
    const std::string vote_policy="stage1_uniform_stage2_heterozygous";
    log("phase evidence stages=two");
    {
        // Copy number per 25 kb window: nearest integer of 2 x (median count of
        // conserved single-copy strings in the window / genome-wide median).
        // Windows with fewer than 20 controls are treated as diploid.
        bool has_sr=false;for(const auto&l:evidence.counts)has_sr|=l.kind=="SR";
        const auto full=hic::DonorMask::full(graph.paths.size());
        std::vector<std::vector<uint64_t>> window_counts(parts.cuts.at(0).back()/copy_window+1);
        std::vector<uint64_t> all;
        for(size_t x=0;x<index.targets.size();++x) {
            const auto&t=index.targets[x];
            if(!t.unique_locus||t.carriers!=full||t.occurrences.size()!=1)continue;
            uint64_t c=0;for(const auto&l:evidence.counts)if(!has_sr||l.kind=="SR")c+=l.counts[x];
            window_counts[parts.cuts.at(0).at(size_t(t.occurrences.front().first))/copy_window].push_back(c);all.push_back(c);
        }
        auto median=[](std::vector<uint64_t>&v){std::nth_element(v.begin(),v.begin()+v.size()/2,v.end());return v[v.size()/2];};
        if(all.empty())throw std::runtime_error("No conserved controls for copy-number windows");
        const long double global=std::max<uint64_t>(1,median(all));
        nondiploid_window.assign(window_counts.size(),0);
        std::ofstream f(out/"copy-number-windows.tsv");f<<"window_start\tcontrols\tmedian\tcopy_number\n";
        for(size_t w=0;w<window_counts.size();++w) {
            auto&v=window_counts[w];if(v.size()<20)continue;
            const uint64_t m=median(v);const long cn=std::lround(2.L*m/global);
            f<<w*copy_window<<'\t'<<v.size()<<'\t'<<m<<'\t'<<cn<<'\n';
            if(cn!=2){nondiploid_window[w]=1;++nondiploid_windows;}
        }
        f.close();if(!f)throw std::runtime_error("Copy-number window write failed");
        log("copy-number windows: non-diploid="+std::to_string(nondiploid_windows)+" of "+std::to_string(window_counts.size()));
    }
    const std::string dosage_prior="global",dosage_weights="logodds";
    // Clip at phred 30 (1000:1 odds), the constant of the phase claims.
    const long double log_odds_cap=30.L*std::log(10.L)/10.L;
    // Linked observations: markers from all strings; panel relations in both stages.
    const std::string marker_policy="all",panel_policy="on";
    auto prepared=hic::build_model(index,raw,int(graph.paths.size()),int(parts.blocks()),budget,threads,evidence.counts,phase_weight,false,
                                   true,false,log_odds_cap,false,true,false);
    log("linkage markers: "+marker_policy+"; panel relations: "+panel_policy);
    log("dosage prior: "+dosage_prior+"; dosage weights: "+dosage_weights);
    auto&m=prepared.model;
    hic::native_factors::Snapshot native_stage1,native_stage2;
    native_stage1=hic::native_factors::prepare(prepared,index,raw,false,false,threads);
    native_stage2=hic::native_factors::prepare(prepared,index,raw,true,false,threads);
    hic::native_factors::apply(native_stage1,m);
    prepared.informative_molecules=native_stage1.molecules;prepared.raw_molecule_units=native_stage1.raw_units;
    hic::native_factors::normalize(prepared,phase_weight);
    {
        std::ofstream receipt(out/"native-factors.stage1.tsv");hic::native_factors::write_snapshot(receipt,native_stage1);
        receipt.close();if(!receipt)throw std::runtime_error("Native stage-1 receipt close failed");
    }
    // Unordered donor pairs (free haplotype exchange); the all-pair budget
    // repair charges donor changes in order.
    m.free_swaps=true;
    const std::string repair_charge="ordered";
    m.ordered_repair=true;
    // Presence phase factor in both stages.
    const std::string factor_policy="presence";
    m.neutral_phase=false;m.presence_phase=true;
    log("phase factor="+factor_policy+"; stage 2="+factor_policy);
    log("recombination charge="+std::string(m.free_swaps?"unordered donor pair (free haplotype exchange)":"ordered v1"));
    {
        std::ofstream f(out/"dosage-fit.tsv");f<<std::setprecision(20);
        const auto&fit=prepared.observed_dosage;
        f<<"metric\tvalue\nmean\t"<<fit.calibration.mean<<"\nshape\t"<<fit.calibration.shape
         <<"\nerror_rate\t"<<fit.error_rate<<"\nerror_positive_mean\t"<<fit.error_positive_mean()<<"\nprior0\t"<<fit.prior[0]
         <<"\nprior1\t"<<fit.prior[1]<<"\nprior2\t"<<fit.prior[2]
         <<"\nobservations\t"<<fit.observations<<"\niterations\t"<<fit.iterations
         <<"\nlog_likelihood\t"<<fit.log_likelihood<<"\nerror_zero\t"<<fit.error_zero
         <<"\ndropout\t"<<fit.dropout<<"\nzero_observations\t"<<fit.zero_observations<<'\n';
        f.close();if(!f)throw std::runtime_error("Dosage fit audit write failed");
        log("complete dosage fit: targets="+std::to_string(fit.observations)+" zeros="+std::to_string(fit.zero_observations)+" iterations="+std::to_string(fit.iterations)+
            " error_positive_mean="+std::to_string(double(fit.error_positive_mean())));
    }
    log("hom="+std::to_string(prepared.hom)+" het="+std::to_string(prepared.het)+" informative molecules="+std::to_string(prepared.informative_molecules)+" contact factors="+std::to_string(m.contacts.size()));
    {
        std::ofstream f(out/"libraries.tsv");
        f<<"library\tkind\traw_pairs\tphase_candidate_pairs\thom_mean\thom_variance\tNB_shape\tdosage_evidence\tadmitted_marker_pairs\n";
        for(size_t l=0;l<libraries.size();++l) {
            const auto&c=prepared.library_calibrations[l];
            f<<libraries[l].name<<'\t'<<libraries[l].kind<<'\t'<<evidence.pairs[l]<<'\t'
             <<evidence.contact_pairs[l]<<'\t'<<c.mean<<'\t'<<c.variance<<'\t'<<c.shape<<'\t'
             <<(std::find(prepared.dosage_libraries.begin(),prepared.dosage_libraries.end(),l)!=prepared.dosage_libraries.end())
             <<'\t'<<evidence.marker_pairs[l]<<'\n';
        }
        f.close();if(!f)throw std::runtime_error("Library calibration write failed");
        std::ofstream counts(out/"library-counts.tsv");counts<<"target";
        for(const auto&l:libraries)counts<<'\t'<<l.name;
        counts<<'\n';
        for(size_t x=0;x<m.targets.size();++x) {
            counts<<x;
            for(const auto&l:evidence.counts)counts<<'\t'<<l.counts[prepared.target_ids[x]];
            counts<<'\n';
        }
        counts.close();if(!counts)throw std::runtime_error("Library count audit write failed");
    }
    {
        std::ofstream f(out/"inputs.tsv");f<<std::setprecision(18);
        f<<"target\tcanonical_word\tread_end_count\thom\tgraph_unique_locus\tfirst_segment\tlast_segment\treward_units\thom_probability\tother_reward_units\tpanel_copy_histogram\tdosage_read_end_count\tabsent_reward_units\tabsent_probability\tcopy_level\tcopy_classes\tcopy_occurrence_count\tpanel_carriers\n";
        for(size_t x=0;x<m.targets.size();++x){
            auto&t=index.targets[prepared.target_ids[x]];
            f<<x<<'\t'<<t.word<<'\t'<<t.count<<'\t'<<m.targets[x].hom<<'\t'<<t.unique_locus<<'\t'<<t.occurrences.front().first<<'\t'<<t.occurrences.back().first
             <<'\t'<<m.targets[x].reward_units<<'\t'<<prepared.hom_probabilities[x]<<'\t'<<m.targets[x].other_reward_units<<'\t';
            for(size_t i=0;i<t.copy_histogram.size();++i)f<<(i?",":"")<<t.copy_histogram[i].first<<':'<<t.copy_histogram[i].second;
            f<<'\t'<<prepared.dosage_counts[x]<<'\t'<<m.targets[x].absent_reward_units<<'\t'<<prepared.absent_probabilities[x]
             <<'\t'<<prepared.copy_level[x]<<'\t';
            if(prepared.copy_level[x]) {
                const uint64_t zero=uint64_t(m.H)-t.carriers.count();
                f<<"0:"<<zero;
                for(auto[c,n]:t.copy_histogram)f<<','<<c<<':'<<n;
            } else f<<'.';
            f<<'\t'<<prepared.occurrence_counts[x]<<'\t'<<t.carriers.count()<<'\n';
        }
        f.close();if(!f)throw std::runtime_error("Input audit write failed");
    }
    {
        // These are the actual fixed factors entering the joint objective,
        // not contacts selected afterwards from the returned haplotypes.
        std::ofstream f(out/"contacts.tsv");f<<"target1\ttarget2\tweight_units\trelation\tsource\n";
        for(const auto&e:m.contacts)f<<e.x<<'\t'<<e.y<<'\t'<<e.units<<'\t'<<e.relation<<'\t'<<(e.panel?"panel":"molecule")<<'\n';
        f.close();if(!f)throw std::runtime_error("Contact audit write failed");
    }
    // Reference coordinate of each target's first segment (linkage distance).
    // word -> (target, copy threshold). Ordinary targets need one copy;
    // copy-level targets need at least `level` copies in the emitted haplotype.
    std::unordered_map<uint64_t,std::vector<std::pair<size_t,uint32_t>>> audit_targets;
    for(size_t x=0;x<prepared.target_ids.size();++x)
        audit_targets[index.targets[prepared.target_ids[x]].word].emplace_back(x,std::max<uint32_t>(1,prepared.copy_level[x]));
    // The old molecule/marker storage is no longer needed by the joint DP.
    raw.markers.clear();raw.markers.shrink_to_fit();index.ids.clear();index.ids.rehash(0);index.targets.clear();index.targets.shrink_to_fit();
    evidence.counts.clear();evidence.counts.shrink_to_fit();
    prepared.linkage_remap.clear();prepared.linkage_remap.shrink_to_fit();
    const size_t search_width=512;
    next_stage("joint_dp");
    auto reduced=hic::fold_constant_presence(m);
    log("constant reward reduction: targets="+std::to_string(reduced.fixed_targets)+" contacts="+std::to_string(reduced.fixed_contacts)+" exact_offset="+hic::score_text(reduced.offset));
    log("joint search: bidirectional decision order; resource-aware and shared-endpoint phase bounds; R="+std::to_string(budget)+"; total_width="+std::to_string(search_width)+"; constants retained; objective unchanged");
    progress.stage("joint DP: segments="+std::to_string(m.blocks)+" histories="+std::to_string(search_width));
    auto joint=hic::solve_presence_bidirectional(reduced.model,search_width,threads,[&](const char*direction,int b,const hic::Result&r){
        log(std::string("DP ")+direction+" segments="+std::to_string(b)+"/"+std::to_string(m.blocks)+" peak_live="+std::to_string(r.peak)+" contact_frontier="+std::to_string(r.max_frontier)+" transitions="+std::to_string(r.generated));
        progress.progress(std::string("DP ")+direction,uint64_t(b),uint64_t(m.blocks));
    });
    if(joint.forward_score!=hic::impossible)joint.forward_score+=reduced.offset;
    if(joint.reverse_score!=hic::impossible)joint.reverse_score+=reduced.offset;
    if(!joint.best.path.empty()){joint.best.score+=reduced.offset;joint.best.upper+=reduced.offset;}
    const auto fixed_targets=reduced.fixed_targets,fixed_contacts=reduced.fixed_contacts;
    const auto fixed_reward=reduced.offset;
    reduced.model=hic::PresenceModel{};
    auto result=std::move(joint.best);
    if(result.path.empty())throw std::runtime_error("No returned path pair; approximation_ratio=NULL");
    const auto initial_joint_score=result.score;
    next_stage("global_phase_and_budget_repair");
    size_t global_moves_accepted=0;
    progress.stage("joint DP complete; global phase proposals and budget repair");
    {
        std::ofstream audit(out/"global-phase-moves.tsv");
        audit<<"maximum_segment_span\tblocks\tedges\tbefore\tunconstrained_score\tunconstrained_R\trepaired_score\trepaired_R\tafter\taccepted\ttrace_bytes\n";
        for(size_t span:{std::numeric_limits<size_t>::max(),size_t(256)}) {
            log("global phase proposal and budget repair; max_segment_span="+std::to_string(span));
            const auto r=hic::global_phase::move(m,result,span,threads);
            global_moves_accepted+=r.accepted;
            audit<<r.span<<'\t'<<r.blocks<<'\t'<<r.edges<<'\t'<<r.before<<'\t'<<r.unconstrained<<'\t'
                 <<r.proposed_recombinations<<'\t'<<r.repaired<<'\t'<<r.repaired_recombinations<<'\t'
                 <<r.after<<'\t'<<r.accepted<<'\t'<<r.trace_bytes<<'\n';
            audit.flush();if(!audit)throw std::runtime_error("Global phase audit write failed");
            log("global move complete; accepted="+std::to_string(r.accepted)+" score="+hic::score_text(r.after)
                +" proposed_R="+std::to_string(r.proposed_recombinations)+" repaired_R="+std::to_string(r.repaired_recombinations));
        }
    }
    next_stage("calibration_and_fixed_genotype_orientation");
    const std::string calibration_policy="distance",calibration_reference="reads";
    hic::Score segment_upper_override=-1;
    {
        progress.stage("phase stage 2: fixed-genotype orientation");
        // Replace stage-1 molecule factors by stage-2 factors; panel relations kept.
        {
            std::ofstream f(out/"contacts.stage1.tsv");f<<"target1\ttarget2\tweight_units\trelation\tsource\n";
            for(const auto&e:m.contacts)f<<e.x<<'\t'<<e.y<<'\t'<<e.units<<'\t'<<e.relation<<'\t'<<(e.panel?"panel":"molecule")<<'\n';
            f.close();if(!f)throw std::runtime_error("Stage-1 contact audit write failed");
        }
        hic::native_factors::apply(native_stage2,m);
        prepared.informative_molecules=native_stage2.molecules;prepared.raw_molecule_units=native_stage2.raw_units;
        native_stage1={};
        prepared.read_contact_units=0;prepared.panel.units=0;
        for(const auto&e:m.contacts)(e.panel?prepared.panel.units:prepared.read_contact_units)+=e.units;
        prepared.contact_units=prepared.read_contact_units+prepared.panel.units;
        m.weights=hic::normalized_weights(prepared.contact_units,prepared.het_reward_units,phase_weight);m.check();
        result.score=m.score(result.path);result.upper=m.universal_bound();segment_upper_override=result.upper;
        log("stage-2 phase factors: molecules="+std::to_string(prepared.informative_molecules)+" factors="+std::to_string(m.contacts.size())
            +" rescored="+hic::score_text(result.score));
    }
    {
        // Two-stage (HapCUT2-style) h-trans calibration: eps per distance class from
        // the stage-1 phased pair, fixed thereafter; orientation is re-optimized.
        // Non-circular reference: the stage-1 pair was oriented with panel
        // relations, so calibrating them against it confirms them by
        // construction. Every class, panel included, is therefore calibrated
        // against the orientation that the read factors alone support
        // (exchanges only; donor pairs and genotypes unchanged).
        std::vector<int> reference_path=result.path;
        {
            hic::PresenceModel reads=m;
            reads.contacts.erase(std::remove_if(reads.contacts.begin(),reads.contacts.end(),[](const auto&e){return e.panel;}),reads.contacts.end());
            hic::Result oriented;oriented.path=result.path;oriented.score=reads.score(oriented.path);oriented.upper=reads.universal_bound();
            const auto report=hic::free_orientation::improve(reads,oriented);
            reference_path=oriented.path;
            log("calibration reference: read-only orientation, accepted flips="+std::to_string(report.accepted));
        }
        {
            const auto&snapshot=native_stage2;
            const auto estimate=hic::native_factors::estimate(snapshot,m,reference_path);
            std::ofstream factors(out/"native-factors.uncalibrated.tsv");hic::native_factors::write_snapshot(factors,snapshot);
            factors.close();if(!factors)throw std::runtime_error("Native factor receipt close failed");
            std::ofstream f(out/"linkage-calibration.tsv");
            f<<std::setprecision(20)<<"class\tagree_units\toppose_units\teps\tllr\tconfidence\n";
            for(const auto&[identity,c]:estimate.classes)f<<hic::native_factors::name(identity)<<'\t'
                <<c.agree<<'\t'<<c.oppose<<'\t'<<c.epsilon<<'\t'<<c.llr<<'\t'<<c.confidence<<'\n';
            f.close();if(!f)throw std::runtime_error("Native calibration receipt write failed");
            std::ofstream reference(out/"native-calibration-presence.tsv");reference<<"target\tpresence\n";
            for(size_t x=0;x<m.targets.size();++x) {
                unsigned char presence=0;for(const auto&[b,mask]:m.targets[x].occurrences)presence|=m.signature(mask,reference_path.at(size_t(b)));
                reference<<x<<'\t'<<unsigned(presence)<<'\n';
            }
            reference.close();if(!reference)throw std::runtime_error("Native calibration reference write failed");
            hic::native_factors::apply(snapshot,m,&estimate);
            hic::native_factors::normalize(prepared,phase_weight);
            result.score=m.score(result.path);result.upper=m.universal_bound();segment_upper_override=result.upper;
            linkage_factor_denominator=hic::native_factors::scale;
            nats_per_unit=double(estimate.maximum_llr/(static_cast<long double>(hic::native_factors::scale)*linkage_factor_denominator));
            log("native linkage calibration: classes="+std::to_string(estimate.classes.size())+" rescored="+hic::score_text(result.score));
        }
    }
    // Persist the final factors.
    {
        std::ofstream f(out/"contacts.tsv");f<<"target1\ttarget2\tweight_units\trelation\tsource\n";
        for(const auto&e:m.contacts)f<<e.x<<'\t'<<e.y<<'\t'<<e.units<<'\t'<<e.relation<<'\t'<<(e.panel?"panel":"molecule")<<'\n';
        f.close();if(!f)throw std::runtime_error("Final native factors write failed");
        native_stage1={};native_stage2={};
    }
    // Donor pairs (and hence genotypes) are fixed from here: only zero-charge
    // orientation exchanges follow.
    std::vector<int> fixed_pairs(result.path.size());
    for(size_t b=0;b<result.path.size();++b){const int z=result.path[b];fixed_pairs[b]=std::min(z/m.H,z%m.H)*m.H+std::max(z/m.H,z%m.H);}
    hic::free_orientation::Report orientation;
    {
        log("exact interval orientation search (zero-charge exchanges)");
        orientation=hic::free_orientation::improve(m,result);
        std::ofstream audit(out/"orientation-search.tsv");
        audit<<"rounds\tcandidates\tevaluated\taccepted\tcontinuity_flips\texchanges_before\texchanges_after\tflipped_segments\tbefore\tafter\n"
             <<orientation.rounds<<'\t'<<orientation.candidates<<'\t'<<orientation.evaluated<<'\t'<<orientation.accepted<<'\t'
             <<orientation.continuity_flips<<'\t'<<orientation.exchanges_before<<'\t'<<orientation.exchanges_after<<'\t'
             <<orientation.flipped_segments<<'\t'<<orientation.before<<'\t'<<orientation.after<<'\n';
        audit.close();if(!audit)throw std::runtime_error("Orientation audit write failed");
        log("orientation search complete; accepted="+std::to_string(orientation.accepted)+" of "+std::to_string(orientation.evaluated)
            +" continuity_flips="+std::to_string(orientation.continuity_flips)+" exchanges "+std::to_string(orientation.exchanges_before)
            +"->"+std::to_string(orientation.exchanges_after)
            +" score="+hic::score_text(orientation.after));
    }
    for(size_t b=0;b<result.path.size();++b){const int z=result.path[b];
        if(fixed_pairs[b]!=std::min(z/m.H,z%m.H)*m.H+std::max(z/m.H,z%m.H))throw std::logic_error("Orientation stage changed a donor pair");}
    next_stage("phase_reporting_reconstruction_and_certificates");
    // Phase threshold: phred 30 (1000:1 odds), chosen empirically; HapCUT2's
    // default pruning threshold is a posterior of 0.8 (about phred 7).
    site_threshold=30.0*std::log(10.0)/10.0;
    if(nats_per_unit>0) {
        std::vector<unsigned char> presence(m.targets.size());
        for(size_t x=0;x<m.targets.size();++x)for(const auto&[b,mask]:m.targets[x].occurrences)presence[x]|=m.signature(mask,result.path[b]);
        std::vector<double> support(m.targets.size());
        auto single=[&](size_t x){return m.targets[x].occurrences.size()==1;};
        for(const auto&e:m.contacts) {
            if(!single(e.x)||!single(e.y)||m.targets[e.x].occurrences[0].first==m.targets[e.y].occurrences[0].first)continue;
            const int s=e.relation*hic::phase_orientation(presence[e.x])*hic::phase_orientation(presence[e.y]);
            if(!s)continue;
            const double v=double(e.units)*s*nats_per_unit;support[e.x]+=v;support[e.y]+=v;
        }
        for(const auto&[w,list]:audit_targets) {
            double total=0;bool any=false;
            for(const auto&[x,need]:list)if(support[x]!=0){total+=support[x];any=true;}
            if(any){word_support.emplace(w,total);++supported_targets;}
        }
        log("per-variant phase support: policy="+site_policy+" threshold_nats="+std::to_string(site_threshold)
            +" nats_per_unit="+std::to_string(nats_per_unit)+" supported_words="+std::to_string(supported_targets));
    }
    const std::string unlinked_policy="separate";
    const size_t confidence_splits=0;
    {
        // A segment whose heterozygous strings share no informative factor with
        // any other segment has an orientation fixed by no evidence (only the
        // continuity tie-break); a segment in a non-diploid copy-number window
        // violates the diploid model its linked evidence assumes. Each such
        // segment is reported as its own phase set: its internal phase (the
        // donor sequences) is kept, its orientation relative to the rest is not
        // claimed, as phasers split blocks where linking fragments are absent.
        std::vector<unsigned char> presence(m.targets.size());
        for(size_t x=0;x<m.targets.size();++x)for(const auto&[b,mask]:m.targets[x].occurrences)presence[x]|=m.signature(mask,result.path[b]);
        auto het=[&](size_t x){return (presence[x]==1||presence[x]==2)&&m.targets[x].occurrences.size()==1;};
        std::vector<unsigned char> has(m.blocks),linked(m.blocks);
        for(size_t x=0;x<m.targets.size();++x)if(het(x))has[m.targets[x].occurrences.front().first]=1;
        for(const auto&e:m.contacts) {
            if(!het(e.x)||!het(e.y))continue;
            const int a=m.targets[e.x].occurrences.front().first,b=m.targets[e.y].occurrences.front().first;
            if(a!=b){linked[a]=1;linked[b]=1;}
        }
        std::ofstream bed(out/"unphased-segments.bed");
        for(int b=0;b<m.blocks;++b) {
            if(!has[b])continue;
            const uint64_t lo=parts.cuts.at(0).at(size_t(b)),hi=parts.cuts.at(0).at(size_t(b)+1);
            const bool cn=nondiploid_window.at(lo/copy_window)!=0;
            if(linked[b]&&!cn)continue;
            unlinked_intervals.push_back({lo,hi});
            if(!linked[b])++unlinked_segments;else ++copy_number_segments;
        }
        for(const auto&[lo,hi]:unlinked_intervals)bed<<graph.paths.at(0).name<<'\t'<<lo<<'\t'<<hi<<'\n';
        bed.close();if(!bed)throw std::runtime_error("Unphased segment BED write failed");
        log("separate phase sets: unlinked segments="+std::to_string(unlinked_segments)+" non-diploid-window segments="+std::to_string(copy_number_segments));
    }
    std::vector<size_t> donor[2];unsigned recombinations=0;
    for(size_t b=0;b<result.path.size();++b){int z=result.path[b];donor[0].push_back(z/m.H);donor[1].push_back(z%m.H);if(b)recombinations+=m.charge(result.path[b-1],z);}
    if(recombinations>unsigned(budget))throw std::logic_error("Traceback exceeds shared budget");
    std::ofstream fasta(out/"all-components.fa");
    std::vector<unsigned char> literal_presence(m.targets.size());
    for(int h=0;h<2;++h) {
        auto& sequence=sequences[h];sequence=parts.mosaic(donor[h]);fasta<<'>'<<sample<<".hap"<<h+1<<'\n';
        std::unordered_map<uint64_t,uint32_t> emitted_copies;
        hic::words(sequence,31,[&](uint64_t w,size_t,bool){if(audit_targets.count(w))++emitted_copies[w];});
        for(const auto&[w,n]:emitted_copies)for(const auto&[x,need]:audit_targets.at(w))if(n>=need)literal_presence[x]|=1<<h;
        for(size_t p=0;p<sequence.size();p+=80)fasta<<sequence.substr(p,80)<<'\n';
    }
    fasta.close();if(!fasta)throw std::runtime_error("FASTA write failed");
    hic::Score literal_score=0;
    size_t hom_satisfied=0,het_satisfied=0;
    for(size_t x=0;x<m.targets.size();++x) {
        unsigned char expected=0;for(const auto&[b,mask]:m.targets[x].occurrences)expected|=m.signature(mask,result.path[b]);
        if(literal_presence[x]!=expected)throw std::logic_error("Emitted FASTA target presence differs from DP");
        const auto covered=hic::PresenceModel::coverage(m.targets[x].hom,literal_presence[x]);
        literal_score+=hic::phase_reward_scale*m.weights.genotype*hic::PresenceModel::genotype_reward(m.targets[x],literal_presence[x]);
        if(covered){if(m.targets[x].hom)++hom_satisfied;else ++het_satisfied;}
    }
    const auto genotype_reward=literal_score;
    hic::Score agreeing_units=0,opposite_units=0,nonheterozygous_units=0;
    for(auto&e:m.contacts) {
        const auto a=literal_presence[e.x],b=literal_presence[e.y];
        if(hic::PresenceModel::agreement(a,b,e.relation))agreeing_units+=e.units;
        else if((a==1||a==2)&&(b==1||b==2))opposite_units+=e.units;
        else nonheterozygous_units+=e.units;
        literal_score+=e.units*m.weights.contact*hic::phase_factor(a,b,e.relation,e.panel,m.neutral_phase,m.presence_phase);
    }
    if(agreeing_units+opposite_units+nonheterozygous_units!=prepared.contact_units)
        throw std::logic_error("Contact outcome masses do not sum to fixed input mass");
    if(literal_score!=result.score)throw std::logic_error("Emitted FASTA reward differs from DP");
    log("independent emitted-DNA objective replay=PASS");
    {
        std::ofstream f(out/"presence.tsv");f<<"target\thaplotype1\thaplotype2\n";
        for(size_t x=0;x<literal_presence.size();++x)
            f<<x<<'\t'<<unsigned(literal_presence[x]&1)<<'\t'<<unsigned((literal_presence[x]>>1)&1)<<'\n';
        f.close();if(!f)throw std::runtime_error("Presence audit write failed");
    }
    std::ofstream trace(out/"pair.trace.tsv");trace<<"segment\thaplotype1_donor\thaplotype2_donor\n";
    for(size_t b=0;b<donor[0].size();++b)trace<<b<<'\t'<<graph.paths[donor[0][b]].name<<'\t'<<graph.paths[donor[1][b]].name<<'\n';
    trace.close();if(!trace)throw std::runtime_error("Trace write failed");
    // Restricted-domain certificate and FULL-graph certificate are separate.
    // The latter must cover recombinations omitted by common-anchor segments.
    // After calibration the DP's restricted bound refers to the stage-1 weights;
    // only the full-graph certificate under the final factors is reported.
    const auto segment_upper=segment_upper_override>=0?segment_upper_override:result.upper;result.upper=m.universal_bound();
    std::ofstream metrics(out/"run.tsv");
    metrics<<"metric\tvalue\nraw_pairs\t"<<raw.pairs<<"\ncontact_pairs\t"<<raw.contact_pairs
           <<"\ncontact_pairs_definition\tcandidates_before_linkage_admission"
           <<"\ntarget_count_policy\tall_exact_panel_targets_once_per_read_end"
           <<"\ntarget_end_observations\t"<<raw.target_end_observations
           <<"\nminimizer_target_end_observations\t"<<raw.minimizer_target_end_observations
           <<"\nhom_targets\t"<<prepared.hom<<"\nhet_targets\t"<<prepared.het
           <<"\ninformative_contact_molecules\t"<<prepared.informative_molecules<<"\ncontact_weight_units\t"<<prepared.contact_units
           <<"\nread_contact_weight_units\t"<<prepared.read_contact_units
           <<"\nraw_molecule_weight_units\t"<<prepared.raw_molecule_units
           <<"\npanel_weight_units\t"<<prepared.panel.units<<"\npanel_phase_links\t"<<prepared.panel.links
           <<"\npanel_supported_boundaries\t"<<prepared.panel.boundaries
           <<"\ncalibration_model\tsingle_segment_panel_copy_NB_cross_segment_presence"
           <<"\npresence_only_targets\t"<<prepared.presence_only_targets
           <<"\npanel_selection_policy\tlocal_window_read_consistency_25000bp"
           <<"\ncopy_number_policy\tsingle_segment_copy_levels_marginal_posterior"
           <<"\nphase_weight_diagnostic\t"<<phase_weight
           <<"\nrecombination_charge\t"<<(m.free_swaps?"unordered_pair":"ordered")
           <<"\nphase_factor\t"<<factor_policy<<"\nphase_factor_stage2\t"<<factor_policy
           <<"\norientation_flips_accepted\t"<<orientation.accepted
           <<"\nphase_vote_allocation\t"<<vote_policy
           <<"\nlinkage_markers\t"<<marker_policy<<"\npanel_relations\t"<<panel_policy
           <<"\nunlinked_phase_policy\t"<<unlinked_policy<<"\nread_linked_blocks\t"<<read_linked_blocks<<"\nconfidence_splits\t"<<confidence_splits<<"\nunlinked_segments\t"<<unlinked_segments
           <<"\nsite_phase_support\t"<<site_policy<<"\nsite_phase_threshold_nats\t"<<site_threshold
           <<"\nnondiploid_windows\t"<<nondiploid_windows<<"\ncopy_number_segments\t"<<copy_number_segments
           <<"\nevidence_mode\tnative\nhic_split\tfive-prime"
           <<"\nlinkage_class_version\t"<<hic::native_factors::version
           <<"\nrepair_charge\t"<<repair_charge<<"\nlinkage_calibration\t"<<calibration_policy<<"\ncalibration_reference\t"<<calibration_reference
           <<"\nlinkage_confidence_scale\t"<<hic::native_factors::scale
           <<"\nlinkage_pseudocount_units\t512"
           <<"\nlinkage_factor_denominator\t"<<linkage_factor_denominator
           <<"\ncontinuity_flips\t"<<orientation.continuity_flips<<"\nexchanges_after\t"<<orientation.exchanges_after
           <<"\ncopy_strings\t"<<prepared.copy_strings
           <<"\ncopy_level_targets\t"<<prepared.copy_level_targets
           <<"\nobjective_model\tzero_aware_dosage_all_library_phase"
           <<"\nerror_aware_targets\t"<<prepared.error_aware_targets
           <<"\nabsent_reward_units\t"<<prepared.absent_reward_units
           <<"\ndosage_library_policy\t"<<(prepared.sr_dosage_only?"SR_only":"all_available")
           <<"\nhic_only_targets_excluded\t"<<prepared.hic_only_targets_excluded
           <<"\nzero_dosage_targets\t"<<prepared.zero_dosage_targets
           <<"\nhic_only_targets_retained\t"<<prepared.hic_only_targets_retained
           <<"\ndosage_prior\t"<<dosage_prior<<"\ndosage_weights\t"<<dosage_weights<<"\nlog_odds_cap\t"<<double(log_odds_cap)<<"\npanel_paths\t"<<m.H
           <<"\nobjective_integer_scale\t"<<hic::phase_reward_scale
           <<"\nconfidence_scale\t"<<hic::PreparedModel::confidence_scale
           <<"\nhom_reward_units\t"<<prepared.hom_reward_units<<"\nhet_reward_units\t"<<prepared.het_reward_units
           <<"\nread_libraries\t"<<libraries.size()
           <<"\ncalibration_controls\t"<<prepared.calibration_controls
           <<"\ncalibrated_hom_mean\t"<<prepared.calibration.mean
           <<"\ncalibrated_hom_sd\t"<<std::sqrt(prepared.calibration.variance)
           <<"\ncalibrated_NB_shape\t"<<prepared.calibration.shape
           <<"\nrecombinations\t"<<recombinations
           <<"\nsearch_width\t"<<search_width
           <<"\nglobal_moves_accepted\t"<<global_moves_accepted
           <<"\ninitial_joint_reward_integer\t"<<initial_joint_score
           <<"\nglobal_move_policy\tcontact_orientation_and_resource_repair_v1"
           <<"\nsearch_policy\tbidirectional_strict_global_top_constant_folding"
           <<"\nsearch_completion_bound\tresource_and_disjoint_shared_endpoint_contacts_v1"
           <<"\nconstant_targets\t"<<fixed_targets<<"\nconstant_contacts\t"<<fixed_contacts
           <<"\nconstant_reward_integer\t"<<fixed_reward
           <<"\nsearch_passes\t"<<joint.passes
           <<"\nforward_reward_integer\t"<<joint.forward_score
           <<"\nreverse_reward_integer\t"<<(joint.reverse_score==hic::impossible?"NULL":hic::score_text(joint.reverse_score))
           <<"\nselected_direction\t"<<(joint.selected_reverse?"reverse":"forward")
           <<"\nsearch_capacity_bound\t"<<search_width
           <<"\npeak_joint_states\t"<<result.peak
           <<"\nsearch_completion_penalty\t"<<result.completion_penalty
           <<"\nsearch_relaxed_switches\t"<<result.completion_relaxed_switches
           <<"\nreward_integer\t"<<result.score<<"\nsegment_upper_integer\t"<<segment_upper
           <<"\ngenotype_reward_integer\t"<<genotype_reward<<"\ncontact_reward_integer\t"<<literal_score-genotype_reward
           <<"\nhom_targets_satisfied\t"<<hom_satisfied<<"\nhet_targets_satisfied\t"<<het_satisfied
           <<"\ncontact_agreeing_units\t"<<agreeing_units<<"\ncontact_opposite_units\t"<<opposite_units
           <<"\ncontact_nonheterozygous_units\t"<<nonheterozygous_units
           <<"\nfull_graph_upper_integer\t"<<result.upper<<"\napproximation_ratio\t"<<result.ratio_text()
           // Split certificate: U_G/G and U_C/C bound each part's gap separately;
           // rho <= max(rho_G, rho_C). Exchanges leave G and hence rho_G unchanged.
           <<"\ngenotype_upper_integer\t"<<m.genotype_bound()<<"\nphase_upper_integer\t"<<m.phase_bound()
           <<"\ngenotype_approximation_ratio\t"<<hic::Result{{0},genotype_reward,m.genotype_bound()}.ratio_text()
           <<"\nphase_approximation_ratio\t"<<hic::Result{{0},literal_score-genotype_reward,m.phase_bound()}.ratio_text()
           <<"\ngenotype_weight_integer\t"<<m.weights.genotype<<"\ncontact_weight_integer\t"<<m.weights.contact
           <<"\nseconds\t"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    metrics.close();if(!metrics)throw std::runtime_error("Metric write failed");
    log("FASTA written; reward="+hic::score_text(result.score)+" full_graph_upper="+hic::score_text(result.upper)+" approximation_ratio="+result.ratio_text());
    {
        // Printed certificate: rho = U/L and its coverage/phase parts (rho <= max).
        const hic::Result g{{0},genotype_reward,m.genotype_bound()},c{{0},literal_score-genotype_reward,m.phase_bound()};
        progress.stage("FASTA ready: rho="+result.ratio_text()+" coverage="+g.ratio_text()+" phase="+c.ratio_text()+" (objective bounds)");
        std::cout<<"approximation_ratio\t"<<result.ratio_text()<<"\ngenotype_approximation_ratio\t"<<g.ratio_text()
                 <<"\nphase_approximation_ratio\t"<<c.ratio_text()<<'\n';
    }
    }
    next_stage(args.count("--output-vcf")?"assembly_calling":"finalization");
    if(args.count("--output-vcf")) {
        progress.stage("calling phased variants with dipcall");
        const auto before=std::chrono::steady_clock::now();
        dipcall_vcf::Config caller{(out/"reference.fa").string(),args["--output-vcf"],sample,unsigned(threads),{},{}};
        caller.separate=unlinked_intervals;
        if(site_policy!="off"&&nats_per_unit>0) {
            if(caller.separate.empty())caller.separate.push_back({0,0});
            const bool strict=site_policy=="phred30";const double threshold=site_threshold;
            caller.site_supported=[&,strict,threshold](const std::array<std::string,2>&haps,const std::array<std::pair<size_t,size_t>,2>&span) {
                // Distinct strings covering the variant allele on either haplotype.
                std::unordered_map<uint64_t,double> seen;
                for(int h=0;h<2;++h) {
                    const auto&s=haps[h];const size_t K=31;
                    const size_t b=span[h].first,e=std::max(span[h].second,b+1);
                    if(s.size()<K)continue;
                    const size_t lo=b+1>=K?b+1-K:0,hi=std::min(e-1,s.size()-K);
                    if(lo>hi)continue;
                    hic::words(std::string_view(s).substr(lo,hi-lo+K),K,[&](uint64_t w,size_t,bool){
                        auto it=word_support.find(w);if(it!=word_support.end())seen.emplace(w,it->second);});
                }
                double total=0;for(const auto&[w,v]:seen){(void)w;total+=v;}
                return strict?total>=threshold:total>0;
            };
        }
        dipcall_vcf::run(caller,sequences);
        std::ofstream timing(out/"caller.tsv");timing<<"caller\tseconds\ndipcall\t"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count()<<'\n';
        timing.close();if(!timing)throw std::runtime_error("Caller timing write failed");
        log("dipcall VCF complete");
        progress.stage("VCF ready: "+args.at("--output-vcf"));
    }
    next_stage("");stages.close();if(!stages)throw std::runtime_error("Stage ledger close failed");
    progress.stage("complete; diagnostics: "+(out/"inference.log").string());
    return 0;
}
