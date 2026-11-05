#pragma once
#include "ReadLibraries.hpp"
#include "WindowSelection.hpp"

#include <filesystem>
#include <iomanip>
#include <queue>
#include <iostream>
namespace hic {
inline uint64_t sketch_hash(uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
// Bottom-K hash sketch with exact counts for every surviving identity.
// The monotonically shrinking threshold never admits a previously rejected
// identity. Counts in the final sketch therefore include its first observation.
class ReadSketch {
    size_t capacity_;std::priority_queue<std::pair<uint64_t,uint64_t>> largest;
public:
    std::unordered_map<uint64_t,uint64_t> counts;
    explicit ReadSketch(size_t n):capacity_(n){if(!n)throw std::invalid_argument("empty sketch capacity");}
    std::pair<uint64_t,uint64_t> threshold()const {
        return counts.size()==capacity_?largest.top():std::make_pair(UINT64_MAX,UINT64_MAX);
    }
    void add(uint64_t word,uint64_t n=1){
        const auto rank=std::make_pair(sketch_hash(word),word);
        if(rank>threshold())return;
        auto found=counts.find(word);if(found!=counts.end()){if(n>UINT64_MAX-found->second)throw std::overflow_error("Sketch count");found->second+=n;return;}
        if(counts.size()==capacity_){if(rank>=largest.top())return;counts.erase(largest.top().second);largest.pop();}
        counts.emplace(word,n);largest.push(rank);
    }
};
// Former global bottom-K sketch ranking. Retained as a tested baseline; it is
// no longer the inference preselection (see select_panel below).
inline void sketch_select_panel(NativeGraph&g,const std::vector<ReadLibrary>&libraries,
                         const std::filesystem::path&out,int threads,size_t cap=51){
    if(g.paths.size()<=cap)return;
    if(cap<3||cap%2!=1||g.paths.size()%2!=1)throw std::runtime_error("Panel selection requires REF plus complete sample pairs");
    ReadSketch sketch(262144);uint64_t raw_pairs=0;RowExecutor pool(threads);
    for(const auto&library:libraries){FastqInput input[2]={FastqInput(library.r1),FastqInput(library.r2)};
        for(;;) {
            std::vector<FastqRead> batch[2];
            pool.for_each(2,[&](size_t side) {
                FastqRead read;batch[side].reserve(4096);
                while(batch[side].size()<4096&&input[side].next(read))batch[side].push_back(std::move(read));
            });
            if(batch[0].size()!=batch[1].size())throw std::runtime_error("Panel ranking mate count mismatch");
            if(batch[0].empty())break;
            const auto threshold=sketch.threshold();
            const size_t count=batch[0].size();std::vector<std::vector<uint64_t>> selected((count+63)/64);
            pool.for_each(selected.size(),[&](size_t task) {
                for(size_t i=task*64;i<std::min(count,(task+1)*64);++i) {
                    (void)paired_molecule_id(batch[0][i],batch[1][i]);
                    for(size_t side=0;side<2;++side) {
                        std::vector<uint64_t> seen;
                        minimizers(batch[side][i].sequence,31,25,[&](uint64_t word,size_t) {
                            if(std::make_pair(sketch_hash(word),word)<=threshold)seen.push_back(word);
                        });
                        std::sort(seen.begin(),seen.end());seen.erase(std::unique(seen.begin(),seen.end()),seen.end());
                        selected[task].insert(selected[task].end(),seen.begin(),seen.end());
                    }
                }
                std::sort(selected[task].begin(),selected[task].end());
            });
            for(const auto&v:selected)for(size_t a=0;a<v.size();) {
                size_t b=a+1;while(b<v.size()&&v[b]==v[a])++b;
                sketch.add(v[a],b-a);a=b;
            }
            const auto before=raw_pairs;raw_pairs+=count;
            if(raw_pairs/1000000!=before/1000000)std::cerr<<"[M::PanelSelection] raw_pairs="<<raw_pairs<<'\n';
        }
    }
    std::vector<uint64_t> words;for(const auto&entry:sketch.counts)words.push_back(entry.first);std::sort(words.begin(),words.end());
    std::unordered_map<uint64_t,size_t> ids;for(size_t i=0;i<words.size();++i)ids.emplace(words[i],i);
    const size_t H=g.paths.size(),N=words.size(),packed=(N+63)/64;
    std::vector<std::vector<uint64_t>> present(H,std::vector<uint64_t>(packed));
    const auto&lookup=ids;
    pool.for_each(H,[&](size_t h){const auto sequence=g.spell(h);minimizers(sequence,31,25,[&](uint64_t word,size_t){auto it=lookup.find(word);if(it!=lookup.end())present[h][it->second/64]|=uint64_t(1)<<(it->second%64);});});
    std::vector<double> weights(N),scores(H);
    const auto&counts=sketch.counts;
    pool.for_each((N+255)/256,[&](size_t block) {
        for(size_t x=block*256;x<std::min(N,(block+1)*256);++x) {
            size_t carriers=0;for(const auto&row:present)carriers+=(row[x/64]>>(x%64))&1;
            if(carriers&&!(present[0][x/64]>>(x%64)&1))weights[x]=std::log1p(double(counts.at(words[x])))*std::log(double(H+1)/double(carriers+1));
        }
    });
    pool.for_each(H-1,[&](size_t donor){const size_t h=donor+1;
        for(size_t x=0;x<N;++x)if(present[h][x/64]>>(x%64)&1)scores[h]+=weights[x];
    });
    auto sample=[](const std::string&name){auto p=name.find('#');return p==std::string::npos?name:name.substr(0,p);};
    std::map<std::string,std::vector<size_t>> groups;for(size_t h=1;h<H;++h)groups[sample(g.paths[h].name)].push_back(h);
    struct Pair{std::string sample;size_t a,b;double score;};std::vector<Pair> rank;
    for(const auto&[name,paths]:groups){if(paths.size()!=2)throw std::runtime_error("Cannot resolve two panel haplotypes for "+name);size_t a=paths[0],b=paths[1];rank.push_back({name,a,b,std::max(scores[a],scores[b])+.5*std::min(scores[a],scores[b])});}
    std::sort(rank.begin(),rank.end(),[](const auto&a,const auto&b){return a.score!=b.score?a.score>b.score:a.sample<b.sample;});
    std::vector<size_t> chosen{0};for(size_t i=0;i<(cap-1)/2;++i){chosen.push_back(rank.at(i).a);chosen.push_back(rank.at(i).b);}std::sort(chosen.begin(),chosen.end());
    std::ofstream report(out/"panel-selection.tsv");report<<"path\tscore\tretained\n";
    for(size_t h=0;h<H;++h)report<<g.paths[h].name<<'\t'<<std::setprecision(17)<<scores[h]<<'\t'<<std::binary_search(chosen.begin(),chosen.end(),h)<<'\n';
    report.close();if(!report)throw std::runtime_error("Cannot write panel selection audit");
    std::vector<NativeGraph::Path> selected;for(auto h:chosen)selected.push_back(std::move(g.paths[h]));g.paths=std::move(selected);
    std::cerr<<"[M::PanelSelection] input_paths="<<H<<" selected_paths="<<g.paths.size()<<" raw_pairs="<<raw_pairs<<" sketch="<<N<<" k=31 w=25 reference_preserved=yes complete_sample_pairs=yes\n";
}
// Raw-read panel preselection, retaining REF plus 25 complete sample pairs,
// by local diploid read consistency (WindowSelection.hpp).
inline void select_panel(NativeGraph&g,const std::vector<ReadLibrary>&libraries,
                         const std::filesystem::path&out,int threads,size_t cap=51,
                         uint64_t window_bases=25000){
    if(g.paths.size()<=cap)return;
    WindowSelectionOptions options;options.cap=cap;options.window_bases=window_bases;
    std::vector<size_t> chosen;std::vector<long double> gains;
    const auto s=window_select_panel(g,libraries,threads,options,&chosen,&gains);
    const size_t H=g.paths.size();
    std::ofstream report(out/"panel-selection.tsv");report<<"path\tscore\tretained\n";
    for(size_t h=0;h<H;++h)report<<g.paths[h].name<<'\t'<<std::setprecision(17)<<double(gains[h])<<'\t'<<std::binary_search(chosen.begin(),chosen.end(),h)<<'\n';
    report.close();if(!report)throw std::runtime_error("Cannot write panel selection audit");
    std::ofstream summary(out/"window-selection.tsv");summary<<std::setprecision(17);
    summary<<"metric\tvalue\nwindow_bases\t"<<window_bases<<"\nitems\t"<<s.items<<"\nwindows\t"<<s.windows<<"\ncontrols\t"<<s.controls
           <<"\ncontrol_mean\t"<<double(s.mean)<<"\ncontrol_shape\t"<<double(s.shape)<<"\nrepeated\t"<<s.repeated<<"\nwindow_ambiguous\t"<<s.ambiguous
           <<"\nunclassified\t"<<s.unclassified<<"\ninitial_cost\t"<<s.initial_cost<<"\nfinal_cost\t"<<s.final_cost<<'\n';
    for(size_t i=0;i<s.steps.size();++i)summary<<"step"<<i<<'\t'<<s.steps[i].first<<':'<<s.steps[i].second<<'\n';
    summary.close();if(!summary)throw std::runtime_error("Cannot write window selection audit");
    std::vector<NativeGraph::Path> selected;for(auto h:chosen)selected.push_back(std::move(g.paths[h]));g.paths=std::move(selected);
    std::cerr<<"[M::PanelSelection] policy=local_window_read_consistency input_paths="<<H<<" selected_paths="<<g.paths.size()<<" items="<<s.items<<" windows="<<s.windows
             <<" control_mean="<<double(s.mean)<<" cost="<<s.initial_cost<<"->"<<s.final_cost<<" reference_preserved=yes complete_sample_pairs=yes\n";
}
}
