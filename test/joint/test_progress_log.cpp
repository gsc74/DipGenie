#include "ProgressLog.hpp"
#include "PresenceDP.hpp"
#include <cassert>
#include <iostream>

int main() {
    std::ostringstream screen,audit;double seconds=0;
    hic::ProgressLog log(screen,audit,[&]{return seconds;});
    log.stage("start");log.detail("large diagnostic arrays");
    log.progress("DP forward",1,100);
    for(int i=2;i<60;++i){seconds=i;log.progress("DP forward",i,100);}
    seconds=60;log.progress("DP forward",60,100);
    seconds=61;log.progress("DP forward",100,100);
    log.stage("FASTA ready: rho=1.123456 coverage=1.000000 phase=1.234567 (objective bounds)");
    const auto s=screen.str();
    assert(s.find("large diagnostic")==std::string::npos);
    assert(audit.str().find("large diagnostic")!=std::string::npos);
    assert(s.find("60.0%")!=std::string::npos&&s.find("59.0%") ==std::string::npos);
    assert(s.find("100.0%")!=std::string::npos&&s.find("rho=1.123456")!=std::string::npos);
    assert(s.find("[M::DipGenie::60.000]")!=std::string::npos);
    bool caught=false;try{log.progress("bad",2,1);}catch(const std::invalid_argument&){caught=true;}
    assert(caught);
    // Unknown-size FASTQ streams report counts, never an invented percentage
    // or a pre-counted total; each library has an independent throttle.
    std::ostringstream reads,read_audit;seconds=0;
    hic::ProgressLog streaming(reads,read_audit,[&]{return seconds;});
    streaming.pairs("reads SR/library1",0,0);
    seconds=1;streaming.pairs("reads SR/library1",4096,3000);
    seconds=60;streaming.pairs("reads SR/library1",8192,6000);
    seconds=61;streaming.pairs("reads SR/library1",8200,6001,true);
    streaming.pairs("reads HIC/library2",0,0);
    streaming.pairs("reads HIC/library2",100,0,true);
    const auto r=reads.str();
    assert(r==read_audit.str()&&r.find('%')==std::string::npos);
    assert(r.find("pairs=4096")==std::string::npos);
    assert(r.find("pairs=8192 linked=6000")!=std::string::npos);
    assert(r.find("pairs=8200 linked=6001 (complete)")!=std::string::npos);
    assert(r.find("reads HIC/library2: pairs=0 linked=0")!=std::string::npos);
    assert(r.find("reads HIC/library2: pairs=100 linked=0 (complete)")!=std::string::npos);
    caught=false;try{streaming.pairs("bad",1,2);}catch(const std::invalid_argument&){caught=true;}
    assert(caught);
    // Compacted empty segments must not leave the displayed pass at B-1/B.
    // Cover empty prefix/interior/suffix, a single occupied segment, and no
    // evidence, with and without pruning and at every supported test width.
    for(const auto&occupied:std::vector<std::vector<int>>{{0,2,4},{2},{0},{4},{}})
        for(int threads:{1,2,4,8})for(size_t beam:{size_t(0),size_t(2)}) {
            hic::PresenceModel model;model.H=2;model.blocks=5;model.budget=2;model.weights={1,1};
            for(int b:occupied)model.targets.push_back({false,{{b,1}},1,0});
            const auto expected=hic::solve_presence(model,beam,threads,{},true);
            std::ostringstream console,details;double now=0;
            hic::ProgressLog dp_log(console,details,[&]{return now;});
            std::vector<int> completed;
            const auto actual=hic::solve_presence(model,beam,threads,[&](int b,const hic::Result&){
                assert(b>0&&b<=model.blocks);
                assert(completed.empty()||b>completed.back());
                completed.push_back(b);dp_log.progress("DP forward",b,model.blocks);
            },true);
            assert(!completed.empty()&&completed.back()==model.blocks);
            assert(console.str().find("100.0% (5/5)")!=std::string::npos);
            assert(actual.path==expected.path&&actual.score==expected.score&&actual.upper==expected.upper);
            assert(actual.generated==expected.generated&&actual.discarded==expected.discarded);
            assert(actual.peak==expected.peak&&actual.max_frontier==expected.max_frontier);
        }
    std::cout<<"PASS concise progress: throttled updates, full diagnostic log, final certificate retained\n";
}
