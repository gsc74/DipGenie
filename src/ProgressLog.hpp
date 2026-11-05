#pragma once
#include <chrono>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <map>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace hic {
// Diagnostics are complete in inference.log; stderr carries only milestones,
// at most one progress update per key/minute, and the final certificates.
// Called by the coordinator (not the worker threads).
class ProgressLog {
    std::ostream&screen_;
    std::ostream&audit_;
    std::function<double()> clock_;
    std::map<std::string,double> last_;
    std::string line(const std::string&message)const {
        std::ostringstream out;out<<"[M::DipGenie::"<<std::fixed<<std::setprecision(3)<<clock_()<<"] "<<message<<'\n';
        return out.str();
    }
    static void emit(std::ostream&out,const std::string&text) {
        out<<text;out.flush();if(!out)throw std::runtime_error("Progress log write failed");
    }
public:
    ProgressLog(std::ostream&screen,std::ostream&audit,std::function<double()> clock={})
        :screen_(screen),audit_(audit),clock_(std::move(clock)) {
        if(!clock_) {
            const auto begin=std::chrono::steady_clock::now();
            clock_=[begin]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();};
        }
    }
    void detail(const std::string&message){emit(audit_,line(message));}
    void stage(const std::string&message) {
        const auto text=line(message);emit(audit_,text);emit(screen_,text);
    }
    void progress(const std::string&key,uint64_t done,uint64_t total) {
        if(!total||done>total)throw std::invalid_argument("Progress count");
        const double now=clock_();auto it=last_.find(key);
        if(it!=last_.end()&&done!=total&&now-it->second<60)return;
        last_[key]=now;
        std::ostringstream out;out<<key<<' '<<std::fixed<<std::setprecision(1)
            <<100.L*done/total<<"% ("<<done<<'/'<<total<<')';
        stage(out.str());
    }
    // FASTQs are streamed once: their total is not known in advance. Linked
    // counts retained paired-marker observations, not calibrated phase factors.
    void pairs(const std::string&key,uint64_t processed,uint64_t linked,bool complete=false) {
        if(linked>processed)throw std::invalid_argument("Linked pair count exceeds processed pairs");
        const double now=clock_();auto it=last_.find(key);
        if(it!=last_.end()&&!complete&&now-it->second<60)return;
        last_[key]=now;
        std::ostringstream out;out<<key<<": pairs="<<processed<<" linked="<<linked;
        if(complete)out<<" (complete)";
        stage(out.str());
    }
};
}
