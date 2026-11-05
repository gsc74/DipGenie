#pragma once
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include "PackedPresence.hpp"

namespace hic {
// Project a sorted frontier onto its surviving entries. Consecutive old
// positions are copied in bulk; new targets start with zero presence bits.
// The plan is shared by all live states at one segment, not rebuilt per state.
class FrontierCopy {
    struct Run {size_t destination,source,length;};
    std::vector<Run> runs_;
    size_t length_=0,source_end_=0;
public:
    explicit FrontierCopy(const std::vector<int>& positions):length_(positions.size()) {
        for(size_t i=0;i<positions.size();++i) {
            if(positions[i]<-1)throw std::invalid_argument("Invalid frontier position");
            if(positions[i]<0)continue;
            const size_t source=size_t(positions[i]);
            if(!runs_.empty()&&runs_.back().destination+runs_.back().length==i&&
               runs_.back().source+runs_.back().length==source)++runs_.back().length;
            else runs_.push_back({i,source,1});
            source_end_=std::max(source_end_,source+1);
        }
    }
    std::vector<unsigned char> gather(const std::vector<unsigned char>& source)const {
        if(source.size()<source_end_)throw std::logic_error("Frontier source is too short");
        std::vector<unsigned char> result(length_,0);
        for(const auto&r:runs_)
            std::copy_n(source.data()+r.source,r.length,result.data()+r.destination);
        return result;
    }
    PackedPresence gather(const PackedPresence&source)const {
        if(source.size()<source_end_)throw std::logic_error("Frontier source is too short");
        PackedPresence result(length_);
        for(const auto&r:runs_)result.copy_zero_range(source,r.source,r.destination,r.length);
        return result;
    }
};
}
