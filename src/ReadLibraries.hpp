#pragma once
#include "BuildModel.hpp"
#include <filesystem>

namespace hic {
struct ReadLibrary {std::string name,kind,r1,r2;};
inline void validate_libraries(std::vector<ReadLibrary>&libraries) {
    if(libraries.empty())throw std::invalid_argument("No read libraries");
    std::set<std::string> names,reads;
    auto path=[](const std::string&s) {
        if(!std::filesystem::is_regular_file(s)||!std::filesystem::file_size(s))
            throw std::invalid_argument("Missing/empty library input: "+s);
        return std::filesystem::canonical(s).string();
    };
    for(auto&l:libraries) {
        if(l.name.empty()||!names.insert(l.name).second||(l.kind!="SR"&&l.kind!="HIC"))
            throw std::invalid_argument("Invalid/duplicate read library");
        l.r1=path(l.r1);l.r2=path(l.r2);
        if(!reads.insert(l.r1).second||!reads.insert(l.r2).second)
            throw std::invalid_argument("A FASTQ is supplied more than once");
    }
}
struct LibraryEvidence {
    RawEvidence combined;
    std::vector<CountLibrary> counts;
    // contact_pairs counts candidates before linkage admission; marker_pairs
    // counts physical pairs with usable markers on both ends after admission.
    std::vector<uint64_t> pairs,contact_pairs,marker_pairs;
};
template<class Read> inline LibraryEvidence read_library_evidence_with(SequenceIndex&index,
        const std::vector<ReadLibrary>&libraries,Read read) {
    if(libraries.size()>=UINT32_MAX)throw std::overflow_error("Too many read libraries");
    LibraryEvidence result;std::vector<uint64_t> totals(index.targets.size());
    for(const auto&l:libraries) {
        for(auto&t:index.targets){t.count=0;t.occurrence_count=0;}
        auto raw=read(l,uint32_t(result.counts.size()));
        if(!raw.pairs)throw std::invalid_argument("Read library has no pairs: "+l.name);
        CountLibrary counts{l.name,l.kind,{},{}};counts.counts.reserve(index.targets.size());
        counts.occurrences.reserve(index.targets.size());
        for(size_t x=0;x<index.targets.size();++x) {
            const auto c=index.targets[x].count;
            if(c>UINT64_MAX-totals[x])throw std::overflow_error("Read-library count overflow");
            totals[x]+=c;counts.counts.push_back(c);counts.occurrences.push_back(index.targets[x].occurrence_count);
        }
        result.counts.push_back(std::move(counts));
        result.pairs.push_back(raw.pairs);result.contact_pairs.push_back(raw.contact_pairs);
        result.marker_pairs.push_back(raw.markers.size());
        result.combined.pairs+=raw.pairs;result.combined.contact_pairs+=raw.contact_pairs;
        if(raw.target_end_observations>UINT64_MAX-result.combined.target_end_observations||
           raw.minimizer_target_end_observations>UINT64_MAX-result.combined.minimizer_target_end_observations)
            throw std::overflow_error("Library target observation audit");
        result.combined.target_end_observations+=raw.target_end_observations;
        result.combined.minimizer_target_end_observations+=raw.minimizer_target_end_observations;
        for(auto&marker:raw.markers) {
            marker.library=uint32_t(result.counts.size()-1);
            result.combined.markers.push_back(std::move(marker));
        }
        // Histogram multiplicities are per-library types, not pooled word counts.
        for(auto[count,n]:raw.count_histogram)result.combined.count_histogram[count]+=n;
    }
    for(size_t x=0;x<index.targets.size();++x)index.targets[x].count=totals[x];
    return result;
}
}
