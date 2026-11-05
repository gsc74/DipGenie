#pragma once
#include "ReadLibraries.hpp"
#include <memory>

namespace hic {
struct ComponentReference {std::string name,sequence;};
struct PairedBinSummary {
    uint64_t pairs=0,bases=0;
};
struct PairedBins {
    std::vector<ReadLibrary> libraries;
    std::vector<PairedBinSummary> summaries;
    uint64_t input_pairs=0,cross_component_pairs=0,unassigned_pairs=0;
    // Cross-component and unassigned pairs are preserved separately. They are
    // never relabelled as within-component links or duplicated among bins.
    std::filesystem::path cross_component,unassigned;
};
// Assigns each read pair to the graph component on which both ends have a
// unique native placement (NativeLocalizer over the stored component references).
class PairedReadBinner {
    class Impl;
    std::unique_ptr<Impl> impl_;
public:
    PairedReadBinner(const std::vector<ComponentReference>&references,int threads);
    ~PairedReadBinner();
    PairedReadBinner(const PairedReadBinner&)=delete;
    PairedReadBinner&operator=(const PairedReadBinner&)=delete;
    PairedBins partition(const ReadLibrary&library,const std::filesystem::path&fresh_directory);
};
}
