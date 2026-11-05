#pragma once

#include "gfa.h"

#include <cstdint>
#include <functional>
#include <string>
#include <unordered_set>
#include <vector>

// Standalone phased-panel VCF input adapted from PanGenie's VariantReader and
// GraphBuilder implementations.
// PanGenie is Copyright (c) 2020 Jana Ebler and distributed under the MIT
// License; see third_party/PanGenie-LICENSE.md.
namespace pangenie_vcf {

struct PanelContig {
    std::string chromosome;
    uint64_t reference_length = 0;
    uint64_t input_variants = 0;
    uint64_t retained_variants = 0;
    uint64_t skipped_variants = 0;
};

struct ReferenceSupportCensus {
    std::vector<uint64_t> read_shared_hashes;
    std::vector<uint64_t> component_unique_hashes;
};

struct ComponentGate {
    const std::vector<uint64_t>* read_spectrum_hashes = nullptr;
    // Disk-backed WGS callers provide this exact query instead of
    // materializing the complete read-unique spectrum in RAM.
    std::function<ReferenceSupportCensus(
        const std::vector<std::vector<uint64_t>>&)> reference_support_query;
    int32_t k = 0;
    int32_t w = 0;
    int32_t num_threads = 1;
    uint64_t minimum_anchor_hashes = 0;
    uint32_t minimum_dominance_basis_points = 0;
    uint32_t minimum_reference_basis_points = 0;
};

struct ComponentSelection {
    uint32_t input_component_id = 0;
    std::string chromosome;
    uint64_t reference_length = 0;
    uint64_t reference_shared_anchors = 0;
    uint64_t reference_unique_anchors = 0;
    uint64_t maximum_unique_anchors = 0;
    uint64_t absolute_threshold = 0;
    uint64_t dominance_threshold = 0;
    uint64_t reference_threshold = 0;
    uint64_t effective_threshold = 0;
    bool admitted = false;
};

struct TabixStatus {
    bool usable = false;
    std::string reason;
};

struct PanelGraph {
    gfa_t* graph = nullptr;
    // Compatibility view for the established single-contig interface.  These
    // fields are populated only when contigs.size() == 1; multi-contig callers
    // must use contigs so a chromosome length can never be applied globally.
    std::string chromosome;
    uint64_t reference_length = 0;
    std::vector<PanelContig> contigs;
    std::vector<std::string> samples;
    uint64_t input_samples = 0;
    uint64_t excluded_samples = 0;
    uint64_t input_variants = 0;
    uint64_t retained_variants = 0;
    uint64_t skipped_variants = 0;
    // In indexed mode, contigs contains only materialized components while
    // component_selection retains the ordered global panel audit.
    std::vector<ComponentSelection> component_selection;
    std::string selection_mode = "legacy-full-import";
};

TabixStatus inspect_tabix(const std::string& vcf_filename);
bool component_gate_enabled(uint64_t minimum_anchor_hashes,
                            uint32_t minimum_dominance_basis_points,
                            uint32_t minimum_reference_basis_points);

PanelGraph read_panel(const std::string& vcf_filename,
                      const std::string& reference_filename,
                      const std::unordered_set<std::string>& excluded_samples = {},
                      const ComponentGate* component_gate = nullptr);

}  // namespace pangenie_vcf
