#pragma once
#include <cstdint>
#include <string>
#include <array>
#include <functional>
#include <utility>
#include <vector>

namespace dipcall_vcf {
struct Config {
    std::string reference_fasta;
    std::string output_vcf;
    std::string sample = "DipGenie";
    uint32_t threads = 1;
    std::vector<uint32_t> cpu_ids;
    // 0-based half-open reference intervals, each its own phase set. When
    // nonempty, every heterozygous call gets PS: 1 outside all intervals (the
    // linked chromosome-wide block), interval start+2 inside an interval.
    // Empty: v1 output (phased GT, no PS).
    std::vector<std::pair<uint64_t,uint64_t>> separate;
    // Read-linked phase blocks: sorted, disjoint 0-based half-open intervals,
    // each with its phase-set id; a heterozygous call gets the id of the
    // interval containing it (its own phase set if none). Takes precedence
    // over `separate` when nonempty.
    struct Block {uint64_t begin=0,end=0,ps=0;};
    std::vector<Block> blocks;
    // Optional per-variant phase support. For each heterozygous call in the
    // linked block (PS 1) it receives both local haplotype sequences (reference
    // with that haplotype's called alleles applied, +-flank bases) and the
    // [begin,end) span of the variant's allele in each; a call with
    // insufficient support gets its own phase set (not claimed).
    std::function<bool(const std::array<std::string,2>&,const std::array<std::pair<size_t,size_t>,2>&)> site_supported;
};
struct SiteRecord {uint64_t pos;std::string ref;std::array<std::string,2> allele;bool linked;};
// Local haplotype windows for record i (exposed for tests).
void local_haplotypes(const std::string& reference,const std::vector<SiteRecord>& records,size_t i,size_t flank,
                      std::array<std::string,2>& haps,std::array<std::pair<size_t,size_t>,2>& span);
// Rewrite linked calls failing the predicate to a singleton phase set. Returns count.
// Linked: PS 1, or (any_block) any phase set below the singleton range.
size_t assign_site_phase_sets(const std::string& source,const std::string& dest,const std::string& reference,
    const std::function<bool(const std::array<std::string,2>&,const std::array<std::pair<size_t,size_t>,2>&)>& supported,
    bool any_block=false);
// Assign PS from read-linked blocks (Config::blocks). Returns the number of
// heterozygous calls that fall in no block (own phase set).
size_t assign_block_phase_sets(const std::string& source,const std::string& dest,const std::vector<Config::Block>& blocks);
// Rewrite a VCF assigning PS as described above. Returns the number of
// heterozygous calls placed in separate (non-chromosome-wide) phase sets.
size_t assign_phase_sets(const std::string& source,const std::string& dest,
                         const std::vector<std::pair<uint64_t,uint64_t>>& intervals);

// Assembly alignment, BAM preparation, calling, normalization, callable BEDs,
// compression and CSI indexing are linked; no external executable is invoked.
// Truth and benchmark BEDs are never used during variant calling.
std::string run(const Config&, const std::vector<std::string>& haplotypes);
} // namespace dipcall_vcf
