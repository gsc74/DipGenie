#pragma once
#include <array>
#include <cstdint>
#include <string>
#include <vector>

// In-process pieces of the pinned dipcall 0.3 workflow. This is assembly
// alignment, not read placement or a new inference objective. The production
// caller is switched only after end-to-end VCF/BED parity is established.
namespace dipcall_vcf::embedded {
struct AlignmentStats {
    uint64_t queries=0, alignments=0, retained_alignments=0;
};
// Produces <prefix>.hap{1,2}.paf and coordinate-sorted .bam. PAF retains all
// alignments; BAM applies dipcall samflt's MAPQ >= 5, MID block length >= 50kb,
// and nonsecondary filters. Supplementary alignments are not discarded.
std::array<AlignmentStats,2> align(const std::string& reference,
    const std::array<std::string,2>& fasta, const std::string& prefix, unsigned threads,
    const std::vector<uint32_t>& cpu_ids={});
// HTSbox pileup -q5 -evcf, linked as a library; no process, shell, or stdout
// redirection. Input BAMs must be sorted and share the reference dictionary.
void pileup(const std::string& reference,const std::array<std::string,2>& bam,
            const std::string& output);
// dipcall-aux.js vcfpair -a -s SAMPLE for the diploid, all-contigs policy used
// by DipGenie. Raw GAP/HET calls remain here; final PASS filtering is separate.
void pair_vcf(const std::string& input,const std::string& output,const std::string& sample);
// paftools.js call's R intervals and dipcall's bedtk isec -m. Uses all primary
// MAPQ>=5 alignments >=10kb to exclude overlap; the sole covering alignment
// must be >=50kb. No truth BED participates in calling.
void callable_bed(const std::array<std::string,2>& paf,const std::string& prefix);
}
