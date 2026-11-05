#pragma once
#include <string>
#include <vector>

namespace dipcall_vcf::output {
// The internal dipcall VCF contract, not a general-purpose VCF normalizer:
// one phased diploid sample, literal DNA alleles or the spanning-deletion
// ALT '*', INFO '.', FORMAT GT:AD.
// Reject unknown fields instead of silently mis-transforming allele annotations.
void normalize(const std::string& reference, const std::string& input,
               const std::string& sorted_vcf);
// Copy an already sorted dipcall VCF (possibly with PS) and build a CSI when
// the destination ends in .gz. The input and every output must be distinct.
void publish(const std::string& input, const std::string& output);
// Join component VCFs without changing GT or PS. Header dictionaries and sample
// identity are validated. All records, including overlapping calls, are kept.
void combine(const std::vector<std::string>& inputs, const std::string& output);
} // namespace dipcall_vcf::output
