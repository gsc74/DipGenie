#pragma once

#include <algorithm>
#include <cstdint>
#include <stdexcept>

namespace dipgenie::component_admission {

// Reference breadth is the union of primary-chain intervals, not their sum:
// many reads mapping to one repeat must not admit an entire chromosome.
constexpr uint32_t DEFAULT_CHAIN_BREADTH_BPS = 500;

struct PreflightThreshold {
    uint64_t absolute = 0;
    uint64_t dominance = 0;
    uint64_t reference_span = 0;
    uint64_t effective = 0;
};

// Compute ceil(value * basis_points / 10,000) without a potentially
// overflowing multiplication.  One percent is exactly 100 basis points.
inline uint64_t ceil_basis_points(uint64_t value, uint32_t basis_points)
{
    constexpr uint64_t denominator = 10000;
    if (basis_points > denominator)
        throw std::invalid_argument(
            "component dominance basis points exceed 10000");
    const uint64_t quotient = value / denominator;
    const uint64_t remainder = value % denominator;
    const uint64_t whole = quotient * basis_points;
    const uint64_t fractional = remainder == 0 || basis_points == 0
        ? 0
        : (remainder * basis_points - 1) / denominator + 1;
    return whole + fractional;
}

inline PreflightThreshold make_preflight_threshold(
    uint64_t maximum_component_unique_hashes, uint64_t absolute_threshold,
    uint32_t dominance_basis_points, uint64_t component_reference_span = 0,
    uint32_t reference_span_basis_points = 0)
{
    // The original dominance gate remains coupled to its absolute floor.
    // The component-local span gate is independent and can therefore be the
    // sole admission rule.
    if (absolute_threshold == 0 && reference_span_basis_points == 0)
        return {};
    const uint64_t dominance = absolute_threshold == 0
        ? 0
        : ceil_basis_points(maximum_component_unique_hashes,
                            dominance_basis_points);
    const uint64_t reference_span = ceil_basis_points(
        component_reference_span, reference_span_basis_points);
    return {absolute_threshold, dominance, reference_span,
            std::max({absolute_threshold, dominance, reference_span})};
}

inline bool passes_preflight(uint64_t component_unique_hashes,
                             const PreflightThreshold& threshold)
{
    return threshold.effective == 0 ||
           component_unique_hashes >= threshold.effective;
}

inline bool passes_chain_breadth(uint64_t covered_bases,
                                uint64_t reference_length,
                                uint32_t basis_points)
{
    if (covered_bases > reference_length)
        throw std::invalid_argument("chain breadth exceeds reference length");
    const uint64_t threshold = ceil_basis_points(reference_length, basis_points);
    if (basis_points == 0) return true;
    if (reference_length == 0)
        throw std::invalid_argument("chain breadth requires a nonempty reference");
    return covered_bases >= threshold;
}

// A zero threshold disables postflight filtering.  Otherwise the component
// must retain at least the configured number of classified anchor hashes.
inline bool passes_postflight(uint64_t retained_anchor_hashes,
                              uint64_t threshold)
{
    return threshold == 0 || retained_anchor_hashes >= threshold;
}

}  // namespace dipgenie::component_admission
