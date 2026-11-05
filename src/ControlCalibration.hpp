#pragma once
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace hic_calibration {
// Moment calibration on single-locus strings present in every panel walk.
// Zeros are included. This estimates a marginal count distribution, not phase.
struct Controls {
    uint64_t n = 0;
    long double sum = 0, squared = 0;
    void add(uint64_t count, uint64_t frequency = 1) {
        if (frequency > UINT64_MAX - n) throw std::overflow_error("Control count");
        n += frequency;
        const auto x = static_cast<long double>(count);
        sum += x * frequency;
        squared += x * x * frequency;
        if (!std::isfinite(sum) || !std::isfinite(squared))
            throw std::overflow_error("Control moments");
    }
};
struct Calibration {
    long double mean, variance, shape;
    static Calibration fit(const Controls& c) {
        if (!c.n || !(c.sum > 0))
            throw std::invalid_argument("No observed conserved single-locus calibration controls");
        const auto mean = c.sum / c.n;
        const auto variance = std::fmax(0.L, c.squared / c.n - mean * mean);
        const auto shape = variance > mean ? mean * mean / (variance - mean) :
                           std::numeric_limits<long double>::infinity();
        return {mean, variance, shape};
    }
    // log NB(c; mean, shape) - log NB(c; mean/2, shape).
    // Equal priors; common factorial/Gamma terms cancel exactly.
    long double hom_log_ratio(uint64_t count) const {
        if (!(mean > 0) || !(shape > 0) || !std::isfinite(mean))
            throw std::invalid_argument("Invalid negative-binomial calibration");
        const auto c = static_cast<long double>(count);
        if (std::isinf(shape)) return c * std::log(2.L) - mean / 2;
        return c * std::log(2.L) - (c + shape) * std::log1p(mean / (2 * shape + mean));
    }
    bool hom(uint64_t count) const {
        if (!count) throw std::invalid_argument("Classify observed targets only");
        // Retain v1's singleton convention, uniformly at every depth.
        return count > 1 && hom_log_ratio(count) > 0;
    }
};
}
