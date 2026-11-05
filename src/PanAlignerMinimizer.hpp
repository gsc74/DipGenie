#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

struct MinimizerHit
{
    uint64_t hash = 0;
    uint32_t start = 0;
    uint32_t end = 0;
    uint8_t strand = 0;
};

#include "CanonicalMinimizer.hpp"

namespace dipgenie {

enum class MinimizerMode {
    panaligner_packed,
    canonical_2bit64_unpacked,
    canonical_murmur,
};

inline MinimizerMode minimizer_mode(int32_t k)
{
    if (k < 1 || k > 31) throw std::invalid_argument("supported k range is 1..31");
    return MinimizerMode::canonical_murmur;
}

inline const char* minimizer_mode_name(MinimizerMode mode)
{
    if (mode == MinimizerMode::canonical_murmur) return "canonical-lexicographic-murmur";
    return mode == MinimizerMode::panaligner_packed
        ? "panaligner-packed"
        : "canonical-2bit64-unpacked";
}

// PanAligner/minigraph's reversible integer hash.  Keeping the implementation
// in one header makes read and completed-walk sketches exactly identical.
inline uint64_t panaligner_hash64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask;
    key ^= key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key ^= key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key ^= key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

inline uint8_t nt4(char base)
{
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return 4;
    }
}

template <typename Sequence, typename Emit>
inline void panaligner_for_each_minimizer(
    const Sequence& sequence, int32_t window, int32_t k, Emit&& emit_hit)
{
    // PanAligner/minigraph emits the final minimum even when a sequence has
    // fewer than w k-mers.  Only a sequence shorter than k has no minimizer.
    if (k <= 0 || k > 28 || window <= 0 || window >= 256 ||
        sequence.size() < static_cast<size_t>(k))
        return;

    struct PackedMinimizer {
        uint64_t x = std::numeric_limits<uint64_t>::max();
        uint64_t y = std::numeric_limits<uint64_t>::max();
    };
    const uint64_t mask = (1ULL << (2 * k)) - 1ULL;
    const uint32_t reverse_shift = 2U * static_cast<uint32_t>(k - 1);
    uint64_t kmer[2] = {0, 0};
    int32_t valid = 0, buffer_position = 0, minimum_position = 0;
    int32_t kmer_span = 0;
    std::vector<PackedMinimizer> buffer(static_cast<size_t>(window));
    PackedMinimizer minimum;

    auto emit = [&](const PackedMinimizer& minimizer) {
        const uint32_t end = static_cast<uint32_t>(minimizer.y) >> 1;
        emit_hit(MinimizerHit{
            minimizer.x, end - static_cast<uint32_t>(k) + 1, end,
            static_cast<uint8_t>(minimizer.y & 1ULL)});
    };

    for (uint32_t position = 0; position < sequence.size(); ++position) {
        const uint8_t base = nt4(sequence[position]);
        PackedMinimizer info;
        if (base < 4) {
            kmer_span = std::min(valid + 1, k);
            kmer[0] = ((kmer[0] << 2) | base) & mask;
            kmer[1] = (kmer[1] >> 2) |
                      (static_cast<uint64_t>(3U ^ base) << reverse_shift);
            if (kmer[0] == kmer[1]) continue;
            const uint32_t strand = kmer[0] < kmer[1] ? 0U : 1U;
            ++valid;
            if (valid >= k && kmer_span < 256) {
                info.x = (panaligner_hash64(kmer[strand], mask) << 8) |
                         static_cast<uint64_t>(kmer_span);
                info.y = (static_cast<uint64_t>(position) << 1) | strand;
            }
        } else {
            valid = 0;
            kmer_span = 0;
        }

        buffer[buffer_position] = info;
        if (valid == window + k - 1 &&
            minimum.x != std::numeric_limits<uint64_t>::max()) {
            for (int32_t j = buffer_position + 1; j < window; ++j)
                if (minimum.x == buffer[j].x && buffer[j].y != minimum.y)
                    emit(buffer[j]);
            for (int32_t j = 0; j < buffer_position; ++j)
                if (minimum.x == buffer[j].x && buffer[j].y != minimum.y)
                    emit(buffer[j]);
        }
        if (info.x <= minimum.x) {
            if (valid >= window + k &&
                minimum.x != std::numeric_limits<uint64_t>::max())
                emit(minimum);
            minimum = info;
            minimum_position = buffer_position;
        } else if (buffer_position == minimum_position) {
            if (valid >= window + k - 1 &&
                minimum.x != std::numeric_limits<uint64_t>::max())
                emit(minimum);
            minimum.x = minimum.y = std::numeric_limits<uint64_t>::max();
            for (int32_t j = buffer_position + 1; j < window; ++j)
                if (minimum.x >= buffer[j].x)
                    minimum = buffer[j], minimum_position = j;
            for (int32_t j = 0; j <= buffer_position; ++j)
                if (minimum.x >= buffer[j].x)
                    minimum = buffer[j], minimum_position = j;
            if (valid >= window + k - 1 &&
                minimum.x != std::numeric_limits<uint64_t>::max()) {
                for (int32_t j = buffer_position + 1; j < window; ++j)
                    if (minimum.x == buffer[j].x && minimum.y != buffer[j].y)
                        emit(buffer[j]);
                for (int32_t j = 0; j <= buffer_position; ++j)
                    if (minimum.x == buffer[j].x && minimum.y != buffer[j].y)
                        emit(buffer[j]);
            }
        }
        if (++buffer_position == window) buffer_position = 0;
    }
    if (minimum.x != std::numeric_limits<uint64_t>::max()) emit(minimum);
}

inline std::vector<MinimizerHit> panaligner_minimizers(
    const std::string& sequence, int32_t window, int32_t k)
{
    std::vector<MinimizerHit> result;
    if (window > 0)
        result.reserve(sequence.size() / static_cast<size_t>(window) + 2);
    panaligner_for_each_minimizer(
        sequence, window, k,
        [&](const MinimizerHit& hit) { result.push_back(hit); });
    return result;
}

// PanAligner's packed key reserves eight low bits for the k-mer span, which
// limits that representation to k<=28.  Ordinary (non-HPC) DNA always has
// span k, so k=29..31 can use the same canonical rolling 2-bit code and the
// same reversible hash without packing the redundant span.  The resulting
// key occupies at most 62 bits at k=31.  Coordinates and strand are carried
// explicitly in MinimizerHit, making this representation lossless here.
inline uint64_t canonical_2bit64_key(uint64_t forward, uint64_t reverse,
                                    int32_t k)
{
    if (k <= 0 || k > 31) return std::numeric_limits<uint64_t>::max();
    const uint64_t mask = (1ULL << (2 * k)) - 1ULL;
    return panaligner_hash64(std::min(forward, reverse), mask);
}

template <typename Sequence, typename Emit>
inline void canonical_2bit64_for_each_minimizer(
    const Sequence& sequence, int32_t window, int32_t k, Emit&& emit_hit)
{
    if (k <= 0 || k > 31 || window <= 0 || window >= 256 ||
        sequence.size() < static_cast<size_t>(k))
        return;

    struct Minimizer {
        uint64_t key = std::numeric_limits<uint64_t>::max();
        uint32_t end = std::numeric_limits<uint32_t>::max();
        uint8_t strand = 0;
    };

    const uint64_t mask = (1ULL << (2 * k)) - 1ULL;
    const uint32_t reverse_shift = 2U * static_cast<uint32_t>(k - 1);
    uint64_t forward = 0, reverse = 0;
    int32_t valid_bases = 0;
    int32_t buffer_position = 0, minimum_position = 0;
    std::vector<Minimizer> buffer(static_cast<size_t>(window));
    Minimizer minimum;

    const auto is_valid = [](const Minimizer& minimizer) {
        return minimizer.key != std::numeric_limits<uint64_t>::max();
    };
    const auto emit = [&](const Minimizer& minimizer) {
        emit_hit(MinimizerHit{
            minimizer.key, minimizer.end - static_cast<uint32_t>(k) + 1,
            minimizer.end, minimizer.strand});
    };
    const auto reset_window = [&]() {
        std::fill(buffer.begin(), buffer.end(), Minimizer{});
        minimum = Minimizer{};
        buffer_position = 0;
        minimum_position = 0;
    };
    const auto finish_run = [&]() {
        if (is_valid(minimum)) emit(minimum);
    };

    for (uint32_t position = 0; position < sequence.size(); ++position) {
        const uint8_t base = nt4(sequence[position]);
        if (base >= 4) {
            if (valid_bases > 0) {
                finish_run();
                reset_window();
            }
            forward = reverse = 0;
            valid_bases = 0;
            continue;
        }

        forward = ((forward << 2) | base) & mask;
        reverse = (reverse >> 2) |
                  (static_cast<uint64_t>(3U ^ base) << reverse_shift);
        ++valid_bases;

        Minimizer info;
        if (valid_bases >= k) {
            info.key = canonical_2bit64_key(forward, reverse, k);
            info.end = position;
            // Canonical ties are deterministic. They cannot occur for odd k,
            // including the k=31 ablation, but choosing the forward strand
            // also makes the helper well-defined for even k.
            info.strand = forward <= reverse ? 0U : 1U;
        }

        buffer[static_cast<size_t>(buffer_position)] = info;
        if (valid_bases == window + k - 1 && is_valid(minimum)) {
            // Emit every tied minimum in the first complete window. This is
            // the symmetric-minimizer tie rule used by PanAligner/minigraph.
            for (int32_t j = buffer_position + 1; j < window; ++j)
                if (minimum.key == buffer[static_cast<size_t>(j)].key &&
                    minimum.end != buffer[static_cast<size_t>(j)].end)
                    emit(buffer[static_cast<size_t>(j)]);
            for (int32_t j = 0; j < buffer_position; ++j)
                if (minimum.key == buffer[static_cast<size_t>(j)].key &&
                    minimum.end != buffer[static_cast<size_t>(j)].end)
                    emit(buffer[static_cast<size_t>(j)]);
            // When the entering k-mer ties the tracked (rightmost) minimum,
            // the update below replaces that tracked occurrence. Emit it now
            // so every tied minimum in the first window appears exactly once.
            if (info.key == minimum.key) emit(minimum);
        }
        if (info.key <= minimum.key) {
            if (valid_bases >= window + k && is_valid(minimum)) emit(minimum);
            minimum = info;
            minimum_position = buffer_position;
        } else if (buffer_position == minimum_position) {
            if (valid_bases >= window + k - 1 && is_valid(minimum))
                emit(minimum);
            minimum = Minimizer{};
            for (int32_t j = buffer_position + 1; j < window; ++j) {
                if (minimum.key >= buffer[static_cast<size_t>(j)].key) {
                    minimum = buffer[static_cast<size_t>(j)];
                    minimum_position = j;
                }
            }
            for (int32_t j = 0; j <= buffer_position; ++j) {
                if (minimum.key >= buffer[static_cast<size_t>(j)].key) {
                    minimum = buffer[static_cast<size_t>(j)];
                    minimum_position = j;
                }
            }
            if (valid_bases >= window + k - 1 && is_valid(minimum)) {
                for (int32_t j = buffer_position + 1; j < window; ++j)
                    if (minimum.key == buffer[static_cast<size_t>(j)].key &&
                        minimum.end != buffer[static_cast<size_t>(j)].end)
                        emit(buffer[static_cast<size_t>(j)]);
                for (int32_t j = 0; j <= buffer_position; ++j)
                    if (minimum.key == buffer[static_cast<size_t>(j)].key &&
                        minimum.end != buffer[static_cast<size_t>(j)].end)
                        emit(buffer[static_cast<size_t>(j)]);
            }
        }
        if (++buffer_position == window) buffer_position = 0;
    }
    finish_run();
}

inline std::vector<MinimizerHit> canonical_2bit64_minimizers(
    const std::string& sequence, int32_t window, int32_t k)
{
    std::vector<MinimizerHit> result;
    if (window > 0)
        result.reserve(sequence.size() / static_cast<size_t>(window) + 2);
    canonical_2bit64_for_each_minimizer(
        sequence, window, k,
        [&](const MinimizerHit& hit) { result.push_back(hit); });
    return result;
}

template <typename Sequence, typename Emit>
inline void for_each_sequence_minimizer(
    const Sequence& sequence, int32_t window, int32_t k, Emit&& emit_hit)
{
    canonical_minimizer::sketch(sequence, window, k, std::forward<Emit>(emit_hit));
}

inline std::vector<uint64_t> sequence_minimizer_hashes(
    const std::string& sequence, int32_t window, int32_t k)
{
    std::vector<uint64_t> result;
    if (window > 0)
        result.reserve(sequence.size() / static_cast<size_t>(window) + 2);
    for_each_sequence_minimizer(
        sequence, window, k,
        [&](const MinimizerHit& hit) { result.push_back(hit.hash); });
    return result;
}

inline std::vector<MinimizerHit> sequence_minimizers(
    const std::string& sequence, int32_t window, int32_t k)
{
    std::vector<MinimizerHit> result;
    for_each_sequence_minimizer(sequence, window, k,
        [&](const MinimizerHit& hit) { result.push_back(hit); });
    return result;
}

}  // namespace dipgenie
