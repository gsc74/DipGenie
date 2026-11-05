#pragma once
#include "MurmurHash.hpp"
#include <cctype>
#include <deque>
#include <stdexcept>
#include <type_traits>

namespace canonical_minimizer {
// Selection preserves v1's lexicographic canonical minimum, rightmost tie,
// full-window requirement, Murmur seed 0, and consecutive-HASH suppression.
// Two-bit comparisons accelerate only A/C/G/T strings; all other symbols use
// the literal v1 rule (unknown complements are unchanged).
inline uint64_t hash(const char *s, int k) {
  uint64_t result[2];
  dipgenie_murmur::MurmurHash3_x64_128(s, k, 0, result);
  return result[0] ^ result[1];
}
template <class Sequence, class Emit>
void sketch(const Sequence &input, int w, int k, Emit &&emit) {
  if constexpr (!std::is_same_v<std::decay_t<Sequence>, std::string>) {
    // A graph spelling view can be forward-only. Materialize it once;
    // selection's DNA validation and canonical scan require two passes.
    std::string sequence;
    sequence.reserve(input.size());
    for (size_t i = 0; i < input.size(); ++i)
      sequence.push_back(input[i]);
    sketch(sequence, w, k, std::forward<Emit>(emit));
    return;
  }
  if (k < 1 || k > 31 || w < 1)
    throw std::invalid_argument("canonical sketch requires 1<=k<=31 and w>=1");
  if (input.size() > UINT32_MAX)
    throw std::overflow_error("v1 sketch coordinate overflow");
  if (input.size() < static_cast<size_t>(k) + w - 1)
    return;
  bool dna = true;
  for (size_t i = 0; i < input.size(); ++i) {
    const unsigned char c = input[i];
    if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
      dna = false;
      break;
    }
  }
  uint64_t previous_hash = UINT64_MAX;
  if (dna) {
    struct Item {
      uint64_t key;
      uint32_t position;
      uint8_t strand;
    };
    std::deque<Item> queue;
    const uint64_t mask = (UINT64_C(1) << (2 * k)) - 1;
    uint64_t forward = 0, reverse = 0, last_key = UINT64_MAX, last_hash = 0;
    char text[31];
    for (uint32_t end = 0; end < input.size(); ++end) {
      const char c = input[end];
      const unsigned base = c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : 3;
      forward = ((forward << 2) | base) & mask;
      reverse =
          (reverse >> 2) | (static_cast<uint64_t>(3 ^ base) << (2 * (k - 1)));
      if (end + 1U < static_cast<unsigned>(k))
        continue;
      const uint32_t start = end + 1 - k;
      const uint64_t key = std::min(forward, reverse);
      while (!queue.empty() && queue.back().key >= key)
        queue.pop_back();
      queue.push_back({key, start, static_cast<uint8_t>(reverse < forward)});
      if (start >= static_cast<unsigned>(w) &&
          queue.front().position <= start - w)
        queue.pop_front();
      if (start + 1U < static_cast<unsigned>(w))
        continue;
      const auto &best = queue.front();
      if (best.key != last_key) {
        auto code = best.key;
        for (int j = k - 1; j >= 0; --j) {
          text[j] = "ACGT"[code & 3];
          code >>= 2;
        }
        last_hash = hash(text, k);
        last_key = best.key;
      }
      if (last_hash != previous_hash) {
        emit(MinimizerHit{last_hash, best.position, best.position + k - 1U,
                          best.strand});
        previous_hash = last_hash;
      }
    }
    return;
  }
  std::string s;
  s.reserve(input.size());
  for (size_t i = 0; i < input.size(); ++i)
    s.push_back(input[i]);
  for (char &c : s)
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  struct Item {
    std::string key;
    uint32_t position;
    uint8_t strand;
  };
  std::deque<Item> queue;
  for (uint32_t start = 0; start <= s.size() - static_cast<size_t>(k);
       ++start) {
    std::string f = s.substr(start, k), r(f.rbegin(), f.rend());
    for (char &c : r) {
      if (c == 'A')
        c = 'T';
      else if (c == 'T')
        c = 'A';
      else if (c == 'C')
        c = 'G';
      else if (c == 'G')
        c = 'C';
    }
    const uint8_t strand = r < f;
    std::string key = std::min(f, r);
    while (!queue.empty() && queue.back().key >= key)
      queue.pop_back();
    queue.push_back({std::move(key), start, strand});
    if (start >= static_cast<unsigned>(w) &&
        queue.front().position <= start - w)
      queue.pop_front();
    if (start + 1U < static_cast<unsigned>(w))
      continue;
    const auto &best = queue.front();
    const uint64_t value = hash(best.key.data(), k);
    if (value != previous_hash) {
      emit(MinimizerHit{value, best.position, best.position + k - 1U,
                        best.strand});
      previous_hash = value;
    }
  }
}
} // namespace canonical_minimizer
