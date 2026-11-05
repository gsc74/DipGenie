// The VCF validation, phased-path parsing, and graph construction in this file
// are adapted from PanGenie's VariantReader and GraphBuilder implementations.
//
// PanGenie is Copyright (c) 2020 Jana Ebler and distributed under the MIT
// License; see third_party/PanGenie-LICENSE.md.

#include "PanGenieVcf.hpp"

#include "ComponentAdmission.hpp"
#include "PanAlignerMinimizer.hpp"
#include "gfa-priv.h"

#include "bgzf.h"
#include "hts.h"
#include "tbx.h"
#include <zlib.h>

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <parallel/algorithm>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace pangenie_vcf {
namespace {

struct VariantRecord {
    std::string chromosome;
    uint64_t position = 0;
    std::string id;
    std::string reference;
    std::vector<std::string> alternatives;
};

struct ParsedVariant {
    VariantRecord record;
    std::vector<std::string> alleles;
    std::vector<uint16_t> path_alleles;
    uint16_t defined_alleles = 0;
};

using TbxHandle = std::unique_ptr<tbx_t, decltype(&tbx_destroy)>;
using BgzfHandle = std::unique_ptr<BGZF, decltype(&bgzf_close)>;
using IteratorHandle = std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)>;

std::vector<std::string> split(const std::string& text, char separator)
{
    std::vector<std::string> result;
    size_t begin = 0;
    for (;;) {
        const size_t end = text.find(separator, begin);
        result.push_back(text.substr(begin, end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return result;
}

bool getline_gz(gzFile stream, std::string& line)
{
    line.clear();
    char buffer[1 << 16];
    for (;;) {
        char* result = gzgets(stream, buffer, sizeof(buffer));
        if (result == nullptr) return !line.empty();
        line += buffer;
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
            if (!line.empty() && line.back() == '\r') line.pop_back();
            return true;
        }
        if (gzeof(stream)) return true;
    }
}

std::string sequence_name(const std::string& header)
{
    const size_t end = header.find_first_of(" \t");
    return header.substr(0, end);
}

std::map<std::string, std::string> read_fasta(const std::string& filename)
{
    gzFile stream = gzopen(filename.c_str(), "rb");
    if (stream == nullptr)
        throw std::runtime_error("PanGenieVcf: cannot open reference FASTA: " + filename);

    std::map<std::string, std::string> sequences;
    auto insert_sequence = [&](std::string sequence_name,
                               std::string sequence_bases) {
        const auto inserted = sequences.emplace(
            std::move(sequence_name), std::move(sequence_bases));
        if (!inserted.second)
            throw std::runtime_error(
                "PanGenieVcf: duplicate FASTA sequence name: " +
                inserted.first->first);
    };
    std::string name, sequence, line;
    try {
        while (getline_gz(stream, line)) {
            if (line.empty()) continue;
            if (line.front() == '>') {
                if (!name.empty()) insert_sequence(name, std::move(sequence));
                name = sequence_name(line.substr(1));
                sequence.clear();
                if (name.empty())
                    throw std::runtime_error(
                        "PanGenieVcf: empty FASTA sequence name");
            } else {
                if (name.empty())
                    throw std::runtime_error(
                        "PanGenieVcf: FASTA sequence before first header");
                for (char base : line)
                    if (!std::isspace(static_cast<unsigned char>(base)))
                        sequence.push_back(static_cast<char>(std::toupper(
                            static_cast<unsigned char>(base))));
            }
        }
        if (!name.empty()) insert_sequence(name, std::move(sequence));
    } catch (...) {
        gzclose(stream);
        throw;
    }
    gzclose(stream);
    if (sequences.empty())
        throw std::runtime_error("PanGenieVcf: reference FASTA contains no sequences");
    return sequences;
}

std::vector<std::string> tabix_sequence_names(tbx_t* tabix)
{
    int metadata_size = 0;
    const uint8_t* metadata = hts_idx_get_meta(tabix->idx, &metadata_size);
    if (metadata == nullptr || metadata_size < 28)
        throw std::runtime_error(
            "PanGenieVcf: tabix index has no sequence dictionary");
    uint32_t names_size = 0;
    std::memcpy(&names_size, metadata + 24, sizeof(names_size));
    if (names_size == 0 || names_size > static_cast<uint32_t>(metadata_size - 28))
        throw std::runtime_error(
            "PanGenieVcf: malformed tabix sequence dictionary");
    std::vector<std::string> names;
    const char* cursor = reinterpret_cast<const char*>(metadata + 28);
    const char* const end = cursor + names_size;
    while (cursor < end) {
        const void* terminator = std::memchr(
            cursor, '\0', static_cast<size_t>(end - cursor));
        if (terminator == nullptr || terminator == cursor)
            throw std::runtime_error(
                "PanGenieVcf: malformed tabix sequence name");
        const char* zero = static_cast<const char*>(terminator);
        names.emplace_back(cursor, zero);
        cursor = zero + 1;
    }
    if (cursor != end || names.empty())
        throw std::runtime_error(
            "PanGenieVcf: malformed tabix sequence dictionary extent");
    for (size_t tid = 0; tid < names.size(); ++tid)
        if (tbx_name2id(tabix, names[tid].c_str()) != static_cast<int>(tid))
            throw std::runtime_error(
                "PanGenieVcf: tabix sequence order is inconsistent");
    return names;
}

std::vector<ComponentSelection> select_components(
    const std::map<std::string, std::string>& reference,
    const std::vector<std::string>& chromosomes, const ComponentGate& gate)
{
    if (gate.read_spectrum_hashes == nullptr &&
        !gate.reference_support_query)
        throw std::invalid_argument(
            "PanGenieVcf: indexed component gate has no read-spectrum query");
    if (gate.read_spectrum_hashes != nullptr &&
        (!std::is_sorted(gate.read_spectrum_hashes->begin(),
                         gate.read_spectrum_hashes->end()) ||
         std::adjacent_find(gate.read_spectrum_hashes->begin(),
                            gate.read_spectrum_hashes->end()) !=
             gate.read_spectrum_hashes->end()))
        throw std::invalid_argument(
            "PanGenieVcf: indexed component gate requires a sorted distinct read spectrum");
    if (gate.k <= 0 || gate.k > 31 || gate.w <= 0 || gate.w >= 256)
        throw std::invalid_argument(
            "PanGenieVcf: indexed component gate has invalid minimizer parameters");
    if (chromosomes.size() >=
        static_cast<size_t>(std::numeric_limits<uint32_t>::max() - 1ULL))
        throw std::overflow_error(
            "PanGenieVcf: too many panel components for compact hash ownership");

    const auto reference_hashes = [&](uint32_t component) {
        const auto found = reference.find(chromosomes[component]);
        if (found == reference.end())
            throw std::runtime_error(
                "VCF chromosome is absent from reference FASTA: " +
                chromosomes[component]);
        if (found->second.empty())
            throw std::runtime_error(
                "VCF chromosome has an empty reference sequence: " +
                chromosomes[component]);
        if (found->second.size() > std::numeric_limits<uint32_t>::max())
            throw std::overflow_error(
                "PanGenieVcf: indexed panel reference exceeds the uint32_t minimizer coordinate space: " +
                chromosomes[component]);
        std::vector<uint64_t> hashes = dipgenie::sequence_minimizer_hashes(
            found->second, gate.w, gate.k);
        if (gate.num_threads > 1 && hashes.size() >= (1ULL << 20))
            __gnu_parallel::sort(hashes.begin(), hashes.end());
        else
            std::sort(hashes.begin(), hashes.end());
        hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
        return hashes;
    };

    std::vector<uint64_t> shared(chromosomes.size(), 0);
    std::vector<uint64_t> unique(chromosomes.size(), 0);
    if (gate.reference_support_query) {
        std::vector<std::vector<uint64_t>> signatures;
        signatures.reserve(chromosomes.size());
        for (uint32_t component = 0; component < chromosomes.size();
             ++component)
            signatures.push_back(reference_hashes(component));
        ReferenceSupportCensus census =
            gate.reference_support_query(signatures);
        if (census.read_shared_hashes.size() != chromosomes.size() ||
            census.component_unique_hashes.size() != chromosomes.size())
            throw std::runtime_error(
                "PanGenieVcf: read-spectrum query returned a malformed component census");
        shared = std::move(census.read_shared_hashes);
        unique = std::move(census.component_unique_hashes);
    } else {
        const auto& reads = *gate.read_spectrum_hashes;
        constexpr uint32_t unseen = std::numeric_limits<uint32_t>::max();
        constexpr uint32_t multi = unseen - 1U;
        std::vector<uint32_t> owner(reads.size(), unseen);
        for (uint32_t component = 0; component < chromosomes.size();
             ++component) {
            const std::vector<uint64_t> hashes =
                reference_hashes(component);
            size_t reference_index = 0;
            size_t read_index = 0;
            while (reference_index < hashes.size() &&
                   read_index < reads.size()) {
                if (hashes[reference_index] < reads[read_index]) {
                    ++reference_index;
                } else if (reads[read_index] < hashes[reference_index]) {
                    ++read_index;
                } else {
                    ++shared[component];
                    uint32_t& hash_owner = owner[read_index];
                    if (hash_owner == unseen)
                        hash_owner = component;
                    else if (hash_owner != component)
                        hash_owner = multi;
                    ++reference_index;
                    ++read_index;
                }
            }
        }
        for (uint32_t hash_owner : owner)
            if (hash_owner < unique.size()) ++unique[hash_owner];
    }
    for (size_t component = 0; component < chromosomes.size(); ++component)
        if (unique[component] > shared[component])
            throw std::runtime_error(
                "PanGenieVcf: component-unique read support exceeds shared support");
    const uint64_t maximum_unique = unique.empty()
        ? 0 : *std::max_element(unique.begin(), unique.end());

    std::vector<ComponentSelection> selection;
    selection.reserve(chromosomes.size());
    for (uint32_t component = 0; component < chromosomes.size(); ++component) {
        const uint64_t span = reference.at(chromosomes[component]).size();
        const auto threshold =
            dipgenie::component_admission::make_preflight_threshold(
                maximum_unique, gate.minimum_anchor_hashes,
                gate.minimum_dominance_basis_points, span,
                gate.minimum_reference_basis_points);
        const bool admitted = dipgenie::component_admission::passes_preflight(
            unique[component], threshold);
        selection.push_back(
            {component, chromosomes[component], span, shared[component],
             unique[component], maximum_unique, threshold.absolute,
             threshold.dominance, threshold.reference_span,
             threshold.effective, admitted});
        fprintf(stderr,
                "[M::ComponentAdmission::preflight] id=%u sequence=%s "
                "reference_read_shared_hashes=%llu "
                "reference_component_unique_hashes=%llu "
                "cross_component_shared_hashes=%llu "
                "maximum_component_unique_hashes=%llu "
                "absolute_threshold=%llu dominance_basis_points=%u "
                "dominance_threshold=%llu admission_span=%llu "
                "admission_span_source=reference reference_basis_points=%u "
                "reference_threshold=%llu effective_threshold=%llu "
                "threshold_reached=%s decision=%s "
                "phase=indexed-reference-FASTA\n",
                component, chromosomes[component].c_str(),
                static_cast<unsigned long long>(shared[component]),
                static_cast<unsigned long long>(unique[component]),
                static_cast<unsigned long long>(
                    shared[component] - unique[component]),
                static_cast<unsigned long long>(maximum_unique),
                static_cast<unsigned long long>(threshold.absolute),
                gate.minimum_dominance_basis_points,
                static_cast<unsigned long long>(threshold.dominance),
                static_cast<unsigned long long>(span),
                gate.minimum_reference_basis_points,
                static_cast<unsigned long long>(threshold.reference_span),
                static_cast<unsigned long long>(threshold.effective),
                admitted ? "yes" : "no", admitted ? "continue" : "skip");
        fprintf(stderr,
                "[M::PanGenieVcf::preselection-component] id=%u "
                "chromosome=%s reference_bp=%llu "
                "reference_read_shared_hashes=%llu "
                "reference_component_unique_hashes=%llu "
                "cross_component_shared_hashes=%llu "
                "maximum_component_unique_hashes=%llu "
                "absolute_threshold=%llu dominance_basis_points=%u "
                "dominance_threshold=%llu reference_basis_points=%u "
                "reference_threshold=%llu effective_threshold=%llu "
                "decision=%s provenance=global-reference-ownership\n",
                component, chromosomes[component].c_str(),
                static_cast<unsigned long long>(span),
                static_cast<unsigned long long>(shared[component]),
                static_cast<unsigned long long>(unique[component]),
                static_cast<unsigned long long>(
                    shared[component] - unique[component]),
                static_cast<unsigned long long>(maximum_unique),
                static_cast<unsigned long long>(threshold.absolute),
                gate.minimum_dominance_basis_points,
                static_cast<unsigned long long>(threshold.dominance),
                gate.minimum_reference_basis_points,
                static_cast<unsigned long long>(threshold.reference_span),
                static_cast<unsigned long long>(threshold.effective),
                admitted ? "materialize" : "skip");
    }
    return selection;
}

bool explicit_dna(const std::string& sequence)
{
    if (sequence.empty()) return false;
    for (char base : sequence) {
        switch (std::toupper(static_cast<unsigned char>(base))) {
            case 'A': case 'C': case 'G': case 'T': break;
            default: return false;
        }
    }
    return true;
}

std::string uppercase(std::string sequence)
{
    for (char& base : sequence)
        base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    return sequence;
}

uint64_t parse_position(const std::string& token)
{
    errno = 0;
    char* end = nullptr;
    const unsigned long long value = std::strtoull(token.c_str(), &end, 10);
    if (errno != 0 || end == token.c_str() || *end != 0 || value == 0)
        throw std::runtime_error("PanGenieVcf: invalid one-based VCF position: " + token);
    return static_cast<uint64_t>(value - 1);
}

uint32_t add_segment(gfa_t* graph, const std::string& name,
                     const std::string& sequence, int32_t stable_name,
                     uint64_t offset)
{
    if (sequence.empty())
        throw std::runtime_error("PanGenieVcf: attempted to create an empty segment");
    if (sequence.size() > static_cast<size_t>(std::numeric_limits<int32_t>::max()))
        throw std::runtime_error("PanGenieVcf: segment exceeds INT32_MAX bases");
    const int64_t id = gfa_add_seg(graph, name.c_str());
    if (id < 0 || static_cast<uint64_t>(id) > UINT32_MAX)
        throw std::runtime_error("PanGenieVcf: component exceeds UINT32_MAX segments");
    gfa_seg_t& segment = graph->seg[id];
    segment.seq = gfa_strdup(sequence.c_str());
    segment.len = static_cast<int32_t>(sequence.size());
    segment.snid = stable_name;
    segment.soff = offset <= static_cast<uint64_t>(INT32_MAX)
        ? static_cast<int32_t>(offset) : -1;
    segment.rank = 0;
    gfa_sseq_update(graph, &segment);
    return static_cast<uint32_t>(id);
}

void connect(gfa_t* graph, const std::vector<uint32_t>& from,
             const std::vector<uint32_t>& to)
{
    for (uint32_t left : from)
        for (uint32_t right : to)
            gfa_add_arc1(graph, static_cast<uint64_t>(left) << 1,
                         static_cast<uint64_t>(right) << 1, 0, 0, -1, 0);
}

void append_all(std::vector<std::vector<uint64_t>>& paths, uint32_t segment)
{
    for (auto& path : paths) path.push_back(static_cast<uint64_t>(segment) << 1);
}

void add_walk(gfa_t* graph, const std::string& sample, int32_t haplotype,
              int32_t stable_name, uint64_t reference_length,
              const std::vector<uint64_t>& path)
{
    if (path.empty()) throw std::runtime_error("PanGenieVcf: empty donor path");
    if (graph->n_walk == graph->m_walk) {
        if (graph->m_walk == UINT32_MAX)
            throw std::overflow_error(
                "PanGenieVcf: graph exceeds UINT32_MAX walks");
        const uint64_t grown = graph->m_walk
            ? static_cast<uint64_t>(graph->m_walk) +
                  std::max<uint32_t>(1, graph->m_walk >> 1)
            : 16;
        graph->m_walk = static_cast<uint32_t>(
            std::min<uint64_t>(grown, UINT32_MAX));
        void* allocation = std::realloc(
            graph->walk,
            static_cast<size_t>(graph->m_walk) * sizeof(gfa_walk_t));
        if (allocation == nullptr)
            throw std::bad_alloc();
        graph->walk = static_cast<gfa_walk_t*>(allocation);
    }
    gfa_walk_t& walk = graph->walk[graph->n_walk++];
    std::memset(&walk, 0, sizeof(walk));
    walk.sample = gfa_strdup(sample.c_str());
    walk.snid = stable_name;
    walk.hap = haplotype;
    walk.st = 0;
    walk.en = reference_length <= static_cast<uint64_t>(INT64_MAX)
        ? static_cast<int64_t>(reference_length) : -1;
    walk.n_v = path.size();
    walk.v = static_cast<uint64_t*>(std::malloc(path.size() * sizeof(uint64_t)));
    if (walk.v == nullptr) throw std::bad_alloc();
    std::copy(path.begin(), path.end(), walk.v);
}

struct ContigBuilder {
    PanelContig metadata;
    const std::string* reference = nullptr;
    int32_t stable_name = -1;
    std::string segment_prefix;
    std::vector<std::vector<uint64_t>> paths;
    std::vector<uint32_t> frontier;
    uint64_t previous_retained_end = 0;
    uint64_t reference_serial = 0;
    uint64_t variant_serial = 0;
    bool has_retained_record = false;
};

std::unique_ptr<ContigBuilder> begin_contig(
    gfa_t* graph, const std::string& chromosome, const std::string& reference,
    size_t sample_count, size_t contig_id)
{
    if (reference.empty())
        throw std::runtime_error(
            "PanGenieVcf: VCF chromosome has an empty reference sequence: " +
            chromosome);
    if (sample_count > (static_cast<size_t>(UINT32_MAX) - 1) / 2)
        throw std::overflow_error(
            "PanGenieVcf: retained sample count exceeds walk capacity");
    const size_t path_count = 1 + 2 * sample_count;
    if (path_count > static_cast<size_t>(UINT32_MAX - graph->n_walk))
        throw std::overflow_error(
            "PanGenieVcf: multi-contig panel exceeds UINT32_MAX walks");

    auto builder = std::make_unique<ContigBuilder>();
    builder->metadata.chromosome = chromosome;
    builder->metadata.reference_length = reference.size();
    builder->reference = &reference;
    builder->stable_name = gfa_sseq_add(graph, chromosome.c_str());
    if (builder->stable_name < 0)
        throw std::runtime_error(
            "PanGenieVcf: cannot register reference sequence: " + chromosome);
    builder->segment_prefix = "pg_c" + std::to_string(contig_id) + "_";
    builder->paths.resize(path_count);
    return builder;
}

void append_variant(gfa_t* graph, ContigBuilder& builder,
                    const ParsedVariant& variant)
{
    const std::string& ref = *builder.reference;
    if (variant.record.position > builder.previous_retained_end) {
        const std::string sequence = ref.substr(
            builder.previous_retained_end,
            variant.record.position - builder.previous_retained_end);
        const uint32_t segment = add_segment(
            graph,
            builder.segment_prefix + "ref_" +
                std::to_string(builder.reference_serial++),
            sequence, builder.stable_name, builder.previous_retained_end);
        connect(graph, builder.frontier, {segment});
        builder.frontier = {segment};
        append_all(builder.paths, segment);
    }

    std::vector<uint32_t> allele_segments;
    allele_segments.reserve(variant.alleles.size());
    for (size_t allele = 0; allele < variant.alleles.size(); ++allele) {
        allele_segments.push_back(add_segment(
            graph,
            builder.segment_prefix + "var_" +
                std::to_string(builder.variant_serial) + "_" +
                std::to_string(allele),
            variant.alleles[allele], builder.stable_name,
            variant.record.position));
    }
    ++builder.variant_serial;
    connect(graph, builder.frontier, allele_segments);
    for (size_t path = 0; path < builder.paths.size(); ++path) {
        const uint16_t allele = variant.path_alleles.at(path);
        builder.paths[path].push_back(
            static_cast<uint64_t>(allele_segments.at(allele)) << 1);
    }
    builder.frontier = std::move(allele_segments);
    builder.previous_retained_end = variant.record.position +
                                    variant.record.reference.size();
    builder.has_retained_record = true;
}

void finish_contig(gfa_t* graph, ContigBuilder& builder,
                   const std::vector<std::string>& samples)
{
    const std::string& ref = *builder.reference;
    if (builder.previous_retained_end < ref.size()) {
        const uint32_t segment = add_segment(
            graph,
            builder.segment_prefix + "ref_" +
                std::to_string(builder.reference_serial++),
            ref.substr(builder.previous_retained_end), builder.stable_name,
            builder.previous_retained_end);
        connect(graph, builder.frontier, {segment});
        append_all(builder.paths, segment);
        builder.frontier = {segment};
    }
    if (builder.frontier.empty()) {
        const uint32_t segment = add_segment(
            graph, builder.segment_prefix + "ref_0", ref,
            builder.stable_name, 0);
        append_all(builder.paths, segment);
    }

    add_walk(graph, "REF", 0, builder.stable_name, ref.size(),
             builder.paths[0]);
    for (size_t sample = 0; sample < samples.size(); ++sample) {
        add_walk(graph, samples[sample], 1, builder.stable_name, ref.size(),
                 builder.paths[1 + 2 * sample]);
        add_walk(graph, samples[sample], 2, builder.stable_name, ref.size(),
                 builder.paths[2 + 2 * sample]);
    }
}

class PanelLineReader {
public:
    // This embedded HTSlib's TBI index has a 2^29-coordinate address space.
    static constexpr int TABIX_COORDINATE_LIMIT = 1 << 29;

    PanelLineReader(gzFile sequential, const std::string& filename,
                    const std::vector<ComponentSelection>* selection)
        : sequential_(sequential), selection_(selection),
          tabix_(nullptr, &tbx_destroy), data_(nullptr, &bgzf_close),
          iterator_(nullptr, &hts_itr_destroy)
    {
        if (selection_ == nullptr) return;
        tabix_.reset(tbx_index_load(filename.c_str()));
        data_.reset(bgzf_open(filename.c_str(), "r"));
        if (!tabix_ || !data_)
            throw std::runtime_error(
                "PanGenieVcf: tabix index became unavailable during indexed import");
    }

    ~PanelLineReader() { std::free(record_.s); }

    bool getline(std::string& line)
    {
        if (selection_ == nullptr) return getline_gz(sequential_, line);
        if (!header_complete_) {
            if (!getline_gz(sequential_, line))
                throw std::runtime_error(
                    "PanGenieVcf: missing VCF header before indexed records");
            if (line.rfind("#CHROM", 0) == 0) header_complete_ = true;
            return true;
        }
        for (;;) {
            if (iterator_) {
                const int status = tbx_itr_next(
                    data_.get(), tabix_.get(), iterator_.get(), &record_);
                if (status >= 0) {
                    line.assign(record_.s, record_.l);
                    if (!line.empty() && line.back() == '\r') line.pop_back();
                    return true;
                }
                if (status < -1)
                    throw std::runtime_error(
                        "PanGenieVcf: tabix record query failed");
                iterator_.reset();
            }
            while (selection_index_ < selection_->size() &&
                   !selection_->at(selection_index_).admitted)
                ++selection_index_;
            if (selection_index_ == selection_->size()) return false;
            const ComponentSelection& component =
                selection_->at(selection_index_++);
            const int tid = tbx_name2id(
                tabix_.get(), component.chromosome.c_str());
            if (tid < 0 || static_cast<uint32_t>(tid) !=
                    component.input_component_id)
                throw std::runtime_error(
                    "PanGenieVcf: selected chromosome is absent from the tabix index: " +
                    component.chromosome);
            iterator_.reset(tbx_itr_queryi(
                tabix_.get(), tid, 0, TABIX_COORDINATE_LIMIT));
            if (!iterator_)
                throw std::runtime_error(
                    "PanGenieVcf: cannot query selected chromosome from tabix index: " +
                    component.chromosome);
        }
    }

private:
    gzFile sequential_ = nullptr;
    const std::vector<ComponentSelection>* selection_ = nullptr;
    TbxHandle tabix_;
    BgzfHandle data_;
    IteratorHandle iterator_;
    kstring_t record_{0, 0, nullptr};
    size_t selection_index_ = 0;
    bool header_complete_ = false;
};

}  // namespace

TabixStatus inspect_tabix(const std::string& vcf_filename)
{
    if (vcf_filename.empty() || vcf_filename == "-")
        return {false, "non-seekable-input"};
    std::error_code filesystem_error;
    if (!std::filesystem::is_regular_file(
            vcf_filename + ".tbi", filesystem_error))
        return {false, "tbi-index-unavailable"};
    TbxHandle tabix(tbx_index_load(vcf_filename.c_str()), &tbx_destroy);
    if (!tabix) return {false, "tbi-index-unreadable"};
    if (tabix->conf.preset != TBX_VCF)
        return {false, "tbi-index-not-vcf"};
    BgzfHandle data(bgzf_open(vcf_filename.c_str(), "r"), &bgzf_close);
    if (!data) return {false, "panel-vcf-unreadable"};
    char byte = 0;
    if (bgzf_read(data.get(), &byte, 1) != 1 || data->errcode != 0)
        return {false, "panel-vcf-not-bgzf"};
    try {
        if (tabix_sequence_names(tabix.get()).empty())
            return {false, "tbi-index-empty"};
    } catch (...) {
        return {false, "tbi-index-malformed"};
    }
    return {true, "usable-tbi"};
}

bool component_gate_enabled(uint64_t minimum_anchor_hashes,
                            uint32_t minimum_dominance_basis_points,
                            uint32_t minimum_reference_basis_points)
{
    // Dominance is coupled to the absolute floor by the established
    // threshold formula; it cannot activate filtering by itself.
    (void)minimum_dominance_basis_points;
    return minimum_anchor_hashes != 0 || minimum_reference_basis_points != 0;
}

PanelGraph read_panel(const std::string& vcf_filename,
                      const std::string& reference_filename,
                      const std::unordered_set<std::string>& excluded_samples,
                      const ComponentGate* component_gate)
{
    PanelGraph result;
    const auto reference = read_fasta(reference_filename);

    if (component_gate != nullptr) {
        TbxHandle tabix(tbx_index_load(vcf_filename.c_str()), &tbx_destroy);
        if (!tabix || tabix->conf.preset != TBX_VCF)
            throw std::runtime_error(
                "PanGenieVcf: indexed preselection requested without a usable VCF tabix index");
        const std::vector<std::string> chromosomes =
            tabix_sequence_names(tabix.get());
        result.component_selection = select_components(
            reference, chromosomes, *component_gate);
        const size_t admitted = static_cast<size_t>(std::count_if(
            result.component_selection.begin(), result.component_selection.end(),
            [](const ComponentSelection& component) {
                return component.admitted;
            }));
        if (admitted == 0)
            throw std::runtime_error(
                "PanGenieVcf: no panel component passed indexed read-support preselection");
        result.selection_mode = "indexed-tabix";
    }

    gzFile stream = gzopen(vcf_filename.c_str(), "rb");
    if (stream == nullptr)
        throw std::runtime_error("PanGenieVcf: cannot open panel VCF: " + vcf_filename);

    std::unique_ptr<gfa_t, decltype(&gfa_destroy)> graph(
        gfa_init(), &gfa_destroy);
    if (!graph) {
        gzclose(stream);
        throw std::bad_alloc();
    }
    if (!result.component_selection.empty()) {
        for (const ComponentSelection& component : result.component_selection) {
            const int32_t stable_name =
                gfa_sseq_add(graph.get(), component.chromosome.c_str());
            if (stable_name < 0 || static_cast<uint32_t>(stable_name) !=
                    component.input_component_id) {
                gzclose(stream);
                throw std::runtime_error(
                    "PanGenieVcf: indexed stable-sequence order disagrees with the panel");
            }
        }
    }
    std::unique_ptr<ContigBuilder> contig;
    std::unordered_set<std::string> finished_contigs;
    std::vector<std::string> input_samples;
    std::vector<size_t> retained_sample_columns;
    std::unordered_set<std::string> seen_samples;
    std::unordered_set<std::string> seen_exclusions;
    std::string line;
    bool saw_header = false;
    size_t line_number = 0;
    std::unique_ptr<PanelLineReader> line_reader;
    try {
        line_reader = std::make_unique<PanelLineReader>(
            stream, vcf_filename,
            component_gate == nullptr ? nullptr : &result.component_selection);
    } catch (...) {
        gzclose(stream);
        throw;
    }

    auto flush_contig = [&]() {
        if (!contig) return;
        finish_contig(graph.get(), *contig, result.samples);
        const PanelContig metadata = contig->metadata;
        fprintf(stderr,
                "[M::PanGenieVcf::component] chromosome=%s reference_bp=%llu "
                "samples=%zu donor_paths=%zu input_variants=%llu "
                "retained_variants=%llu skipped_variants=%llu\n",
                metadata.chromosome.c_str(),
                static_cast<unsigned long long>(metadata.reference_length),
                result.samples.size(), contig->paths.size(),
                static_cast<unsigned long long>(metadata.input_variants),
                static_cast<unsigned long long>(metadata.retained_variants),
                static_cast<unsigned long long>(metadata.skipped_variants));
        finished_contigs.insert(metadata.chromosome);
        result.contigs.push_back(metadata);
        contig.reset();
    };

    try {
        while (line_reader->getline(line)) {
            ++line_number;
            if (line.empty() || line.rfind("##", 0) == 0) continue;
            const std::vector<std::string> fields = split(line, '\t');
            if (!line.empty() && line.front() == '#') {
                if (saw_header)
                    throw std::runtime_error("duplicate VCF column header");
                static const char* required[] = {
                    "#CHROM", "POS", "ID", "REF", "ALT",
                    "QUAL", "FILTER", "INFO", "FORMAT"};
                if (fields.size() < 10)
                    throw std::runtime_error("panel VCF has no phased donor samples");
                for (size_t i = 0; i < 9; ++i)
                    if (fields[i] != required[i])
                        throw std::runtime_error("malformed VCF header");
                input_samples.assign(fields.begin() + 9, fields.end());
                result.input_samples = input_samples.size();
                for (size_t sample = 0; sample < input_samples.size(); ++sample) {
                    const std::string& name = input_samples[sample];
                    if (name.empty() || !seen_samples.insert(name).second)
                        throw std::runtime_error(
                            "empty or duplicate panel sample name: " + name);
                    if (excluded_samples.count(name) != 0) {
                        seen_exclusions.insert(name);
                        ++result.excluded_samples;
                    } else {
                        retained_sample_columns.push_back(sample);
                        result.samples.push_back(name);
                    }
                }
                if (result.samples.size() >
                    (static_cast<size_t>(UINT32_MAX) - 1) / 2)
                    throw std::runtime_error(
                        "PanGenieVcf: retained sample count exceeds walk capacity");
                if (result.samples.empty())
                    throw std::runtime_error(
                        "panel VCF has no samples after exclusions");
                for (const std::string& excluded : excluded_samples)
                    if (seen_exclusions.count(excluded) == 0)
                        throw std::runtime_error(
                            "excluded panel sample is absent from the VCF: " +
                            excluded);
                saw_header = true;
                continue;
            }
            ++result.input_variants;
            if (!saw_header || fields.size() != 9 + input_samples.size())
                throw std::runtime_error("malformed VCF record or inconsistent sample count");

            const auto reference_it = reference.find(fields[0]);
            if (reference_it == reference.end())
                throw std::runtime_error(
                    "VCF chromosome is absent from reference FASTA: " + fields[0]);
            if (reference_it->second.empty())
                throw std::runtime_error(
                    "VCF chromosome has an empty reference sequence: " + fields[0]);
            if (!contig || contig->metadata.chromosome != fields[0]) {
                flush_contig();
                if (finished_contigs.count(fields[0]) != 0)
                    throw std::runtime_error(
                        "VCF chromosome block reappears after completion: " +
                        fields[0]);
                size_t input_contig_id = result.contigs.size();
                if (!result.component_selection.empty()) {
                    const auto selected = std::find_if(
                        result.component_selection.begin(),
                        result.component_selection.end(),
                        [&](const ComponentSelection& component) {
                            return component.chromosome == fields[0];
                        });
                    if (selected == result.component_selection.end() ||
                        !selected->admitted)
                        throw std::runtime_error(
                            "tabix returned a chromosome outside indexed component preselection: " +
                            fields[0]);
                    input_contig_id = selected->input_component_id;
                }
                contig = begin_contig(
                    graph.get(), fields[0], reference_it->second,
                    result.samples.size(), input_contig_id);
            }
            ++contig->metadata.input_variants;

            ParsedVariant variant;
            variant.record.chromosome = fields[0];
            variant.record.position = parse_position(fields[1]);
            variant.record.id = fields[2];
            variant.record.reference = uppercase(fields[3]);
            variant.record.alternatives = split(fields[4], ',');
            for (std::string& allele : variant.record.alternatives)
                allele = uppercase(std::move(allele));

            if (variant.record.reference.size() >
                std::numeric_limits<uint64_t>::max() - variant.record.position)
                throw std::runtime_error("VCF REF interval overflows uint64_t");
            const uint64_t end = variant.record.position +
                                 variant.record.reference.size();
            if (end > reference_it->second.size())
                throw std::runtime_error("VCF REF allele extends beyond reference");
            if (contig->has_retained_record &&
                variant.record.position < contig->previous_retained_end)
                throw std::runtime_error(
                    "PanGenieVcf: overlapping or unsorted VCF records on " +
                    fields[0] + " do not define independent bubbles");
            if (uppercase(reference_it->second.substr(
                    variant.record.position, variant.record.reference.size())) !=
                variant.record.reference)
                throw std::runtime_error(
                    "PanGenieVcf: VCF REF allele does not match the reference FASTA");

            bool defined = explicit_dna(variant.record.reference);
            for (const auto& allele : variant.record.alternatives)
                defined = defined && explicit_dna(allele);
            if (!defined) {
                ++result.skipped_variants;
                ++contig->metadata.skipped_variants;
                continue;
            }

            variant.alleles.push_back(variant.record.reference);
            variant.alleles.insert(variant.alleles.end(),
                                   variant.record.alternatives.begin(),
                                   variant.record.alternatives.end());
            if (variant.alleles.size() > 65534)
                throw std::runtime_error("PanGenieVcf: more than 65,534 alleles at one variant");
            variant.defined_alleles = static_cast<uint16_t>(variant.alleles.size());

            const std::vector<std::string> format = split(fields[8], ':');
            const auto gt_it = std::find(format.begin(), format.end(), "GT");
            if (gt_it == format.end())
                throw std::runtime_error("PanGenieVcf: GT is absent from FORMAT");
            const size_t gt_index = static_cast<size_t>(gt_it - format.begin());
            variant.path_alleles.reserve(1 + 2 * result.samples.size());
            variant.path_alleles.push_back(0);  // added reference path
            uint16_t undefined_index = variant.defined_alleles;
            bool undefined_added = false;
            for (size_t retained = 0;
                 retained < retained_sample_columns.size(); ++retained) {
                const size_t sample = retained_sample_columns[retained];
                const std::vector<std::string> sample_fields = split(
                    fields[9 + sample], ':');
                if (gt_index >= sample_fields.size())
                    throw std::runtime_error("PanGenieVcf: sample has no GT field");
                const std::string& genotype = sample_fields[gt_index];
                if (genotype.find('/') != std::string::npos ||
                    genotype.find('|') == std::string::npos)
                    throw std::runtime_error("PanGenieVcf: every panel genotype must be phased");
                const std::vector<std::string> alleles = split(genotype, '|');
                if (alleles.size() != 2)
                    throw std::runtime_error("PanGenieVcf: panel genotypes must be diploid");
                for (const std::string& allele : alleles) {
                    if (allele == ".") {
                        if (!undefined_added) {
                            variant.alleles.push_back("N");
                            undefined_added = true;
                        }
                        variant.path_alleles.push_back(undefined_index);
                        continue;
                    }
                    errno = 0;
                    char* allele_end = nullptr;
                    const unsigned long value = std::strtoul(allele.c_str(), &allele_end, 10);
                    if (errno != 0 || allele_end == allele.c_str() || *allele_end != 0 ||
                        value >= variant.defined_alleles)
                        throw std::runtime_error("PanGenieVcf: invalid genotype allele index");
                    variant.path_alleles.push_back(static_cast<uint16_t>(value));
                }
            }
            append_variant(graph.get(), *contig, variant);
            ++contig->metadata.retained_variants;
            ++result.retained_variants;
        }
        flush_contig();
    } catch (const std::exception& error) {
        gzclose(stream);
        throw std::runtime_error("PanGenieVcf line " + std::to_string(line_number) +
                                 ": " + error.what());
    }
    gzclose(stream);
    if (!saw_header) throw std::runtime_error("PanGenieVcf: missing VCF header");
    if (result.contigs.empty())
        throw std::runtime_error("PanGenieVcf: panel has no VCF records");
    if (result.retained_variants == 0)
        throw std::runtime_error("PanGenieVcf: panel has no usable variants");
    gfa_finalize(graph.get());

    if (result.contigs.size() == 1) {
        result.chromosome = result.contigs.front().chromosome;
        result.reference_length = result.contigs.front().reference_length;
    }
    fprintf(stderr,
            "[M::PanGenieVcf::read] selection_mode=%s panel_contigs=%zu "
            "contigs=%zu input_samples=%llu "
            "excluded_samples=%llu retained_samples=%zu donor_paths=%u "
            "input_variants=%llu retained_variants=%llu skipped_variants=%llu "
            "segments=%llu arcs=%llu\n",
            result.selection_mode.c_str(),
            result.component_selection.empty()
                ? result.contigs.size() : result.component_selection.size(),
            result.contigs.size(),
            static_cast<unsigned long long>(result.input_samples),
            static_cast<unsigned long long>(result.excluded_samples),
            result.samples.size(), graph->n_walk,
            static_cast<unsigned long long>(result.input_variants),
            static_cast<unsigned long long>(result.retained_variants),
            static_cast<unsigned long long>(result.skipped_variants),
            static_cast<unsigned long long>(graph->n_seg),
            static_cast<unsigned long long>(graph->n_arc));
    if (!result.component_selection.empty()) {
        const size_t admitted = static_cast<size_t>(std::count_if(
            result.component_selection.begin(), result.component_selection.end(),
            [](const ComponentSelection& component) {
                return component.admitted;
            }));
        fprintf(stderr,
                "[M::PanGenieVcf::preselection] mode=indexed-tabix "
                "index=%s.tbi gate=enabled panel_contigs=%zu "
                "admitted_contigs=%zu skipped_contigs=%zu "
                "materialized_records=%llu marker=end\n",
                vcf_filename.c_str(), result.component_selection.size(),
                admitted, result.component_selection.size() - admitted,
                static_cast<unsigned long long>(result.input_variants));
    }
    result.graph = graph.release();
    return result;
}

}  // namespace pangenie_vcf
