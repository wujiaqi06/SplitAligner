// Independent bounded C ABI reuse benchmark. No CLI parsing or file output is timed.
// Run only the final Release candidate on the coordinator's pinned logical CPU.
#include "splitaligner/fast_c_api.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
constexpr size_t kTaxa = 302;
constexpr size_t kRepeats = 3;
constexpr size_t kCalls = 10000;

[[noreturn]] void fail(const char *operation) {
    std::fprintf(stderr, "memory_benchmark: %s: %s\n", operation, saf_last_error());
    std::exit(1);
}
void check(saf_status status, const char *operation) {
    if (status != SAF_OK) fail(operation);
}
double seconds(Clock::duration elapsed) {
    return std::chrono::duration<double>(elapsed).count();
}
struct Handles {
    saf_reference *reference = nullptr;
    saf_coverage *coverage = nullptr;
    saf_mapping *mapping = nullptr;
    ~Handles() {
        saf_mapping_destroy(mapping);
        saf_coverage_destroy(coverage);
        saf_reference_destroy(reference);
    }
};
// Deterministic rooted balanced binary array fixture; IDs are preorder node IDs.
size_t append_subtree(std::vector<int64_t> &parent, std::vector<int64_t> &taxon,
                      int64_t ancestor, size_t first, size_t last) {
    size_t node = parent.size();
    parent.push_back(ancestor);
    taxon.push_back(-1);
    if (last - first == 1) {
        taxon[node] = static_cast<int64_t>(first);
    } else {
        size_t middle = first + (last - first) / 2;
        append_subtree(parent, taxon, static_cast<int64_t>(node), first, middle);
        append_subtree(parent, taxon, static_cast<int64_t>(node), middle, last);
    }
    return node;
}
struct Result {
    double total_seconds = 0;
    double apply_seconds = 0;
    double checksum = 0;
    double expected_checksum = 0;
};
}

int main(int argc, char **argv) {
    if (argc != 1) {
        if (argc == 2 && std::strcmp(argv[1], "--help") == 0) {
            std::puts("memory_benchmark: fixed balanced302 full-coverage synthetic fixture; "
                      "3 repetitions x 10000 C ABI calls; all three output families.");
            return 0;
        }
        std::fputs("usage: memory_benchmark [--help]\n", stderr);
        return 2;
    }
    const auto preparation_start = Clock::now();
    Handles handles;
    std::vector<int64_t> parent, taxon;
    parent.reserve(2 * kTaxa - 1);
    taxon.reserve(2 * kTaxa - 1);
    append_subtree(parent, taxon, -1, 0, kTaxa);
    const auto fixture_ready = Clock::now();
    check(saf_reference_create(kTaxa, "memory-benchmark-balanced302-ids-0-through-301-v1",
                               parent.data(), taxon.data(), parent.size(), &handles.reference),
          "reference_create");
    std::vector<uint64_t> mask((kTaxa + 63) / 64, UINT64_MAX);
    mask.back() = (UINT64_C(1) << (kTaxa % 64)) - 1;
    check(saf_coverage_create(handles.reference, mask.data(), mask.size(), &handles.coverage),
          "coverage_create");
    check(saf_mapping_create(handles.coverage, parent.data(), taxon.data(), parent.size(),
                            &handles.mapping), "mapping_create");
    const auto plans_ready = Clock::now();
    saf_mapping_info info{};
    check(saf_mapping_get_info(handles.mapping, &info), "mapping_get_info");
    if (info.n_axes != 601 || info.n_projections != 601 || info.n_gene_splits != 601 ||
        info.n_reference_absent != 0 || info.n_gene_unmatched != 0) {
        std::fputs("memory_benchmark: unexpected synthetic fixture shape\n", stderr);
        return 1;
    }
    std::vector<double> lengths(parent.size());
    for (size_t i = 1; i < lengths.size(); ++i)
        lengths[i] = 0.125 + static_cast<double>(i % 127) / 1024.0;
    lengths[0] = std::numeric_limits<double>::quiet_NaN(); // ignored root slot
    std::vector<double> primitive(info.n_axes), projection(info.n_projections), gene(info.n_gene_splits);
    std::vector<uint8_t> primitive_valid(info.n_axes), projection_valid(info.n_projections), gene_valid(info.n_gene_splits);
    const saf_value_buffer outputs[3] = {
        {primitive.data(), primitive_valid.data(), primitive.size()},
        {projection.data(), projection_valid.data(), projection.size()},
        {gene.data(), gene_valid.data(), gene.size()}
    };
    // All raw states are D: each source appears once per output family, including
    // both halves of the root edge. This independent scalar invariant checks use
    // of changed inputs and consumption of every output, without a second matcher.
    double fixed_source_sum = 0;
    for (size_t i = 2; i < lengths.size(); ++i) fixed_source_sum += lengths[i];
    const auto preparation_end = Clock::now();
    Result results[kRepeats];
    for (size_t repeat = 0; repeat < kRepeats; ++repeat) {
        const auto loop_start = Clock::now();
        Clock::duration apply_duration = Clock::duration::zero();
        double checksum = 0, expected_checksum = 0;
        for (size_t iteration = 0; iteration < kCalls; ++iteration) {
            lengths[1] = 0.25 + static_cast<double>((repeat * kCalls + iteration) % 1024) / 4096.0;
            const auto apply_start = Clock::now();
            saf_status status = saf_mapping_apply(handles.mapping, lengths.data(), nullptr,
                                                 lengths.size(), outputs[0], outputs[1], outputs[2]);
            const auto apply_end = Clock::now();
            apply_duration += apply_end - apply_start;
            check(status, "mapping_apply");
            for (const auto &output : outputs) {
                for (size_t i = 0; i < output.capacity; ++i) {
                    if (output.valid[i] != 1 || !std::isfinite(output.values[i])) {
                        std::fputs("memory_benchmark: unexpected unavailable output\n", stderr);
                        return 1;
                    }
                    checksum += output.values[i];
                }
            }
            expected_checksum += 3.0 * (fixed_source_sum + lengths[1]);
        }
        const auto loop_end = Clock::now();
        if (checksum != expected_checksum) {
            std::fprintf(stderr, "memory_benchmark: checksum mismatch %.17g != %.17g\n",
                         checksum, expected_checksum);
            return 1;
        }
        results[repeat] = {seconds(loop_end - loop_start), seconds(apply_duration),
                           checksum, expected_checksum};
    }
    // All formatting/I/O happens after measured runs. Capture stdout plus GNU
    // time, affinity, host load and the executable SHA in the coordinator ledger.
    std::printf("fixture\tbalanced302_full_coverage_arrays_v1\n");
    std::printf("abi\t%s\n", saf_version());
    std::printf("n_taxa\t%zu\nn_nodes\t%zu\nn_axes\t%zu\nn_projections\t%zu\nn_gene_splits\t%zu\n",
                kTaxa, lengths.size(), info.n_axes, info.n_projections, info.n_gene_splits);
    std::printf("requested_output_slots_per_call\t%zu\noutput_bytes_per_call\t%zu\n",
                info.n_axes + info.n_projections + info.n_gene_splits,
                (info.n_axes + info.n_projections + info.n_gene_splits) * (sizeof(double) + sizeof(uint8_t)));
    std::printf("fixture_prepare_seconds\t%.9f\nplans_prepare_seconds\t%.9f\nbuffer_prepare_seconds\t%.9f\nprepare_total_seconds\t%.9f\n",
                seconds(fixture_ready - preparation_start), seconds(plans_ready - fixture_ready),
                seconds(preparation_end - plans_ready), seconds(preparation_end - preparation_start));
    std::puts("timing_note\tapply timing excludes input changes/checksum; two steady_clock reads per call, timer overhead not subtracted; loop timing includes input changes/checksum; no warmup apply calls");
    std::puts("repeat\tapply_calls\tapply_seconds\tapply_ns_per_call\tloop_seconds\tloop_ns_per_call\tchecksum\texpected_checksum");
    for (size_t repeat = 0; repeat < kRepeats; ++repeat) {
        const auto &r = results[repeat];
        std::printf("%zu\t%zu\t%.9f\t%.3f\t%.9f\t%.3f\t%.17g\t%.17g\n", repeat + 1,
                    kCalls, r.apply_seconds, 1e9 * r.apply_seconds / kCalls,
                    r.total_seconds, 1e9 * r.total_seconds / kCalls, r.checksum, r.expected_checksum);
    }
    return 0;
}
