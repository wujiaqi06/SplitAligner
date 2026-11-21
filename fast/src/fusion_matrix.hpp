// Private matrix-assembly layer. Projection and gene-edge IDs refer to one input row.
#pragma once
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>
#include <vector>

namespace fusion_matrix {

using Members = std::vector<size_t>;

struct Component {
    Members members;
    std::string split;
    size_t edge;
    Members source_nodes;
    double length;
};

struct Row {
    std::string gene;
    std::vector<Component> components;
    std::vector<size_t> axis_to_component;
};

// Store the per-gene observations, since later genes may introduce new columns.
// This private spool never reopens the user's input (which may be a pipe).
// POSIX unlinks it immediately; Windows removes it after closing, before publication.
class Assembly {
    std::fstream spool;
    std::filesystem::path pending_spool;
    std::set<Members> targets;
    size_t row_count = 0;
    size_t axis_count;

    template<class T> void put(T value) {
        spool.write(reinterpret_cast<const char*>(&value), sizeof(value));
    }
    template<class T> T get() {
        T value;
        spool.read(reinterpret_cast<char*>(&value), sizeof(value));
        return value;
    }
    template<class Container> void put_array(const Container& values) {
        put(values.size());
        if (!values.empty())
            spool.write(reinterpret_cast<const char*>(values.data()),
                        sizeof(typename Container::value_type) * values.size());
    }
    template<class Container> Container get_array() {
        Container values;
        values.resize(get<size_t>());
        if (!values.empty())
            spool.read(reinterpret_cast<char*>(values.data()),
                       sizeof(typename Container::value_type) * values.size());
        return values;
    }

    Row read_row() {
        Row row;
        row.gene = get_array<std::string>();
        row.axis_to_component.assign(axis_count, SAF_NONE);
        const auto count = get<size_t>();
        row.components.reserve(count);
        for (size_t p = 0; p < count; ++p) {
            auto members = get_array<Members>();
            for (auto axis : members) row.axis_to_component.at(axis) = p;
            row.components.push_back({std::move(members), get_array<std::string>(),
                                      get<size_t>(), get_array<Members>(), get<double>()});
        }
        return row;
    }

    static std::string name(const Members& target) {
        std::string result;
        for (auto axis : target) {
            if (!result.empty()) result += '|';
            result += 'A' + std::to_string(axis);
        }
        return result;
    }

    static void write_component(std::ostream& out, size_t id, const Component& component) {
        out << "{\"projection_id\":" << id << ",\"members\":";
        array(out, component.members);
        out << ",\"split_hex\":" << json(component.split) << ",\"gene_edge_id\":";
        if (component.edge == SAF_NONE) out << "null"; else out << component.edge;
        out << ",\"source_nodes\":";
        array(out, component.source_nodes);
        out << ",\"value\":";
        if (std::isfinite(component.length)) out << component.length; else out << "null";
        out << '}';
    }

    static void write_cell(std::ostream& matrix, std::ostream& audit, size_t gene_idx,
                           const Row& row, const Members& target, const std::string& label,
                           bool direct_only) {
        Members projections;
        bool deleted = false;
        for (auto axis : target) {
            const auto p = row.axis_to_component.at(axis);
            if (p == SAF_NONE) deleted = true;
            else projections.push_back(p);
        }
        std::sort(projections.begin(), projections.end());
        projections.erase(std::unique(projections.begin(), projections.end()), projections.end());
        bool partial = false, recovered = true, complete_lengths = true;
        for (auto p : projections) {
            const auto& component = row.components[p];
            if (!std::includes(target.begin(), target.end(), component.members.begin(),
                               component.members.end())) partial = true;
            if (component.edge == SAF_NONE) recovered = false;
            if (!std::isfinite(component.length)) complete_lengths = false;
        }

        const bool representable = !deleted && !partial;
        const bool active = representable && projections.size() == 1;
        const char* structure = representable ? (active ? "active" : "decomposable") : "unavailable";
        const char* recovery = !representable ? "not_applicable" : recovered ? "recovered" : "not_recovered";
        const char* numeric = "available";
        const char* mode = !representable ? "none" : active ? "single_edge_lookup" : "selected_component_edge_sum";
        double value = NAN;
        if (!representable) numeric = "not_applicable";
        else if (!recovered) numeric = "not_recovered";
        // Suppression is a view decision, not an attempted numeric sum.
        // Retain each component's length and source in the audit below.
        else if (!active && direct_only) numeric = "direct_only";
        else if (!complete_lengths) numeric = "missing_length";
        else {
            // Distinct projection IDs in an exact cover have distinct matched gene edges.
            // Apply already aggregated a binary display root into one unrooted edge.
            if (active) {
                value = row.components[projections.front()].length;
            } else {
                std::vector<double> terms;
                for (auto p : projections) terms.push_back(row.components[p].length);
                // Give all front ends the same binary64 addition order, independent
                // of the reference axis numbering and empirical node order.
                std::sort(terms.begin(), terms.end(), [](double a, double b) {
                    if (std::abs(a) == std::abs(b)) return a < b;
                    return std::abs(a) < std::abs(b);
                });
                value = 0;
                for (auto term : terms) value += term;
            }
            if (!std::isfinite(value)) numeric = "nonfinite_sum";
        }
        matrix << '\t';
        if (std::isfinite(value)) matrix << value; else matrix << "NA";
        audit << "{\"gene_idx\":" << gene_idx << ",\"gene_id\":" << json(row.gene)
              << ",\"column\":" << json(label) << ",\"members\":";
        array(audit, target);
        audit << ",\"reference_state\":\"" << structure << "\",\"deleted_member\":"
              << (deleted ? "true" : "false") << ",\"partial_component\":" << (partial ? "true" : "false")
              << ",\"recovery_status\":\"" << recovery << "\",\"numeric_status\":\"" << numeric
              << "\",\"construction_mode\":\"" << mode << "\",\"value\":";
        if (std::isfinite(value)) audit << value; else audit << "null";
        audit << ",\"components\":[";
        for (size_t i = 0; i < projections.size(); ++i) {
            if (i) audit << ',';
            const auto p = projections[i];
            write_component(audit, p, row.components[p]);
        }
        audit << "]}\n";
    }

public:
    Assembly(const std::filesystem::path& stage, size_t axes) : axis_count(axes) {
        const auto path = stage / ".fusion-spool";
        spool.exceptions(std::ios::badbit | std::ios::failbit);
        spool.open(path, std::ios::binary | std::ios::out | std::ios::in | std::ios::trunc);
#if defined(_WIN32) || defined(SAF_TEST_NAMED_SPOOL)
        pending_spool = path;
#else
        std::filesystem::remove(path);
#endif
    }

    ~Assembly() noexcept {
        // Error paths retain the failed staging directory, but try to remove
        // the private spool after releasing the Windows file handle.
        try {
            spool.exceptions(std::ios::goodbit);
            if (spool.is_open()) spool.close();
            if (!pending_spool.empty()) {
                std::error_code error;
                std::filesystem::remove(pending_spool, error);
                if (error) std::cerr << "Fusion spool cleanup failed: " << pending_spool.u8string()
                                     << ": " << error.message() << '\n';
            }
        } catch (...) { /* Never replace the original failure while unwinding. */ }
    }

    template<class Sources>
    void record(const std::string& gene, const Coverage& coverage,
                const std::vector<size_t>& matched_edges, const Values& values, Sources sources) {
        put_array(gene);
        put(coverage.groups.size());
        for (size_t p = 0; p < coverage.groups.size(); ++p) {
            const auto& group = coverage.groups[p];
            if (group.members.size() > 1) targets.insert(group.members);
            put_array(group.members);
            put_array(group.hex);
            const auto edge = matched_edges[p];
            put(edge);
            if (edge == SAF_NONE) put_array(Members{}); else put_array(sources(edge));
            put(values.a[p] ? values.v[p] : double(NAN));
        }
        ++row_count;
    }

    void finish(std::ostream& matrix, std::ostream& audit, bool direct_only) {
        std::vector<std::string> labels;
        matrix << "gene_id";
        for (const auto& target : targets) {
            labels.push_back(name(target));
            matrix << '\t' << labels.back();
        }
        matrix << '\n';
        spool.flush();
        spool.seekg(0);
        for (size_t i = 0; i < row_count; ++i) {
            auto row = read_row();
            matrix << row.gene;
            size_t j = 0;
            for (const auto& target : targets)
                write_cell(matrix, audit, i, row, target, labels[j++], direct_only);
            matrix << '\n';
        }
        spool.close();
        if (!pending_spool.empty()) {
            if (!std::filesystem::remove(pending_spool))
                throw std::runtime_error("fusion spool disappeared before cleanup");
            pending_spool.clear();
        }
    }
};

} // namespace fusion_matrix
