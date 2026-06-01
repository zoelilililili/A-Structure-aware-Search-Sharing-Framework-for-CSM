#include "query_bctree.h"
#include "graph/graph.h"
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
namespace {

// n1/n2 can be adjusted directly here.
constexpr uint32_t kSubstructureSizeN1 = 3;

static void sort_unique_u32(std::vector<uint32_t>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

static std::vector<std::vector<uint32_t>> build_adjacency(const Graph* query_graph) {
    const uint32_t n = query_graph->getVerticesCount();
    std::vector<std::vector<uint32_t>> adj(n);
    for (uint32_t u = 0; u < n; ++u) {
        ui cnt = 0;
        const ui* nbr = query_graph->getVertexNeighbors(static_cast<VertexID>(u), cnt);
        adj[u].reserve(cnt);
        for (ui i = 0; i < cnt; ++i) {
            const uint32_t v = nbr[i];
            if (v < n) {
                adj[u].push_back(v);
            }
        }
        sort_unique_u32(adj[u]);
    }
    return adj;
}

static bool is_connected_subset(const std::vector<uint32_t>& subset,
                                const std::vector<std::vector<uint32_t>>& adj,
                                const std::vector<char>& in_subset) {
    if (subset.empty()) {
        return false;
    }

    std::queue<uint32_t> q;
    std::vector<char> visited(in_subset.size(), 0);

    q.push(subset[0]);
    visited[subset[0]] = 1;

    uint32_t reached = 0;
    while (!q.empty()) {
        const uint32_t u = q.front();
        q.pop();
        ++reached;

        for (uint32_t v : adj[u]) {
            if (!in_subset[v] || visited[v]) {
                continue;
            }
            visited[v] = 1;
            q.push(v);
        }
    }

    return reached == subset.size();
}

static bool has_exactly_one_boundary_vertex(const std::vector<uint32_t>& subset,
                                            const std::vector<std::vector<uint32_t>>& adj,
                                            const std::vector<char>& in_subset,
                                            uint32_t& boundary_vertex) {
    boundary_vertex = std::numeric_limits<uint32_t>::max();
    uint32_t count = 0;

    for (uint32_t u : subset) {
        bool touches_outside = false;
        for (uint32_t v : adj[u]) {
            if (!in_subset[v]) {
                touches_outside = true;
                break;
            }
        }

        if (touches_outside) {
            ++count;
            boundary_vertex = u;
            if (count > 1) {
                return false;
            }
        }
    }

    return count == 1;
}

static std::vector<bool> compute_articulation_points(const std::vector<std::vector<uint32_t>>& adj) {
    const uint32_t n = static_cast<uint32_t>(adj.size());
    std::vector<bool> is_cut(n, false);
    std::vector<int> tin(n, -1);
    std::vector<int> low(n, -1);
    int timer = 0;

    std::function<void(uint32_t, int)> dfs = [&](uint32_t u, int parent) {
        tin[u] = low[u] = timer++;
        int children = 0;

        for (uint32_t v : adj[u]) {
            if (static_cast<int>(v) == parent) {
                continue;
            }

            if (tin[v] != -1) {
                low[u] = std::min(low[u], tin[v]);
            } else {
                dfs(v, static_cast<int>(u));
                low[u] = std::min(low[u], low[v]);
                if (parent != -1 && low[v] >= tin[u]) {
                    is_cut[u] = true;
                }
                ++children;
            }
        }

        if (parent == -1 && children > 1) {
            is_cut[u] = true;
        }
    };

    for (uint32_t u = 0; u < n; ++u) {
        if (tin[u] == -1) {
            dfs(u, -1);
        }
    }

    return is_cut;
}

struct CandidateBlock {
    std::vector<uint32_t> vertices;
    uint32_t boundary_vertex = std::numeric_limits<uint32_t>::max();
};

static bool has_label_overlap_inside_outside(
    const std::vector<uint32_t>& subset,
    const std::vector<LabelID>& vertex_labels,
    const std::unordered_map<LabelID, uint32_t>& total_label_count) {
    std::unordered_map<LabelID, uint32_t> inside_label_count;
    inside_label_count.reserve(subset.size() * 2 + 1);

    for (uint32_t u : subset) {
        if (u < vertex_labels.size()) {
            inside_label_count[vertex_labels[u]] += 1;
        }
    }

    for (const auto& kv : inside_label_count) {
        const LabelID label = kv.first;
        const uint32_t inside_cnt = kv.second;
        auto it = total_label_count.find(label);
        if (it != total_label_count.end() && it->second > inside_cnt) {
            // Same label appears outside the subset as well.
            return true;
        }
    }
    return false;
}

static std::vector<CandidateBlock> enumerate_candidate_blocks(const std::vector<std::vector<uint32_t>>& adj,
                                                              const std::vector<LabelID>& vertex_labels,
                                                              const std::unordered_map<LabelID, uint32_t>& total_label_count,
                                                              uint32_t n1,
                                                              uint32_t n2) {
    const uint32_t n = static_cast<uint32_t>(adj.size());
    std::vector<CandidateBlock> candidates;

    if (n == 0) {
        return candidates;
    }

    const uint32_t low_k = std::min(n1, n2);
    const uint32_t high_k = std::min<uint32_t>(n, std::max(n1, n2));
    if (low_k == 0 || low_k > high_k) {
        return candidates;
    }

    std::vector<char> in_subset(n, 0);
    std::vector<uint32_t> subset;

    for (uint32_t k = low_k; k <= high_k; ++k) {
        subset.clear();
        subset.reserve(k);
        std::fill(in_subset.begin(), in_subset.end(), 0);

        std::function<void(uint32_t, uint32_t)> dfs_choose = [&](uint32_t start, uint32_t need) {
            if (need == 0) {
                if (!is_connected_subset(subset, adj, in_subset)) {
                    return;
                }

                if (has_label_overlap_inside_outside(subset, vertex_labels, total_label_count)) {
                    return;
                }

                uint32_t boundary_vertex = std::numeric_limits<uint32_t>::max();
                if (!has_exactly_one_boundary_vertex(subset, adj, in_subset, boundary_vertex)) {
                    return;
                }

                candidates.push_back({subset, boundary_vertex});
                return;
            }

            if (start >= n || (n - start) < need) {
                return;
            }

            for (uint32_t u = start; u < n; ++u) {
                subset.push_back(u);
                in_subset[u] = 1;
                dfs_choose(u + 1, need - 1);
                in_subset[u] = 0;
                subset.pop_back();
            }
        };

        dfs_choose(0, k);
    }

    return candidates;
}

static std::vector<CandidateBlock> select_disjoint_blocks(std::vector<CandidateBlock> candidates,
                                                           uint32_t vertex_count) {
    std::sort(candidates.begin(), candidates.end(), [](const CandidateBlock& a, const CandidateBlock& b) {
        if (a.vertices.size() != b.vertices.size()) {
            return a.vertices.size() > b.vertices.size();
        }
        if (a.boundary_vertex != b.boundary_vertex) {
            return a.boundary_vertex < b.boundary_vertex;
        }
        return a.vertices < b.vertices;
    });

    std::vector<char> used(vertex_count, 0);
    std::vector<CandidateBlock> selected;

    for (const CandidateBlock& c : candidates) {
        bool overlap = false;
        for (uint32_t u : c.vertices) {
            if (u < used.size() && used[u]) {
                overlap = true;
                break;
            }
        }

        if (overlap) {
            continue;
        }

        selected.push_back(c);
        for (uint32_t u : c.vertices) {
            if (u < used.size()) {
                used[u] = 1;
            }
        }
    }

    return selected;
}

static std::vector<std::vector<uint32_t>> build_block_neighbors(
    const std::vector<std::vector<uint32_t>>& adj,
    const std::vector<std::vector<uint32_t>>& blocks) {
    const uint32_t n = static_cast<uint32_t>(adj.size());
    const uint32_t B = static_cast<uint32_t>(blocks.size());

    std::vector<std::vector<uint32_t>> block_neighbors(B);
    std::vector<int> vertex_to_block(n, -1);

    for (uint32_t bi = 0; bi < B; ++bi) {
        for (uint32_t u : blocks[bi]) {
            if (u < n) {
                vertex_to_block[u] = static_cast<int>(bi);
            }
        }
    }

    for (uint32_t u = 0; u < n; ++u) {
        const int bu = vertex_to_block[u];
        if (bu < 0) {
            continue;
        }

        for (uint32_t v : adj[u]) {
            const int bv = vertex_to_block[v];
            if (bv < 0 || bv == bu) {
                continue;
            }
            block_neighbors[static_cast<size_t>(bu)].push_back(static_cast<uint32_t>(bv));
            block_neighbors[static_cast<size_t>(bv)].push_back(static_cast<uint32_t>(bu));
        }
    }

    for (auto& nbrs : block_neighbors) {
        sort_unique_u32(nbrs);
    }

    return block_neighbors;
}

} // namespace

enum class GraphTextFormat {
    kUnknown,
    kWithMeta,
    kWithoutMeta
};

static GraphTextFormat detect_graph_text_format(const std::string& graph_path) {
    std::ifstream in(graph_path);
    if (!in.is_open()) {
        return GraphTextFormat::kUnknown;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        std::istringstream iss(line);
        char type = 0;
        iss >> type;
        if (type == 't') {
            return GraphTextFormat::kWithMeta;
        }
        if (type == 'v' || type == 'e') {
            return GraphTextFormat::kWithoutMeta;
        }
    }

    return GraphTextFormat::kUnknown;
}

static bool detect_edge_labeled_text_graph(const std::string& graph_path) {
    std::ifstream in(graph_path);
    if (!in.is_open()) {
        return false;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        std::istringstream iss(line);
        char type = 0;
        iss >> type;
        if (type != 'e') {
            continue;
        }

        uint32_t v = 0;
        int value_count = 0;
        while (iss >> v) {
            ++value_count;
        }
        // e u v label -> 3 numbers; e u v -> 2 numbers.
        return value_count >= 3;
    }

    return false;
}

void QueryBCTreeDecomposition::clear() {
    query_vertex_count = 0;
    block_count = 0;
    blocks.clear();
    block_neighbor_block_indices.clear();
    is_cut_vertex.clear();
    vertex_block_membership_count.clear();
    blocks.clear();
    is_inner_non_cut.clear();
}

QueryBCTreeDecomposition build_query_bctree_decomposition(const Graph* query_graph) {
    QueryBCTreeDecomposition out;
    if (!query_graph || query_graph->getVerticesCount() == 0) {
        return out;
    }

    out.query_vertex_count = query_graph->getVerticesCount();
    const auto adj = build_adjacency(query_graph);

    std::vector<LabelID> vertex_labels(out.query_vertex_count, 0);
    std::unordered_map<LabelID, uint32_t> total_label_count;
    total_label_count.reserve(out.query_vertex_count * 2 + 1);
    for (uint32_t u = 0; u < out.query_vertex_count; ++u) {
        LabelID label = query_graph->getVertexLabel(static_cast<VertexID>(u));
        vertex_labels[u] = label;
        total_label_count[label] += 1;
    }

    // 判断图节点数量，多则跑三分之一，少则跑二分之一
    if (out.query_vertex_count >= 10) {
        out.block_count = out.query_vertex_count / 3;
    } else {
        out.block_count = out.query_vertex_count / 2;
    }
    const uint32_t dynamic_n2 = out.block_count;

    auto candidates = enumerate_candidate_blocks(adj, vertex_labels, total_label_count,
                                                 kSubstructureSizeN1, dynamic_n2);
    auto selected = select_disjoint_blocks(std::move(candidates), out.query_vertex_count);

    out.blocks.clear();
    out.blocks.reserve(selected.size());
    for (auto& b : selected) {
        sort_unique_u32(b.vertices);
        out.blocks.push_back(std::move(b.vertices));
    }

    out.block_count = static_cast<uint32_t>(out.blocks.size());
    out.block_neighbor_block_indices = build_block_neighbors(adj, out.blocks);

    out.is_cut_vertex.assign(out.query_vertex_count, false);
    out.vertex_block_membership_count.assign(out.query_vertex_count, 0);
    out.blocks.assign(out.query_vertex_count, {});
    out.is_inner_non_cut.assign(out.query_vertex_count, 0);

    for (uint32_t bi = 0; bi < out.block_count; ++bi) {
        for (uint32_t u : out.blocks[bi]) {
            if (u < out.vertex_block_membership_count.size()) {
                out.vertex_block_membership_count[u] += 1;
            }
            if (u < out.blocks.size()) {
                out.blocks[u].push_back(bi);
            }
        }
    }

    {
        std::vector<char> in_block(out.query_vertex_count, 0);
        for (uint32_t bi = 0; bi < out.block_count; ++bi) {
            std::fill(in_block.begin(), in_block.end(), 0);
            for (uint32_t u : out.blocks[bi]) {
                if (u < in_block.size()) in_block[u] = 1;
            }
            for (uint32_t u : out.blocks[bi]) {
                bool touches_outside = false;
                for (uint32_t v : adj[u]) {
                    if (v < in_block.size() && !in_block[v]) {
                        touches_outside = true;
                        break;
                    }
                }
                if (touches_outside) {
                    out.is_cut_vertex[u] = true;
                }
            }
        }
    }

    for (uint32_t u = 0; u < out.query_vertex_count; ++u) {
        out.is_inner_non_cut[u] = out.is_cut_vertex[u] ? 0 : 1;
    }

    return out;
}

void print_query_bctree_decomposition(const QueryBCTreeDecomposition& d, std::ostream& os) {
    os << "[Phase3][Block Decomposition] vertices=" << d.query_vertex_count
       << " blocks=" << d.block_count << '\n';
    for (uint32_t bi = 0; bi < d.block_count; ++bi) {
        os << "  Block " << bi << " vertices:";
        for (uint32_t u : d.blocks[bi]) {
            os << ' ' << u;
        }
        os << " | neighbor blocks:";
        for (uint32_t nb : d.block_neighbor_block_indices[bi]) {
            os << ' ' << nb;
        }
        os << '\n';
    }
    os << "  Block boundary vertices (block reuse entry):";
    for (uint32_t u = 0; u < d.query_vertex_count; ++u) {
        if (d.is_cut_vertex[u]) {
            os << ' ' << u;
        }
    }
    os << '\n';
    os << "  Vertex -> #blocks containing it:";
    for (uint32_t u = 0; u < d.query_vertex_count; ++u) {
        if (d.vertex_block_membership_count[u] != 0) {
            os << " (" << u << ":" << d.vertex_block_membership_count[u] << ")";
        }
    }
    os << '\n';
}


