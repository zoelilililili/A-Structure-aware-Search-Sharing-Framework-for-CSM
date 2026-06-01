#include "subgraph_automorphism_lib.h"

#include <algorithm>
#include <cctype>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "graph/graph.h"

namespace subgraph_automorphism_lib {
namespace {

struct Subgraph {
    uint32_t k = 0;
    std::vector<uint32_t> labels;        
    std::vector<uint32_t> degree;        
    std::vector<uint8_t> adj;            
    std::vector<uint32_t> edge_label;    
};

static Subgraph build_induced_subgraph(const Graph& g,
                                       const std::vector<uint32_t>& vertices_sorted,
                                       const std::vector<std::vector<uint8_t>>& full_adj,
                                       const std::vector<std::vector<uint32_t>>& full_elabel,
                                       const std::vector<uint32_t>& full_vlabel) {
    Subgraph sg;
    sg.k = static_cast<uint32_t>(vertices_sorted.size());
    uint32_t k = sg.k;
    sg.labels.resize(k);
    sg.degree.assign(k, 0);
    sg.adj.assign((size_t)k * k, 0);
    sg.edge_label.assign((size_t)k * k, std::numeric_limits<uint32_t>::max());

    std::unordered_map<uint32_t, uint32_t> old_to_new;
    old_to_new.reserve(k * 2);
    for (uint32_t i = 0; i < k; ++i) {
        old_to_new[vertices_sorted[i]] = i;
        sg.labels[i] = full_vlabel[vertices_sorted[i]];
    }

    for (uint32_t i = 0; i < k; ++i) {
        uint32_t u_old = vertices_sorted[i];
        for (uint32_t j = 0; j < k; ++j) {
            if (i == j) continue;
            uint32_t v_old = vertices_sorted[j];
            if (full_adj[u_old][v_old]) {
                sg.adj[(size_t)i * k + j] = 1;
                sg.edge_label[(size_t)i * k + j] = full_elabel[u_old][v_old];
                sg.degree[i] += 1;
            }
        }
    }
    (void)g;
    return sg;
}

static bool is_isomorphism_preserving(const Subgraph& sg,
                                      const std::vector<int32_t>& map_u_to_v,
                                      uint32_t u, uint32_t v) {
    if (sg.labels[u] != sg.labels[v]) return false;
    if (sg.degree[u] != sg.degree[v]) return false;

    uint32_t k = sg.k;
    for (uint32_t uu = 0; uu < k; ++uu) {
        int32_t vv = map_u_to_v[uu];
        if (vv < 0) continue;
        uint32_t vvu = static_cast<uint32_t>(vv);

        uint8_t a1 = sg.adj[(size_t)u * k + uu];
        uint8_t a2 = sg.adj[(size_t)v * k + vvu];
        if (a1 != a2) return false;
        if (a1) {
            uint32_t l1 = sg.edge_label[(size_t)u * k + uu];
            uint32_t l2 = sg.edge_label[(size_t)v * k + vvu];
            if (l1 != l2) return false;
        }

        uint8_t b1 = sg.adj[(size_t)uu * k + u];
        uint8_t b2 = sg.adj[(size_t)vvu * k + v];
        if (b1 != b2) return false;
        if (b1) {
            uint32_t l1r = sg.edge_label[(size_t)uu * k + u];
            uint32_t l2r = sg.edge_label[(size_t)vvu * k + v];
            if (l1r != l2r) return false;
        }
    }
    return true;
}

static void build_candidates(const Subgraph& sg,
                             std::vector<std::vector<uint32_t>>& cand) {
    uint32_t k = sg.k;
    cand.assign(k, {});
    for (uint32_t u = 0; u < k; ++u) {
        for (uint32_t v = 0; v < k; ++v) {
            if (sg.labels[u] == sg.labels[v] && sg.degree[u] == sg.degree[v]) {
                cand[u].push_back(v);
            }
        }
    }
}

static void order_vertices(const std::vector<std::vector<uint32_t>>& cand,
                           std::vector<uint32_t>& order) {
    uint32_t k = static_cast<uint32_t>(cand.size());
    order.resize(k);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
        return cand[a].size() < cand[b].size();
    });
}

static void enumerate_automorphisms_dfs(const Subgraph& sg,
                                        const std::vector<std::vector<uint32_t>>& cand,
                                        const std::vector<uint32_t>& order,
                                        uint32_t depth,
                                        std::vector<int32_t>& map_u_to_v,
                                        std::vector<int32_t>& map_v_to_u,
                                        std::vector<std::vector<uint32_t>>& autos,
                                        uint64_t max_autos) {
    if (autos.size() >= max_autos) return;
    uint32_t k = sg.k;
    if (depth == k) {
        std::vector<uint32_t> perm(k);
        for (uint32_t u = 0; u < k; ++u) perm[u] = static_cast<uint32_t>(map_u_to_v[u]);
        autos.push_back(std::move(perm));
        return;
    }

    uint32_t u = order[depth];
    for (uint32_t v : cand[u]) {
        if (map_v_to_u[v] >= 0) continue;
        if (!is_isomorphism_preserving(sg, map_u_to_v, u, v)) continue;

        map_u_to_v[u] = static_cast<int32_t>(v);
        map_v_to_u[v] = static_cast<int32_t>(u);

        enumerate_automorphisms_dfs(sg, cand, order, depth + 1, map_u_to_v, map_v_to_u, autos, max_autos);

        map_u_to_v[u] = -1;
        map_v_to_u[v] = -1;
        if (autos.size() >= max_autos) return;
    }
}

static bool has_nontrivial_automorphism(const Subgraph& sg, uint64_t max_autos_to_search = 2) {
    std::vector<std::vector<uint32_t>> cand;
    build_candidates(sg, cand);
    std::vector<uint32_t> order;
    order_vertices(cand, order);

    std::vector<int32_t> map_u_to_v(sg.k, -1);
    std::vector<int32_t> map_v_to_u(sg.k, -1);
    std::vector<std::vector<uint32_t>> autos;
    autos.reserve((size_t)max_autos_to_search);
    enumerate_automorphisms_dfs(sg, cand, order, 0, map_u_to_v, map_v_to_u, autos, max_autos_to_search);

    for (auto& p : autos) {
        bool identity = true;
        for (uint32_t i = 0; i < sg.k; ++i) {
            if (p[i] != i) { identity = false; break; }
        }
        if (!identity) return true;
    }
    return false;
}

static void enumerate_all_automorphisms(const Subgraph& sg,
                                       std::vector<std::vector<uint32_t>>& autos,
                                       uint64_t max_autos = 500000) {
    std::vector<std::vector<uint32_t>> cand;
    build_candidates(sg, cand);
    std::vector<uint32_t> order;
    order_vertices(cand, order);

    std::vector<int32_t> map_u_to_v(sg.k, -1);
    std::vector<int32_t> map_v_to_u(sg.k, -1);

    autos.clear();
    autos.reserve(1024);
    enumerate_automorphisms_dfs(sg, cand, order, 0, map_u_to_v, map_v_to_u, autos, max_autos);
}

static bool induced_subgraph_connected(uint64_t mask, uint32_t n,
                                       const std::vector<std::vector<uint8_t>>& adj) {
    const int cnt = __builtin_popcountll(mask);
    if (cnt <= 1) return true;

    uint32_t start = n;
    for (uint32_t i = 0; i < n; ++i) {
        if (mask & (1ULL << i)) { start = i; break; }
    }
    if (start >= n) return false;

    std::vector<uint8_t> vis(n, 0);
    std::vector<uint32_t> stack;
    stack.push_back(start);
    vis[start] = 1;
    uint32_t seen = 0;
    while (!stack.empty()) {
        uint32_t u = stack.back();
        stack.pop_back();
        ++seen;
        for (uint32_t v = 0; v < n; ++v) {
            if (!(mask & (1ULL << v))) continue;
            if (!adj[u][v]) continue;
            if (!vis[v]) {
                vis[v] = 1;
                stack.push_back(v);
            }
        }
    }
    return seen == static_cast<uint32_t>(cnt);
}

static bool is_subset(uint64_t a, uint64_t b) {
    return (a & b) == a;
}

static std::vector<uint32_t> vertices_from_mask(uint64_t mask, uint32_t n) {
    std::vector<uint32_t> v;
    for (uint32_t i = 0; i < n; ++i) if (mask & (1ULL << i)) v.push_back(i);
    return v;
}

static bool is_identity_perm(const std::vector<uint32_t>& p) {
    for (size_t i = 0; i < p.size(); ++i) if (p[i] != i) return false;
    return true;
}

static void sort_automorphisms_identity_first(std::vector<std::vector<uint32_t>>& autos) {
    std::sort(autos.begin(), autos.end(), [](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
        bool ia = is_identity_perm(a);
        bool ib = is_identity_perm(b);
        if (ia != ib) return ia;
        return a < b;
    });
}

static uint64_t moved_old_vertices_mask_from_automorphisms(const std::vector<uint32_t>& old_vertices_sorted,
                                                          const std::vector<std::vector<uint32_t>>& autos) {
    uint64_t mask = 0;
    if (old_vertices_sorted.size() > 63) return 0;
    for (const auto& perm : autos) {
        if (is_identity_perm(perm)) continue;
        const size_t k = perm.size();
        for (size_t i = 0; i < k; ++i) {
            if (perm[i] == i) continue;
            uint32_t old_id = old_vertices_sorted[i];
            if (old_id < 63) mask |= (1ULL << old_id);
        }
    }
    return mask;
}

}

void build_maximal_automorphic_induced_subgraphs_in_memory(const Graph& query_graph,
                                                          uint64_t max_autos,
                                                          std::vector<SubAutomorphismResult>& out) {
    out.clear();

    const uint32_t n = query_graph.getVerticesCount();
    if (n == 0) return;
    if (n > 63) {
        return;
    }

    std::vector<uint32_t> full_vlabel(n);
    for (uint32_t u = 0; u < n; ++u) full_vlabel[u] = query_graph.getVertexLabel(u);

    std::vector<std::vector<uint8_t>> full_adj(n, std::vector<uint8_t>(n, 0));
    std::vector<std::vector<uint32_t>> full_elabel(n, std::vector<uint32_t>(n, std::numeric_limits<uint32_t>::max()));
    for (uint32_t u = 0; u < n; ++u) {
        uint32_t cnt = 0;
        const ui* nbrs = query_graph.getVertexNeighbors(u, cnt);
        for (uint32_t i = 0; i < cnt; ++i) {
            uint32_t v = nbrs[i];
            full_adj[u][v] = 1;
            full_elabel[u][v] = query_graph.getEdgeLabelByLocalOffset(u, i);
        }
    }

    std::vector<uint64_t> automorphic_masks;
    automorphic_masks.reserve(1024);

    const uint32_t min_size = 3;
    const uint32_t max_size = (n >= 10) ? (n / 3) : (n / 2);

    for (uint64_t mask = 1; mask < (1ULL << n) - 1ULL; ++mask) {
        uint32_t pop = (uint32_t)__builtin_popcountll(mask);
        if (pop < min_size || pop > max_size) continue;
        if (!induced_subgraph_connected(mask, n, full_adj)) continue;

        auto verts = vertices_from_mask(mask, n);
        Subgraph sg = build_induced_subgraph(query_graph, verts, full_adj, full_elabel, full_vlabel);
        if (has_nontrivial_automorphism(sg)) {
            automorphic_masks.push_back(mask);
        }
    }

    std::sort(automorphic_masks.begin(), automorphic_masks.end(), [](uint64_t a, uint64_t b) {
        int pa = __builtin_popcountll(a);
        int pb = __builtin_popcountll(b);
        if (pa != pb) return pa > pb;
        return a < b;
    });

    std::vector<uint64_t> maximal;
    maximal.reserve(automorphic_masks.size());
    for (uint64_t m : automorphic_masks) {
        bool contained = false;
        for (uint64_t mm : maximal) {
            if (is_subset(m, mm)) { contained = true; break; }
        }
        if (!contained) maximal.push_back(m);
    }

    out.reserve(maximal.size());
    std::unordered_map<uint64_t, size_t> support_to_out_idx;
    support_to_out_idx.reserve(maximal.size() * 2);
    for (size_t idx = 0; idx < maximal.size(); ++idx) {
        uint64_t mask = maximal[idx];
        auto verts = vertices_from_mask(mask, n);
        Subgraph sg = build_induced_subgraph(query_graph, verts, full_adj, full_elabel, full_vlabel);

        std::vector<std::vector<uint32_t>> autos;
        enumerate_all_automorphisms(sg, autos, max_autos);
        if (autos.empty()) continue;
        sort_automorphisms_identity_first(autos);

        SubAutomorphismResult r;
        r.old_vertices_sorted = std::move(verts);
        r.moved_old_vertices_mask = moved_old_vertices_mask_from_automorphisms(r.old_vertices_sorted, autos);
        r.automorphisms = std::move(autos);

        r.vertex_list.reserve(sg.k);
        for (uint32_t u = 0; u < sg.k; ++u) {
            r.vertex_list.emplace_back(u, sg.labels[u]);
        }
        for (uint32_t u = 0; u < sg.k; ++u) {
            for (uint32_t v = u + 1; v < sg.k; ++v) {
                if (sg.adj[(size_t)u * sg.k + v]) {
                    uint32_t el = sg.edge_label[(size_t)u * sg.k + v];
                    r.edge_list.emplace_back(u, v, el);
                }
            }
        }

        auto it = support_to_out_idx.find(r.moved_old_vertices_mask);
        if (it == support_to_out_idx.end()) {
            r.sub_id = static_cast<uint32_t>(out.size());
            support_to_out_idx.emplace(r.moved_old_vertices_mask, out.size());
            out.push_back(std::move(r));
        } else {
            size_t existing_idx = it->second;
            // Prefer the result with fewer vertices (less trivial "padding").
            if (out[existing_idx].old_vertices_sorted.size() > r.old_vertices_sorted.size()) {
                r.sub_id = out[existing_idx].sub_id;
                out[existing_idx] = std::move(r);
            }
        }
    }
}

} 

