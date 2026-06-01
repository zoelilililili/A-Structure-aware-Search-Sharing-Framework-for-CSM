#include "subgraph_isomorphism_lib.h"

#include <algorithm>
#include <cctype>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <map>
#include <set>

#include "graph/graph.h"

namespace subgraph_isomorphism_lib {
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

static bool induced_subgraph_connected(uint64_t mask, uint32_t n,
                                       const std::vector<std::vector<uint8_t>>& adj) {
    int cnt = __builtin_popcountll(mask);
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

static std::vector<uint32_t> vertices_from_mask(uint64_t mask, uint32_t n) {
    std::vector<uint32_t> v;
    for (uint32_t i = 0; i < n; ++i)
        if (mask & (1ULL << i)) v.push_back(i);
    return v;
}

static uint64_t mask_from_mapping(const std::vector<uint32_t>& mapping) {
    uint64_t mask = 0;
    for (uint32_t v : mapping) {
        mask |= (1ULL << v);
    }
    return mask;
}

static bool is_subset(uint64_t a, uint64_t b) {
    return (a & b) == a;
}

struct CanonicalEncoding {
    std::vector<uint32_t> vertex_labels_seq; 
    std::vector<uint8_t> adj_flat;            
    std::vector<uint32_t> edge_labels_seq;    

    bool operator==(const CanonicalEncoding& other) const {
        return vertex_labels_seq == other.vertex_labels_seq &&
               adj_flat == other.adj_flat &&
               edge_labels_seq == other.edge_labels_seq;
    }
};

struct CanonicalEncodingHash {
    size_t operator()(const CanonicalEncoding& enc) const {
        size_t h = 0;
        for (auto v : enc.vertex_labels_seq) h ^= std::hash<uint32_t>()(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
        for (auto a : enc.adj_flat) h ^= std::hash<uint8_t>()(a) + 0x9e3779b9 + (h << 6) + (h >> 2);
        for (auto e : enc.edge_labels_seq) h ^= std::hash<uint32_t>()(e) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

static void compute_canonical_form(const Subgraph& sg,
                                   CanonicalEncoding& out_enc,
                                   std::vector<uint32_t>& out_canon_to_local) {
    uint32_t k = sg.k;
    std::vector<uint32_t> perm(k);
    std::iota(perm.begin(), perm.end(), 0);

    CanonicalEncoding best_enc;
    bool first = true;
    std::vector<uint32_t> best_perm;

    do {
        CanonicalEncoding cur;
        cur.vertex_labels_seq.resize(k);
        for (uint32_t i = 0; i < k; ++i) {
            cur.vertex_labels_seq[i] = sg.labels[perm[i]];
        }

        cur.adj_flat.resize((size_t)k * k, 0);
        cur.edge_labels_seq.resize((size_t)k * k, std::numeric_limits<uint32_t>::max());
        for (uint32_t i = 0; i < k; ++i) {
            for (uint32_t j = 0; j < k; ++j) {
                uint32_t u = perm[i];
                uint32_t v = perm[j];
                uint8_t a = sg.adj[(size_t)u * k + v];
                cur.adj_flat[(size_t)i * k + j] = a;
                if (a) {
                    cur.edge_labels_seq[(size_t)i * k + j] = sg.edge_label[(size_t)u * k + v];
                }
            }
        }

        if (first) {
            best_enc = std::move(cur);
            best_perm = perm;
            first = false;
        } else {
            bool less = false;
            if (cur.vertex_labels_seq != best_enc.vertex_labels_seq) {
                less = cur.vertex_labels_seq < best_enc.vertex_labels_seq;
            } else if (cur.adj_flat != best_enc.adj_flat) {
                less = cur.adj_flat < best_enc.adj_flat;
            } else {
                less = cur.edge_labels_seq < best_enc.edge_labels_seq;
            }
            if (less) {
                best_enc = std::move(cur);
                best_perm = perm;
            }
        }
    } while (std::next_permutation(perm.begin(), perm.end()));

    out_enc = std::move(best_enc);
    out_canon_to_local = std::move(best_perm);
}


}

void build_maximal_isomorphic_induced_subgraphs_in_memory(const Graph& query_graph,
                                                          std::vector<SubIsomorphismGroupResult>& out,
                                                          const std::vector<uint64_t>* phase1_moved_masks) {
    out.clear();

    const uint32_t n = query_graph.getVerticesCount();
    if (n < 3) return;
    if (n > 63) {
        return;
    }

    std::vector<uint32_t> full_vlabel(n);
    for (uint32_t u = 0; u < n; ++u)
        full_vlabel[u] = query_graph.getVertexLabel(u);

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

    struct InstanceInfo {
        uint64_t mask;
        std::vector<uint32_t> old_vertices;       
        std::vector<uint32_t> canon_to_old;        
    };

    struct Group {
        CanonicalEncoding encoding;
        Subgraph canonical_subgraph;                  
        std::vector<InstanceInfo> instances;
    };
    std::unordered_map<CanonicalEncoding, Group, CanonicalEncodingHash> groups;

    uint64_t total_masks = (n == 64) ? UINT64_MAX : ((1ULL << n) - 1);
    const uint32_t max_size = n / 4;
    for (uint64_t mask = 1; mask <= total_masks; ++mask) {
        uint32_t pop = __builtin_popcountll(mask);
        if (pop < 3 || pop > max_size) continue;
        if (!induced_subgraph_connected(mask, n, full_adj)) continue;

        if (phase1_moved_masks) {
            bool conflict = false;
            for (uint64_t pm : *phase1_moved_masks) {
                if (mask & pm) { conflict = true; break; }
            }
            if (conflict) continue;
        }

        std::vector<uint32_t> verts = vertices_from_mask(mask, n);
        Subgraph sg = build_induced_subgraph(query_graph, verts, full_adj, full_elabel, full_vlabel);

        CanonicalEncoding enc;
        std::vector<uint32_t> canon_to_local; 
        compute_canonical_form(sg, enc, canon_to_local);

        std::vector<uint32_t> canon_to_old(sg.k);
        for (uint32_t i = 0; i < sg.k; ++i) {
            canon_to_old[i] = verts[canon_to_local[i]];
        }

        auto it = groups.find(enc);
        if (it == groups.end()) {
            Group new_group;
            new_group.encoding = enc;
            new_group.canonical_subgraph = std::move(sg);
            new_group.instances.push_back({mask, std::move(verts), std::move(canon_to_old)});
            groups.emplace(std::move(enc), std::move(new_group));
        } else {
            it->second.instances.push_back({mask, std::move(verts), std::move(canon_to_old)});
        }
    }

    if (groups.empty()) return;

    std::vector<Group*> sorted_groups;
    sorted_groups.reserve(groups.size());
    for (auto& kv : groups) sorted_groups.push_back(&kv.second);
    std::sort(sorted_groups.begin(), sorted_groups.end(),
              [](const Group* a, const Group* b) {
                  return a->canonical_subgraph.k > b->canonical_subgraph.k;
              });
    std::set<uint64_t> subsumed_masks;
    std::vector<Group*> maximal_groups;

    for (Group* g : sorted_groups) {
        bool has_maximal = false;
        for (auto& inst : g->instances) {
            if (subsumed_masks.count(inst.mask) == 0) {
                has_maximal = true;
                break;
            }
        }
        if (!has_maximal) continue;

        Group maximal_group;
        maximal_group.encoding = g->encoding;
        maximal_group.canonical_subgraph = g->canonical_subgraph;
        for (auto& inst : g->instances) {
            if (subsumed_masks.count(inst.mask) == 0) {
                maximal_group.instances.push_back(inst);
            }
        }
        if (maximal_group.instances.size() <= 1) continue;
        for (auto& inst : maximal_group.instances) {
            subsumed_masks.insert(inst.mask);
        }
        maximal_groups.push_back(new Group(std::move(maximal_group)));
    }

    std::vector<std::pair<uint64_t, Group*>> all_instance_masks;
    
    for (Group* g : maximal_groups) {
        for (auto& inst : g->instances) {
            uint64_t mask = 0;
            for (uint32_t v : inst.old_vertices) {
                mask |= (1ULL << v);
            }
            all_instance_masks.push_back({mask, g});
        }
    }

    std::set<uint64_t> truly_maximal_masks;
    
    for (size_t i = 0; i < all_instance_masks.size(); ++i) {
        uint64_t mask_i = all_instance_masks[i].first;
        bool is_maximal_i = true;
        
        for (size_t j = 0; j < all_instance_masks.size(); ++j) {
            if (i == j) continue;
            
            uint64_t mask_j = all_instance_masks[j].first;
            
            if (is_subset(mask_i, mask_j) && mask_i != mask_j) {
                is_maximal_i = false;
                break;
            }
        }
        
        if (is_maximal_i) {
            truly_maximal_masks.insert(mask_i);
        }
    }
    
    for (Group* g : maximal_groups) {
        delete g;
    }
    maximal_groups.clear();
    
    std::unordered_map<CanonicalEncoding, Group*, CanonicalEncodingHash> group_map;
    
    for (Group* g : sorted_groups) {
        Group* filtered_group = nullptr;
        
        for (auto& inst : g->instances) {
            uint64_t mask = 0;
            for (uint32_t v : inst.old_vertices) {
                mask |= (1ULL << v);
            }
            
            if (truly_maximal_masks.count(mask) > 0) {
                if (!filtered_group) {
                    filtered_group = new Group();
                    filtered_group->encoding = g->encoding;
                    filtered_group->canonical_subgraph = g->canonical_subgraph;
                }
                filtered_group->instances.push_back(inst);
            }
        }
        
        if (filtered_group && filtered_group->instances.size() > 1) {
            maximal_groups.push_back(filtered_group);
        }
    }
    for (Group* g : maximal_groups) {
        if (g->instances.size() <= 1) continue;
        std::vector<InstanceInfo> disjoint_instances;
        std::set<uint64_t> used_vertices;
        for (auto& inst : g->instances) {
            bool overlap = false;
            for (uint32_t v : inst.old_vertices) {
                if (used_vertices.count(1ULL << v)) { overlap = true; break; }
            }
            if (overlap) continue;
            for (uint32_t v : inst.old_vertices) {
                used_vertices.insert(1ULL << v);
            }
            disjoint_instances.push_back(std::move(inst));
        }
        g->instances = std::move(disjoint_instances);
    }

    out.reserve(maximal_groups.size());
    uint32_t group_id = 0;
    for (Group* g : maximal_groups) {
        SubIsomorphismGroupResult res;
        res.group_id = group_id++;

        uint32_t k = g->canonical_subgraph.k;
        const Subgraph& sg = g->canonical_subgraph;

        const CanonicalEncoding& enc = g->encoding;
        for (uint32_t i = 0; i < k; ++i) {
            res.vertex_list.emplace_back(i, enc.vertex_labels_seq[i]);
        }
        for (uint32_t i = 0; i < k; ++i) {
            for (uint32_t j = i + 1; j < k; ++j) {
                if (enc.adj_flat[(size_t)i * k + j]) {
                    uint32_t el = enc.edge_labels_seq[(size_t)i * k + j];
                    res.edge_list.emplace_back(i, j, el);
                }
            }
        }

        for (auto& inst : g->instances) {
            res.instances_mapping.push_back(inst.canon_to_old);
        }

        out.push_back(std::move(res));
    }

    for (Group* g : maximal_groups) delete g;
}

} 