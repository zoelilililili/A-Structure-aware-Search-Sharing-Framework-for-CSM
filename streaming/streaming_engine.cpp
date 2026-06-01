#include <cctype>                   
#ifdef _WIN32
#include <windows.h>
#include <dbghelp.h>
#endif
#include <future>                     
#include <fstream>                  
#include <memory>                      
#include <sstream>                     
#include <filesystem>    
#include <map>                        
#include "subgraph_automorphism_lib.h" 
#include "subgraph_isomorphism_lib.h"  
#include <unordered_map>               
#include <unordered_set>              
#include <limits>                      
#include "streaming_engine.h"        
#include "streaming_config.h"
#include "graphoperations.h"           
#include "utility/simple_command_parser.h"  
#include "utility/mem_usage.h"        

namespace {
void build_mapped_new_orbit_peers_by_key_from_index(
        uint32_t k,
        const std::unordered_map<uint32_t, std::vector<uint32_t>> &mapping_ids_index,
        const std::unordered_map<uint32_t, std::vector<uint32_t>> &vertex_permutations,
        std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>> &out_orbit_peers) {
    out_orbit_peers.clear();
    if (k == 0) return;
    out_orbit_peers.reserve(mapping_ids_index.size());
    for (const auto &key_item : mapping_ids_index) {
        const auto &candidate_mapping_ids = key_item.second;
        std::vector<uint32_t> parent(k);
        for (uint32_t i = 0; i < k; ++i)
            parent[i] = i;
        auto find_root = [&parent](uint32_t x) {
            while (parent[x] != x)
                x = parent[x];
            return x;
        };
        auto unite = [&parent, &find_root](uint32_t a, uint32_t b) {
            a = find_root(a);
            b = find_root(b);
            if (a != b)
                parent[a] = b;
        };
        for (uint32_t mapping_id : candidate_mapping_ids) {
            auto it_perm = vertex_permutations.find(mapping_id);
            if (it_perm == vertex_permutations.end())
                continue;
            const auto &perm = it_perm->second;
            if (perm.size() != k)
                continue;
            for (uint32_t i = 0; i < k; ++i)
                unite(i, perm[i]);
        }
        std::unordered_map<uint32_t, std::vector<uint32_t>> comp;
        for (uint32_t i = 0; i < k; ++i) {
            uint32_t r = find_root(i);
            comp[r].push_back(i);
        }
        std::vector<std::vector<uint32_t>> peers(k);
        for (auto &pr : comp) {
            std::sort(pr.second.begin(), pr.second.end());
            for (uint32_t id : pr.second)
                peers[id] = pr.second;
        }
        out_orbit_peers.emplace(key_item.first, std::move(peers));
    }
}
} 

#ifdef MEASURE_UPDATE_COST
    std::vector<uint64_t> g_measure_time_cost_per_update;
#endif
#ifdef MEASURE_INDEXING_COST
    uint64_t g_update_count;
#endif
void StreamingEngine::initialize(const Graph *query_graph, const Graph *data_graph, uint64_t target_embedding_num) {
    data_graph_ = data_graph;                
    gvm_.initialize(query_graph);            
    lvm_.initialize(query_graph, data_graph); 
    om_.initialize(query_graph);              
    sm_.initialize(query_graph, data_graph);  
    sm_.target_number = target_embedding_num; 
}

void StreamingEngine::release() {
    gvm_.release();          
    lvm_.release();          
    om_.release();            
    sm_.release();             
    subquery_cache_.clear();   
    isosubquery_cache_.clear();
    query_bctree_.clear();
}

void StreamingEngine::preprocess(const Graph *query_graph, const Graph *data_graph) {
    if (phase3_enabled_) {                                     
        preprocess_phase3_on(query_graph, data_graph);
    } else if (phase4_enabled_) {                              
        preprocess_phase4_on(query_graph, data_graph);
    } else if (phase2_enabled_) {                             
        preprocess_phase2_on(query_graph, data_graph);
    } else if (phase5_enabled_) {                             
        preprocess_phase5_on(query_graph, data_graph);
    }else if (phase1_enabled_) {                              
        preprocess_phase1_on(query_graph, data_graph);
    } else {                                                   
        preprocess_phase1_off(query_graph, data_graph);
    }
}

void StreamingEngine::preprocess_phase1_off(const Graph *query_graph, const Graph *data_graph) {
    double initial_vm_mem_usage_KB = 0;
    double initial_rss_mem_usage_KB = 0;
    double vm_mem_usage_KB = 0;
    double rss_mem_usage_KB = 0;
    mem_usage::process_mem_usage(initial_vm_mem_usage_KB, initial_rss_mem_usage_KB);
    auto start = std::chrono::high_resolution_clock::now();
    gvm_.create_views(query_graph, data_graph);
    auto end = std::chrono::high_resolution_clock::now();
    global_view_initialize_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    mem_usage::process_mem_usage(vm_mem_usage_KB, rss_mem_usage_KB);
    start = std::chrono::high_resolution_clock::now();
    om_.create_orders(query_graph);
    end = std::chrono::high_resolution_clock::now();
    order_generation_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
}
//To be revised
void StreamingEngine::preprocess_phase2_on(const Graph *query_graph, const Graph *data_graph) {
    build_subquery_cache_in_memory(query_graph, data_graph);
    build_isosubquery_cache_in_memory(query_graph, data_graph);
    preprocess_phase1_grouping_tags(query_graph);
    preprocess_phase2_grouping_tags(query_graph);
    double initial_vm_mem_usage_KB = 0; 
    double initial_rss_mem_usage_KB = 0; 
    double vm_mem_usage_KB = 0;          
    double rss_mem_usage_KB = 0;         
    mem_usage::process_mem_usage(initial_vm_mem_usage_KB, initial_rss_mem_usage_KB);
    auto start = std::chrono::high_resolution_clock::now();
    gvm_.create_views(query_graph, data_graph);
    auto end = std::chrono::high_resolution_clock::now();
    global_view_initialize_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    mem_usage::process_mem_usage(vm_mem_usage_KB, rss_mem_usage_KB);
    start = std::chrono::high_resolution_clock::now();
    om_.create_orders(query_graph);
    if (!subquery_cache_.empty()) {
        std::vector<std::vector<int32_t>> sub_new2old(subquery_cache_.size());
        for (const auto& kv : subquery_cache_) {
            const uint32_t sub_id = kv.first;
            const SubqueryContext& sc = kv.second;
            if (sub_id >= sub_new2old.size()) continue;
            sub_new2old[sub_id] = sc.new2old; 
        }
        om_.create_special_orders(query_graph, sub_new2old);
    }

    if (!isosubquery_cache_.empty()) {
        uint32_t max_iso_id = 0;
        for (const auto& kv : isosubquery_cache_) {
            if (kv.first > max_iso_id) max_iso_id = kv.first;
        }
        std::vector<std::vector<std::vector<int32_t>>> iso_instance_vertices(max_iso_id + 1);
        for (const auto& kv : isosubquery_cache_) {
            const uint32_t group_id = kv.first;
            const isoSubqueryContext& sc = kv.second;
            if (group_id >= iso_instance_vertices.size()) continue;
            iso_instance_vertices[group_id].reserve(sc.mappings.size());
            for (const auto& mapping : sc.mappings) {
                std::vector<int32_t> old_vertices;
                old_vertices.reserve(mapping.size());
                for (uint32_t old_u : mapping) {
                    if (old_u < query_graph->getVerticesCount()) old_vertices.push_back(static_cast<int32_t>(old_u));
                }
                iso_instance_vertices[group_id].push_back(std::move(old_vertices));
            }
        }
        om_.create_special_orders2(query_graph, iso_instance_vertices);
    }
    end = std::chrono::high_resolution_clock::now();
    order_generation_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
}

void StreamingEngine::preprocess_phase4_on(const Graph *query_graph, const Graph *data_graph) {

    build_subquery_cache_in_memory(query_graph, data_graph);
    preprocess_phase1_grouping_tags(query_graph);
    double initial_vm_mem_usage_KB = 0;  
    double initial_rss_mem_usage_KB = 0; 
    double vm_mem_usage_KB = 0;          
    double rss_mem_usage_KB = 0;         

    mem_usage::process_mem_usage(initial_vm_mem_usage_KB, initial_rss_mem_usage_KB);
    auto start = std::chrono::high_resolution_clock::now();
    gvm_.create_views(query_graph, data_graph);

    auto end = std::chrono::high_resolution_clock::now();
    global_view_initialize_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    mem_usage::process_mem_usage(vm_mem_usage_KB, rss_mem_usage_KB);
    start = std::chrono::high_resolution_clock::now();
    om_.create_orders(query_graph);
    if (!subquery_cache_.empty()) {
        std::vector<std::vector<int32_t>> sub_new2old(subquery_cache_.size());
        for (const auto& kv : subquery_cache_) {
            const uint32_t sub_id = kv.first;
            const SubqueryContext& sc = kv.second;
            if (sub_id >= sub_new2old.size()) continue;
            sub_new2old[sub_id] = sc.new2old; 
        }
        om_.create_special_orders(query_graph, sub_new2old);
    }

    end = std::chrono::high_resolution_clock::now();
    order_generation_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    query_bctree_.clear();
    query_bctree_ = build_query_bctree_decomposition(query_graph);
    om_.create_block_reuse_orders(query_graph, query_bctree_);
}
//To be revised
void StreamingEngine::preprocess_phase5_on(const Graph *query_graph, const Graph *data_graph) {
    build_subquery_cache_in_memory(query_graph, data_graph);
    build_isosubquery_cache_in_memory(query_graph, data_graph);
    preprocess_phase1_grouping_tags(query_graph);
    preprocess_phase2_grouping_tags(query_graph);

    double initial_vm_mem_usage_KB = 0; 
    double initial_rss_mem_usage_KB = 0;
    double vm_mem_usage_KB = 0;         
    double rss_mem_usage_KB = 0;        

    mem_usage::process_mem_usage(initial_vm_mem_usage_KB, initial_rss_mem_usage_KB);
    auto start = std::chrono::high_resolution_clock::now();
    gvm_.create_views(query_graph, data_graph);
    auto end = std::chrono::high_resolution_clock::now();
    global_view_initialize_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    mem_usage::process_mem_usage(vm_mem_usage_KB, rss_mem_usage_KB);
    start = std::chrono::high_resolution_clock::now();
    om_.create_orders(query_graph);
    if (!subquery_cache_.empty()) {
        std::vector<std::vector<int32_t>> sub_new2old(subquery_cache_.size());
        for (const auto& kv : subquery_cache_) {
            const uint32_t sub_id = kv.first;
            const SubqueryContext& sc = kv.second;
            if (sub_id >= sub_new2old.size()) continue;
            sub_new2old[sub_id] = sc.new2old; 
        }
        om_.create_special_orders(query_graph, sub_new2old);
    }
     if (!isosubquery_cache_.empty()) {
        uint32_t max_iso_id = 0;
        for (const auto& kv : isosubquery_cache_) {
            if (kv.first > max_iso_id) max_iso_id = kv.first;
        }
        std::vector<std::vector<std::vector<int32_t>>> iso_instance_vertices(max_iso_id + 1);
        for (const auto& kv : isosubquery_cache_) {
            const uint32_t group_id = kv.first;
            const isoSubqueryContext& sc = kv.second;
            if (group_id >= iso_instance_vertices.size()) continue;
            iso_instance_vertices[group_id].reserve(sc.mappings.size());
            for (const auto& mapping : sc.mappings) {
                std::vector<int32_t> old_vertices;
                old_vertices.reserve(mapping.size());
                for (uint32_t old_u : mapping) {
                    if (old_u < query_graph->getVerticesCount()) old_vertices.push_back(static_cast<int32_t>(old_u));
                }
                iso_instance_vertices[group_id].push_back(std::move(old_vertices));
            }
        }
        om_.create_special_orders2(query_graph, iso_instance_vertices);
    }
    end = std::chrono::high_resolution_clock::now();
    order_generation_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    query_bctree_.clear();
    query_bctree_ = build_query_bctree_decomposition(query_graph);
    om_.create_block_reuse_orders(query_graph, query_bctree_);
}

void StreamingEngine::preprocess_phase3_on(const Graph *query_graph, const Graph *data_graph) {
    double initial_vm_mem_usage_KB = 0;
    double initial_rss_mem_usage_KB = 0;
    double vm_mem_usage_KB = 0;
    double rss_mem_usage_KB = 0;
    mem_usage::process_mem_usage(initial_vm_mem_usage_KB, initial_rss_mem_usage_KB);
    auto start = std::chrono::high_resolution_clock::now();
    gvm_.create_views(query_graph, data_graph);
    auto end = std::chrono::high_resolution_clock::now();

    global_view_initialize_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    mem_usage::process_mem_usage(vm_mem_usage_KB, rss_mem_usage_KB);
    start = std::chrono::high_resolution_clock::now();

    om_.create_orders(query_graph);
    end = std::chrono::high_resolution_clock::now();
    order_generation_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    query_bctree_.clear();
    query_bctree_ = build_query_bctree_decomposition(query_graph);
    om_.create_block_reuse_orders(query_graph, query_bctree_);
}

void StreamingEngine::preprocess_phase1_on(const Graph *query_graph, const Graph *data_graph) {
    build_subquery_cache_in_memory(query_graph, data_graph);
    preprocess_phase1_grouping_tags(query_graph);
    double initial_vm_mem_usage_KB = 0; 
    double initial_rss_mem_usage_KB = 0; 
    double vm_mem_usage_KB = 0;          
    double rss_mem_usage_KB = 0;       

    mem_usage::process_mem_usage(initial_vm_mem_usage_KB, initial_rss_mem_usage_KB);
    auto start = std::chrono::high_resolution_clock::now();
    gvm_.create_views(query_graph, data_graph);
    auto end = std::chrono::high_resolution_clock::now();
    global_view_initialize_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    mem_usage::process_mem_usage(vm_mem_usage_KB, rss_mem_usage_KB);
    start = std::chrono::high_resolution_clock::now();
    om_.create_orders(query_graph);
    if (!subquery_cache_.empty()) {
        std::vector<std::vector<int32_t>> sub_new2old(subquery_cache_.size());
        for (const auto& kv : subquery_cache_) {
            const uint32_t sub_id = kv.first;
            const SubqueryContext& sc = kv.second;
            if (sub_id >= sub_new2old.size()) continue;
            sub_new2old[sub_id] = sc.new2old; // new(sub) -> old(query) 映射，后续 special matching order 用它来确定 forced old 顶点集合
        }
        om_.create_special_orders(query_graph, sub_new2old);
    }
    end = std::chrono::high_resolution_clock::now();
    order_generation_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
}

void StreamingEngine::preprocess_phase1_grouping_tags(const Graph* query_graph) {
    phase1_edge_reuse_tag_.clear();
    phase1_edge_tag_atoms_.clear();
    phase1_edge_tag_atoms_.clear();
    phase1_edge_group_id_.clear();
    phase1_group_chosen_sub_id_.clear();
    phase1_group_chosen_sub_id_.push_back(std::numeric_limits<uint32_t>::max()); 

    if (!query_graph) return;
    if (subquery_cache_.empty()) return;

    const uint32_t qn = query_graph->getVerticesCount();
    if (qn == 0) return;

    auto build_atoms_for_directed_edge = [&](uint32_t u0, uint32_t u1, uint32_t lbl) {
        std::vector<Phase1TagAtom> atoms;
        for (const auto& kv : subquery_cache_) {
            const uint32_t sub_id = kv.first;
            const SubqueryContext& sc = kv.second;
            if (!sc.initialized || !sc.sub_query) continue;
            if (u0 >= sc.old2new.size() || u1 >= sc.old2new.size()) continue;

            const int32_t n0 = sc.old2new[u0];
            const int32_t n1 = sc.old2new[u1];
            if (n0 < 0 || n1 < 0) continue; 
            const uint32_t nu = static_cast<uint32_t>(std::min(n0, n1));
            const uint32_t nv = static_cast<uint32_t>(std::max(n0, n1));
            const uint32_t el = sc.sub_query->getEdgeLabelByVertex(nu, nv);
            if (el == std::numeric_limits<uint32_t>::max()) continue;
            (void)lbl; 

            bool movable = false;
            const uint32_t nn0 = static_cast<uint32_t>(n0);
            const uint32_t nn1 = static_cast<uint32_t>(n1);
            for (const auto& pv : sc.vertex_permutations) {
                const auto& perm = pv.second;
                if (nn0 >= perm.size() || nn1 >= perm.size()) continue;
                if (perm[nn0] != nn0 || perm[nn1] != nn1) {
                    movable = true;
                    break;
                }
            }
            if (!movable) continue;

            if (u0 >= sc.old_vertex_rep_new.size() || u1 >= sc.old_vertex_rep_new.size()) continue;
            const int32_t r0 = sc.old_vertex_rep_new[u0];
            const int32_t r1 = sc.old_vertex_rep_new[u1];
            if (r0 < 0 || r1 < 0) continue;

            atoms.push_back(Phase1TagAtom{sub_id, static_cast<uint32_t>(r0), static_cast<uint32_t>(r1)});
        }

        std::sort(atoms.begin(), atoms.end(), Phase1TagAtomLess{});
        atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
        return atoms;
    };

    uint64_t undirected_edge_count = 0;
    uint64_t directed_edge_count = 0;
    uint64_t shareable_directed_edge_count = 0;
    uint64_t total_atom_count = 0;

    for (uint32_t a = 0; a < qn; ++a) {
        uint32_t cnt = 0;
        const ui* nbrs = query_graph->getVertexNeighbors(a, cnt);
        for (uint32_t i = 0; i < cnt; ++i) {
            const uint32_t b = nbrs[i];
            if (b >= qn) continue;
            if (a >= b) continue; 

            const uint32_t lbl = 0;
            undirected_edge_count += 1;
            {
                DirectedEdgeKey key{a, b, lbl};
                auto atoms = build_atoms_for_directed_edge(a, b, lbl);
                Phase1EdgeReuseTag tag;
                tag.shareable = !atoms.empty();
                // tag_set_id 先占位为 0；后续如需稳定 id，可对 atoms hash 后赋值。
                phase1_edge_reuse_tag_[key] = tag;
                total_atom_count += atoms.size();
                if (tag.shareable) shareable_directed_edge_count += 1;
                directed_edge_count += 1;
                phase1_edge_tag_atoms_.emplace(key, std::move(atoms));
            }
            {
                DirectedEdgeKey key{b, a, lbl};
                auto atoms = build_atoms_for_directed_edge(b, a, lbl);
                Phase1EdgeReuseTag tag;
                tag.shareable = !atoms.empty();
                phase1_edge_reuse_tag_[key] = tag;
                total_atom_count += atoms.size();
                if (tag.shareable) shareable_directed_edge_count += 1;
                directed_edge_count += 1;
                phase1_edge_tag_atoms_.emplace(key, std::move(atoms));
            }
        }
    }

    auto intersect_atoms = [&](const std::vector<Phase1TagAtom>& A, const std::vector<Phase1TagAtom>& B) {
        std::vector<Phase1TagAtom> out;                       
        out.reserve(std::min(A.size(), B.size()));             
        size_t i = 0, j = 0;                                    
        Phase1TagAtomLess less;                                
        while (i < A.size() && j < B.size()) {                
            if (less(A[i], B[j])) { ++i; continue; }          
            if (less(B[j], A[i])) { ++j; continue; }           
            out.push_back(A[i]);                               
            ++i; ++j;                                          
        }
        return out;                                             
    };

    struct TempGroup {
        std::vector<DirectedEdgeKey> edges;                   
        std::vector<Phase1TagAtom> common;                     
    };
    std::unordered_map<uint32_t, std::vector<DirectedEdgeKey>> by_label;
    by_label.reserve(64);
  
    for (const auto& kv : phase1_edge_tag_atoms_) {             
        const DirectedEdgeKey& key = kv.first;                 
        const auto& atoms = kv.second;                          
        if (atoms.empty()) continue;                          
        by_label[key.label].push_back(key);                     
    }

    uint32_t next_group_id = 1;                                 
    for (auto& lb : by_label) {                                 
        auto& keys = lb.second;                                
        std::vector<TempGroup> groups;                          
        groups.reserve(keys.size() / 2 + 1);                    

        for (const auto& key : keys) {                          
            const auto it = phase1_edge_tag_atoms_.find(key);   
            if (it == phase1_edge_tag_atoms_.end()) continue;   
            const auto& tags = it->second;                      
            if (tags.empty()) continue;                        

            bool merged = false;                               
            for (auto& g : groups) {                            
                auto inter = intersect_atoms(g.common, tags);  
                if (!inter.empty()) {                          
                    g.edges.push_back(key);                     
                    g.common = std::move(inter);                
                    merged = true;                             
                    break;                                      
                }
            }
            if (!merged) {                                     
                TempGroup g;                                   
                g.edges.push_back(key);                        
                g.common = tags;                               
                groups.push_back(std::move(g));                 
            }
        }

        for (auto& g : groups) {                                
            if (g.edges.size() <= 1 || g.common.empty()) {       
                continue;
            }
            uint32_t chosen_sub_id = std::numeric_limits<uint32_t>::max(); 
            for (const auto& a : g.common) {                    
                chosen_sub_id = std::min(chosen_sub_id, a.sub_id);
            }
            const uint32_t gid = next_group_id++;                
            if (phase1_group_chosen_sub_id_.size() <= gid) {      
                phase1_group_chosen_sub_id_.resize(gid + 1, std::numeric_limits<uint32_t>::max());
            }
            phase1_group_chosen_sub_id_[gid] = chosen_sub_id;     
            for (const auto& key : g.edges) {                    
                phase1_edge_group_id_[key] = gid;
            }
        }
    }
}
//first band
void StreamingEngine::preprocess_phase2_grouping_tags(const Graph* query_graph) {

    phase2_edge_group_id_.clear();
    phase2_group_chosen_sub_id_.clear();
    phase2_group_chosen_sub_id_.push_back(std::numeric_limits<uint32_t>::max()); // group_id=0 占位

    if (!query_graph) return;
    if (isosubquery_cache_.empty()) return;

    const uint32_t qn = query_graph->getVerticesCount();
    if (qn == 0) return;

    uint32_t next_group_id = 1;

    for (const auto& kv : isosubquery_cache_) {
        const uint32_t iso_group_id = kv.first;
        const isoSubqueryContext& sc = kv.second;
        if (!sc.initialized || !sc.sub_query) continue;
        if (sc.mappings.size() <= 1) continue;

        const uint32_t k = sc.sub_query->getVerticesCount();

        for (uint32_t cu = 0; cu < k; ++cu) {
            uint32_t deg = 0;
            const uint32_t* nbrs = sc.sub_query->getVertexNeighbors(cu, deg);
            for (uint32_t ei = 0; ei < deg; ++ei) {
                const uint32_t cv = nbrs[ei];
                if (cv >= k) continue;
                if (cu >= cv) continue;

                const uint32_t edge_label = sc.sub_query->getEdgeLabelByVertex(cu, cv);

                // 方向1：canon (cu → cv) → old (mapping[cu] → mapping[cv])
                {
                    std::vector<DirectedEdgeKey> fwd_edges;
                    fwd_edges.reserve(sc.mappings.size());
                    for (const auto& mapping : sc.mappings) {
                        if (cu >= mapping.size() || cv >= mapping.size()) continue;
                        uint32_t old_u = mapping[cu];
                        uint32_t old_v = mapping[cv];
                        if (old_u >= qn || old_v >= qn) continue;
                        fwd_edges.push_back(DirectedEdgeKey{old_u, old_v, edge_label});
                    }
                    std::sort(fwd_edges.begin(), fwd_edges.end(),
                              [](const DirectedEdgeKey& a, const DirectedEdgeKey& b) {
                                  if (a.src != b.src) return a.src < b.src;
                                  if (a.dst != b.dst) return a.dst < b.dst;
                                  return a.label < b.label;
                              });
                    fwd_edges.erase(std::unique(fwd_edges.begin(), fwd_edges.end()), fwd_edges.end());
                    if (fwd_edges.size() > 1) {
                        const uint32_t gid = next_group_id++;
                        if (phase2_group_chosen_sub_id_.size() <= gid) {
                            phase2_group_chosen_sub_id_.resize(gid + 1, std::numeric_limits<uint32_t>::max());
                        }
                        phase2_group_chosen_sub_id_[gid] = iso_group_id;
                        for (const auto& key : fwd_edges) {
                            phase2_edge_group_id_[key] = gid;
                        }
                    }
                }
                {
                    std::vector<DirectedEdgeKey> rev_edges;
                    rev_edges.reserve(sc.mappings.size());
                    for (const auto& mapping : sc.mappings) {
                        if (cu >= mapping.size() || cv >= mapping.size()) continue;
                        uint32_t old_u = mapping[cu];
                        uint32_t old_v = mapping[cv];
                        if (old_u >= qn || old_v >= qn) continue;
                        rev_edges.push_back(DirectedEdgeKey{old_v, old_u, edge_label});
                    }
                    std::sort(rev_edges.begin(), rev_edges.end(),
                              [](const DirectedEdgeKey& a, const DirectedEdgeKey& b) {
                                  if (a.src != b.src) return a.src < b.src;
                                  if (a.dst != b.dst) return a.dst < b.dst;
                                  return a.label < b.label;
                              });
                    rev_edges.erase(std::unique(rev_edges.begin(), rev_edges.end()), rev_edges.end());
                    if (rev_edges.size() > 1) {
                        const uint32_t gid = next_group_id++;
                        if (phase2_group_chosen_sub_id_.size() <= gid) {
                            phase2_group_chosen_sub_id_.resize(gid + 1, std::numeric_limits<uint32_t>::max());
                        }
                        phase2_group_chosen_sub_id_[gid] = iso_group_id;
                        for (const auto& key : rev_edges) {
                            phase2_edge_group_id_[key] = gid;
                        }
                    }
                }
            }
        }
    }
}

StreamingEngine::SubqueryContext* StreamingEngine::get_subquery_context(uint32_t sub_id) {
    auto it = subquery_cache_.find(sub_id);
    if (it == subquery_cache_.end()) return nullptr;
    return &it->second;
}

StreamingEngine::isoSubqueryContext* StreamingEngine::get_isosubquery_context(uint32_t group_id) {
    auto it = isosubquery_cache_.find(group_id);
    if (it == isosubquery_cache_.end()) return nullptr;
    return &it->second;
}

//first band
void StreamingEngine::build_isosubquery_cache_in_memory(const Graph* query_graph, const Graph* data_graph, uint64_t max_isos) {

    isosubquery_cache_.clear();
    if (!query_graph || !data_graph) return;
    if (query_graph->getVerticesCount() > 63) {
        return;
    }
    std::vector<subgraph_isomorphism_lib::SubIsomorphismGroupResult> subs;
    subgraph_isomorphism_lib::build_maximal_isomorphic_induced_subgraphs_in_memory(*query_graph, subs);
    if (subs.empty()) return;

    const uint32_t qn = query_graph->getVerticesCount();

    for (const auto& r : subs) {
        const uint32_t group_id = r.group_id;
        auto emplace_result = isosubquery_cache_.emplace(std::piecewise_construct,
                                                         std::forward_as_tuple(group_id),
                                                         std::forward_as_tuple());
        isoSubqueryContext& ctx = emplace_result.first->second;

        const uint32_t k = static_cast<uint32_t>(r.vertex_list.size());
        if (k == 0 || r.instances_mapping.empty()) continue;

        const auto& first_mapping = r.instances_mapping[0];

        ctx.old2new.assign(qn, -1);
        for (uint32_t canon_id = 0; canon_id < k; ++canon_id) {
            uint32_t old_id = first_mapping[canon_id];
            if (old_id >= qn) continue;
            ctx.old2new[old_id] = static_cast<int32_t>(canon_id);
        }

        ctx.mappings = r.instances_mapping;

        ctx.sub_query = std::make_unique<Graph>(false);
        ctx.sub_query->is_edge_labeled = true;
        {
            auto vertex_list = r.vertex_list;
            auto edge_list = r.edge_list;
            ctx.sub_query->loadGraphFromMemory(vertex_list, edge_list);
        }

        ctx.om.initialize(ctx.sub_query.get());
        ctx.om.create_orders_no_automorphism(ctx.sub_query.get());
        ctx.lvm.initialize(ctx.sub_query.get(), data_graph);
        ctx.sm.initialize(ctx.sub_query.get(), data_graph);
        ctx.sm.target_number = sm_.target_number;
        ctx.initialized = true;
    }
}

void StreamingEngine::build_subquery_cache_in_memory(const Graph* query_graph, const Graph* data_graph, uint64_t max_autos) {
    subquery_cache_.clear();
    if (!query_graph || !data_graph) return;

    if (query_graph->getVerticesCount() > 63) {
        return;
    }
    std::vector<subgraph_automorphism_lib::SubAutomorphismResult> subs;
    subgraph_automorphism_lib::build_maximal_automorphic_induced_subgraphs_in_memory(*query_graph, max_autos, subs);
    if (subs.empty()) return;

    const uint32_t qn = query_graph->getVerticesCount();

    for (const auto& r : subs) {
        const uint32_t sub_id = r.sub_id;
        auto emplace_result = subquery_cache_.emplace(std::piecewise_construct,
                                                      std::forward_as_tuple(sub_id),
                                                      std::forward_as_tuple());
        SubqueryContext& ctx = emplace_result.first->second;
        const uint32_t k = static_cast<uint32_t>(r.old_vertices_sorted.size());
        if (k == 0) continue;
        ctx.old2new.assign(qn, -1);
        ctx.new2old.assign(k, -1);
        for (uint32_t new_id = 0; new_id < k; ++new_id) {
            uint32_t old_id = r.old_vertices_sorted[new_id];
            if (old_id >= qn) continue;
            ctx.old2new[old_id] = static_cast<int32_t>(new_id);
            ctx.new2old[new_id] = static_cast<int32_t>(old_id);
        }
        ctx.old_vertex_rep_new.assign(qn, -1);
        for (uint32_t new_id = 0; new_id < k; ++new_id) {
            uint32_t rep = new_id;
            for (const auto& perm : r.automorphisms) {
                if (perm.size() != k) continue;
                rep = std::min(rep, perm[new_id]);
            }
            uint32_t old_id = r.old_vertices_sorted[new_id];
            if (old_id < qn) {
                ctx.old_vertex_rep_new[old_id] = static_cast<int32_t>(rep);
            }
        }
        ctx.vertex_permutations.clear();
        for (uint32_t aid = 0; aid < r.automorphisms.size(); ++aid) {
            if (r.automorphisms[aid].size() == k) {
                ctx.vertex_permutations[aid] = r.automorphisms[aid];
            }
        }

        ctx.sub_old_vertices.clear();                                                   
        for (uint32_t old_u = 0; old_u < ctx.old2new.size(); ++old_u) {                  
            if (ctx.old2new[old_u] >= 0) {                                                 
                ctx.sub_old_vertices.push_back(old_u);                                    
            }
        }
        ctx.forced_pairs_template.clear();                                                
        ctx.forced_pairs_template.reserve(ctx.vertex_permutations.size());               
        for (const auto& mp : ctx.vertex_permutations) {                                   
            const uint32_t mapping_id = mp.first;                                          
            const auto& perm = mp.second;                                                   
            if (perm.size() != k) continue;                                                
            std::vector<std::pair<uint32_t, uint32_t>> tpl;                                
            tpl.reserve(ctx.sub_old_vertices.size());                                      
            for (uint32_t old_u : ctx.sub_old_vertices) {                                   
                const uint32_t new_u = static_cast<uint32_t>(ctx.old2new[old_u]);           
                const uint32_t mapped_new_u = perm[new_u];                                  
                tpl.emplace_back(old_u, mapped_new_u);                                     
            }
            ctx.forced_pairs_template.emplace(mapping_id, std::move(tpl));                  
        }
    
        ctx.mapping_ids_index.clear();                                                     
        ctx.mapping_ids_index.reserve(ctx.vertex_permutations.size() * 4);                
        ctx.mapped_new_set_by_key_and_old.clear();
        ctx.mapped_new_orbit_peers_by_key.clear();
        auto pack4 = [](uint32_t a, uint32_t b, uint32_t ns, uint32_t nt) -> uint32_t {   
            return (a & 0xFFu) | ((b & 0xFFu) << 8) | ((ns & 0xFFu) << 16) | ((nt & 0xFFu) << 24);
        };
        for (const auto& mp : ctx.vertex_permutations) {                                   
            const uint32_t mapping_id = mp.first;                                          
            const auto& perm = mp.second;                                                 
            if (perm.size() != k) continue;                                                
            for (uint32_t task_u0 : r.old_vertices_sorted) {                              
                for (uint32_t task_u1 : r.old_vertices_sorted) {
                    if (task_u0 == task_u1) continue;                                       
                    if (task_u0 >= ctx.old2new.size() || task_u1 >= ctx.old2new.size()) continue;
                    const int32_t task_n0_i = ctx.old2new[task_u0];
                    const int32_t task_n1_i = ctx.old2new[task_u1];
                    if (task_n0_i < 0 || task_n1_i < 0) continue;                         
                    const uint32_t task_n0 = static_cast<uint32_t>(task_n0_i);
                    const uint32_t task_n1 = static_cast<uint32_t>(task_n1_i);
                    const uint32_t ns = perm[task_n0];                                     
                    const uint32_t nt = perm[task_n1];                                     
                    const uint32_t key = pack4(task_u0, task_u1, ns, nt);                  
                    ctx.mapping_ids_index[key].push_back(mapping_id);                       
                }
            }
        }
        ctx.mapped_new_set_by_key_and_old.reserve(ctx.mapping_ids_index.size());
        for (const auto& key_item : ctx.mapping_ids_index) {
            const uint32_t key = key_item.first;
            const auto& candidate_mapping_ids = key_item.second;
            auto& per_old = ctx.mapped_new_set_by_key_and_old[key];
            per_old.assign(qn, {});
            for (uint32_t old_u : ctx.sub_old_vertices) {
                const int32_t new_u_i = ctx.old2new[old_u];
                if (new_u_i < 0) continue;
                const uint32_t new_u = static_cast<uint32_t>(new_u_i);
                std::vector<uint8_t> seen(k, 0); 
                auto& mapped_set = per_old[old_u];
                for (uint32_t mapping_id : candidate_mapping_ids) {
                    auto it_perm = ctx.vertex_permutations.find(mapping_id);
                    if (it_perm == ctx.vertex_permutations.end()) continue;
                    const auto& perm = it_perm->second;
                    if (new_u >= perm.size()) continue;
                    const uint32_t mapped_new_u = perm[new_u];
                    if (mapped_new_u >= k) continue;
                    if (!seen[mapped_new_u]) {
                        seen[mapped_new_u] = 1;
                        mapped_set.push_back(mapped_new_u);
                    }
                }
            }
        }

        build_mapped_new_orbit_peers_by_key_from_index(k, ctx.mapping_ids_index, ctx.vertex_permutations,
                                                     ctx.mapped_new_orbit_peers_by_key);
        ctx.sub_query = std::make_unique<Graph>(false);
        ctx.sub_query->is_edge_labeled = true;
        {
            auto vertex_list = r.vertex_list;
            auto edge_list = r.edge_list;
            ctx.sub_query->loadGraphFromMemory(vertex_list, edge_list);
        }
        ctx.om.initialize(ctx.sub_query.get());
        ctx.om.create_orders_no_automorphism(ctx.sub_query.get());
        ctx.lvm.initialize(ctx.sub_query.get(), data_graph);
        ctx.sm.initialize(ctx.sub_query.get(), data_graph);
        ctx.sm.target_number = sm_.target_number;
        ctx.initialized = true;
    }
}

uint64_t StreamingEngine::execute(const Graph *query_graph, const Update &update, bool enable_local_view,
                                  bool enable_search) {
    if (phase3_enabled_) {
        return execute_phase3_on(query_graph, update, enable_local_view, enable_search);
    }
    if (phase4_enabled_) {
        return execute_phase4_on(query_graph, update, enable_local_view, enable_search);
    }
    if (phase2_enabled_) {
        return execute_phase2_on(query_graph, update, enable_local_view, enable_search);
    }
    if (phase5_enabled_) {
        return execute_phase5_on(query_graph, update, enable_local_view, enable_search);
    }
    if (phase1_enabled_) {
        return execute_phase1_on(query_graph, update, enable_local_view, enable_search);
    }
    return execute_phase1_off(query_graph, update, enable_local_view, enable_search);
}

//To be revised
uint64_t StreamingEngine::execute_phase2_on(const Graph *query_graph, const Update &update, bool enable_local_view,
                                           bool enable_search) {
   
    LabelTriple label_triple = update.labels_; 
    Edge data_edge = update.edge_; 
    
    const LabelTriple label_triple_initial = label_triple;
    uint64_t result_count = 0; 
    is_relevant_ = false; 
    is_searched_ = false;
    if (update.op_ == '+') {       
        auto it = gvm_.get_mapped_views(label_triple); 
        if (it != nullptr) { 
            is_relevant_ = true; 
            relevant_update_count_ += 1; 
            auto t_gvm_begin = std::chrono::high_resolution_clock::now();
            for (uint32_t i = 0; i < 2; ++i) { 
                it = gvm_.get_mapped_views(label_triple); 
                uint32_t view_id = it->second; 
                gvm_.update_view(update.op_, view_id, data_edge, label_triple); 
                // 处理无向边的反向键：(v,u) 与交换后的标签三元组
                std::swap(data_edge.first, data_edge.second);
                std::swap(label_triple.src_label_, label_triple.dst_label_); 
            }
            std::swap(data_edge.first, data_edge.second); 
            std::swap(label_triple.src_label_, label_triple.dst_label_); 
            auto t_gvm_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_views_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_gvm_end - t_gvm_begin).count();
        }
        if (is_relevant_) { 
            auto t_nlf_begin = std::chrono::high_resolution_clock::now();
            gvm_.update_nlf_view(update.op_, data_edge, label_triple); 
            auto t_nlf_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_nlf_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_nlf_end - t_nlf_begin).count();
        }
   
        std::vector<Edge> *mapped_query_edges = om_.get_mapped_edges(label_triple);
        const size_t phase1_matching_query_edge_count = mapped_query_edges ? mapped_query_edges->size() : 0;
        auto mapped_automorphism = om_.get_mapped_automorphism(label_triple);
        if (mapped_automorphism != nullptr) { 
            is_relevant_ = true; 
            if (enable_local_view) { 
                if (phase1_matching_query_edge_count <= 1) { 
                    for (auto automorphism_id: *mapped_automorphism) { 
                        auto meta = om_.get_automorphism_meta(automorphism_id); 
                        OrdersPerEdge *orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta); 
                        auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                        bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge); 
                        auto t_lvm_end = std::chrono::high_resolution_clock::now();
                        time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_end - t_lvm_begin).count();
                        if (is_valid) { 
                            search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_; 
                            search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                        } else { 
                            non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_; 
                            non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_; 
                            if (lvm_.generate_visited_neighbor_count_ == 0) 
                                direct_rejection_count_ += 1; 
                        }
                        first_indexing_vertex_ += lvm_.first_vertex_neighbor_; 
                        if (enable_search && is_valid) { 
                            auto t_plain_begin = std::chrono::high_resolution_clock::now();
                            uint64_t local_result_count = sm_.search_on_reduced_query(query_graph, *orders, lvm_, gvm_); 
                            auto t_plain_end = std::chrono::high_resolution_clock::now();
                            time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_plain_end - t_plain_begin).count();
                            result_count += local_result_count * automorphism_size;
                            is_searched_ = true; 
                        }
                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view(); 
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_destroy_end - t_lvm_destroy_begin).count();

                        if (result_count >= sm_.target_number) break;
                        if (g_exit) break; 
                    }
                }else {
                    struct PreGroupTask {
                        OrdersPerEdge* orders = nullptr;
                        uint32_t automorphism_size = 1;
                        uint32_t automorphism_id = 0;
                    };
                    struct PreGroupBucket {
                        uint32_t group_id = 0;
                        std::vector<PreGroupTask> tasks;
                        uint32_t chosen_sub_id = std::numeric_limits<uint32_t>::max();
                    };

                    std::vector<PreGroupTask> plain_search_pre;
                    std::unordered_map<uint32_t, PreGroupBucket> grouped_pre;
                    std::unordered_map<uint32_t, PreGroupBucket> iso_grouped_pre;
                    plain_search_pre.reserve(mapped_automorphism->size());
                    grouped_pre.reserve(32);
                    iso_grouped_pre.reserve(32);

                    for (auto automorphism_id : *mapped_automorphism) {
                        auto meta = om_.get_automorphism_meta(automorphism_id);
                        OrdersPerEdge* orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta);
                        if (orders == nullptr) continue;

                        const uint32_t u0 = orders->indexing_order_[0];
                        const uint32_t u1 = orders->indexing_order_[1];

                        if (!gvm_.nlf_check(u0, data_edge.first) || !gvm_.nlf_check(u1, data_edge.second)) {
                            continue;
                        }

                        DirectedEdgeKey key{u0, u1, 0};
                        uint32_t gid = 0;
                        auto it_gid = phase1_edge_group_id_.find(key);
                        if (it_gid != phase1_edge_group_id_.end()) {
                             gid = it_gid->second;
                        }

                        if (gid == 0) {
                            uint32_t iso_gid = 0;
                            auto it_iso = phase2_edge_group_id_.find(key);
                            if (it_iso != phase2_edge_group_id_.end()) {
                                iso_gid = it_iso->second;
                            }
                            if (iso_gid >= 1) {
                                auto& ib = iso_grouped_pre[iso_gid];
                                ib.group_id = iso_gid;
                                ib.tasks.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                                ib.chosen_sub_id = phase2_group_chosen_sub_id_[iso_gid];
                                continue;
                            }
                            plain_search_pre.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                            continue;
                        }

                        auto& b = grouped_pre[gid];
                        b.group_id = gid;
                        b.tasks.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                        b.chosen_sub_id = phase1_group_chosen_sub_id_[gid]; 
                    }

                    for (auto it = grouped_pre.begin(); it != grouped_pre.end(); ) {
                        if (it->second.tasks.size() <= 1) {
                            for (auto& t : it->second.tasks)
                                plain_search_pre.push_back(t);
                            it = grouped_pre.erase(it);
                        } else {
                            ++it;
                        }
                    }
                    for (auto it = iso_grouped_pre.begin(); it != iso_grouped_pre.end(); ) {
                        if (it->second.tasks.size() <= 1) {
                            for (auto& t : it->second.tasks)
                                plain_search_pre.push_back(t);
                            it = iso_grouped_pre.erase(it);
                        } else {
                            ++it;
                        }
                    }

                    const auto t_plain_total_begin = std::chrono::high_resolution_clock::now();
                    for (const auto &pt : plain_search_pre) {
                        if (pt.orders == nullptr) continue;
                        // 这里走“单边组”回溯：直接按当前 orders 构建局部视图并搜索（不依赖 pending candidates）。
                        auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                        bool is_valid = lvm_.create_view(query_graph, *pt.orders, gvm_, data_edge);
                        auto t_lvm_end = std::chrono::high_resolution_clock::now();
                        time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                t_lvm_end - t_lvm_begin).count();
                        if (enable_search && is_valid) {
                            auto t_plain_begin = std::chrono::high_resolution_clock::now();
                            uint64_t local_result_count =
                                    sm_.search_on_reduced_query(query_graph, *pt.orders, lvm_, gvm_);
                            auto t_plain_end = std::chrono::high_resolution_clock::now();
                            time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                    t_plain_end - t_plain_begin).count();
                            result_count += local_result_count * pt.automorphism_size;
                            is_searched_ = true;
                        }
                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view();
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                t_lvm_destroy_end - t_lvm_destroy_begin).count();
                        if (result_count >= sm_.target_number) break;
                        if (g_exit) break;
                    }
                    const auto t_plain_total_end = std::chrono::high_resolution_clock::now();
                    time_phase1on_plain_total_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            t_plain_total_end - t_plain_total_begin).count();


                    for (auto &kv : iso_grouped_pre) {
                        auto &bucket = kv.second;
                        // TODO: 同构组复用搜索逻辑
                        const uint32_t group_sub_id = bucket.chosen_sub_id;
                        auto *subctx = get_isosubquery_context(group_sub_id);
                        if (subctx == nullptr) {
                            continue;
                        }

                        OrdersPerEdge* group_seed_orders =
                                std::get<1>(om_.get_automorphism_meta_special2(bucket.tasks[0].automorphism_id, group_sub_id));
                        if (group_seed_orders == nullptr) {
                            group_seed_orders = bucket.tasks[0].orders;
                        }
                        if (group_seed_orders == nullptr) {
                            continue;
                        }
                        const uint32_t u0 = group_seed_orders->indexing_order_[0];
                        const uint32_t u1 = group_seed_orders->indexing_order_[1];
                
                        const int32_t n0 = subctx->old2new[u0];
                        const int32_t n1 = subctx->old2new[u1];
                
                        const uint32_t ns = static_cast<uint32_t>(n0);
                        const uint32_t nt = static_cast<uint32_t>(n1);

                        const uint32_t sub_step0_key =
                                (u0 & 0xFFu) | ((u1 & 0xFFu) << 8) | ((ns & 0xFFu) << 16) | ((nt & 0xFFu) << 24);
                    
                        OrdersPerEdge *isosub_orders = subctx->om.get_orders({ns, nt});

                        if (isosub_orders == nullptr) {
                            continue;
                        }
                        uint32_t sub_k = subctx->sub_query->getVerticesCount();
                        bool sub_valid = subctx->lvm.create_candidates_only_for_subgraph_iso(
                                subctx->sub_query.get(),
                                query_graph,
                                *isosub_orders,
                                gvm_,
                                data_edge,
                                &subctx->old2new,
                                sub_k);
                        if (!sub_valid) {
                            continue;
                        }
                        LocalCandidatesSnapshot sub_candidates_snapshot =
                                subctx->lvm.extract_candidates_snapshot(subctx->sub_query.get());

                        LocalCandidatesSnapshot sub_candidates_snapshot_for_tasks = sub_candidates_snapshot;
                        subctx->lvm.restore_candidates_from(std::move(sub_candidates_snapshot));

                        std::vector<std::vector<uint32_t>> sub_embeddings_for_group;
                        bool sub_embeddings_ready = false;
                        static const std::vector<std::vector<uint32_t>> k_empty_mapped_new_by_old;
                        for (auto &pt : bucket.tasks) {
                            OrdersPerEdge* orders_for_group =
                                    std::get<1>(om_.get_automorphism_meta_special2(pt.automorphism_id, group_sub_id));
                            const uint32_t task_u0 = orders_for_group->indexing_order_[0];
                            const uint32_t task_u1 = orders_for_group->indexing_order_[1];

                            int32_t chosen_instance = -1;
                            for (size_t inst = 0; inst < subctx->mappings.size(); ++inst) {
                                bool has_u0 = false, has_u1 = false;
                                for (uint32_t old_u : subctx->mappings[inst]) {
                                    if (old_u == task_u0) has_u0 = true;
                                    if (old_u == task_u1) has_u1 = true;
                                }
                                if (has_u0 && has_u1) { chosen_instance = static_cast<int32_t>(inst); break; }
                            }
                            if (chosen_instance < 0) {
                                if (subctx->mappings.empty()) continue;
                                chosen_instance = 0;
                            }

                            std::vector<int32_t> local_new2old(sub_k, -1);
                            const auto& mapping = subctx->mappings[chosen_instance];
                            for (uint32_t c = 0; c < mapping.size(); ++c) {
                                local_new2old[c] = static_cast<int32_t>(mapping[c]);
                            }

                            uint32_t old_limit = static_cast<uint32_t>(subctx->old2new.size());
                            std::vector<std::vector<uint32_t>> mapped_for_task(old_limit);
                            for (uint32_t c = 0; c < mapping.size(); ++c) {
                                uint32_t old_u = mapping[c];
                                if (old_u < old_limit) {
                                    mapped_for_task[old_u].push_back(c);
                                }
                            }

                            bool task_pruned_valid = lvm_.create_candidates_only_with_sub_constraint(
                                    query_graph,
                                    *orders_for_group,
                                    gvm_,
                                    data_edge,
                                    sub_candidates_snapshot_for_tasks,
                                    mapped_for_task,
                                    subctx->sub_query.get(),
                                    &local_new2old);

                            if (!task_pruned_valid) {
                                continue;
                            }

                            if (!sub_embeddings_ready) {
                                subctx->lvm.build_view_from_candidates_for_subgraph(
                                        subctx->sub_query.get(), *isosub_orders, gvm_, &local_new2old);
                                sub_embeddings_for_group = subctx->sm.special_search_on_reduced_query(subctx->sub_query.get(),
                                                                                                           *isosub_orders,
                                                                                                           subctx->lvm,
                                                                                                           gvm_,
                                                                                                           subctx->sm.target_number,
                                                                                                           &local_new2old);
                                subctx->lvm.destroy_view();
                                sub_embeddings_ready = true;
                            }

                            if (sub_embeddings_for_group.empty()) {
                                continue;
                            }

                                
                            std::vector<std::pair<uint32_t, uint32_t>> forced_pairs;
                            forced_pairs.reserve(local_new2old.size());
                            for (size_t c = 0; c < local_new2old.size(); ++c) {
                                if (local_new2old[c] >= 0) {
                                    forced_pairs.emplace_back(static_cast<uint32_t>(local_new2old[c]), static_cast<uint32_t>(c));
                                }
                            }

                            lvm_.build_view_from_candidates(query_graph, *orders_for_group, gvm_);
                            for (const auto& emb : sub_embeddings_for_group) {
                                uint64_t local_result_count =
                                        sm_.constrained_search_with_subemb(query_graph,
                                                                            *orders_for_group,
                                                                            lvm_,
                                                                            gvm_,
                                                                            forced_pairs,
                                                                            emb);
                                result_count += local_result_count * pt.automorphism_size;
                                is_searched_ = true;

                                if (result_count >= sm_.target_number) break;
                                if (g_exit) break;
                            }
                            lvm_.destroy_view();
                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;


                        }



                    }

                    for (auto &kv : grouped_pre) {
                        auto &bucket = kv.second;
                        const uint32_t group_sub_id = bucket.chosen_sub_id;

                        auto *subctx = get_subquery_context(group_sub_id);
                        if (subctx == nullptr) {
                            continue;
                        }

                        OrdersPerEdge* group_seed_orders =
                                std::get<1>(om_.get_automorphism_meta_special(bucket.tasks[0].automorphism_id, group_sub_id));
                        if (group_seed_orders == nullptr) {
                            group_seed_orders = bucket.tasks[0].orders;
                        }
                        if (group_seed_orders == nullptr) {
                            continue;
                        }
                        const uint32_t u0 = group_seed_orders->indexing_order_[0];
                        const uint32_t u1 = group_seed_orders->indexing_order_[1];
                
                        const int32_t n0 = subctx->old2new[u0];
                        const int32_t n1 = subctx->old2new[u1];
                
                        const uint32_t ns = static_cast<uint32_t>(n0);
                        const uint32_t nt = static_cast<uint32_t>(n1);

                        const uint32_t sub_step0_key =
                                (u0 & 0xFFu) | ((u1 & 0xFFu) << 8) | ((ns & 0xFFu) << 16) | ((nt & 0xFFu) << 24);
                        const std::vector<std::vector<uint32_t>>* sub_mapped_new_orbit_peers = nullptr;
                        if (n0 >= 0 && n1 >= 0) {
                            auto it_orbit = subctx->mapped_new_orbit_peers_by_key.find(sub_step0_key);
                            if (it_orbit != subctx->mapped_new_orbit_peers_by_key.end())
                                sub_mapped_new_orbit_peers = &it_orbit->second;
                        }

                        OrdersPerEdge *sub_orders = subctx->om.get_orders({ns, nt});

                        if (sub_orders == nullptr) {
                            continue;
                        }
                        bool sub_valid = subctx->lvm.create_candidates_only_for_subgraph(subctx->sub_query.get(),
                                                                                               query_graph,
                                                                                               *sub_orders,
                                                                                               gvm_,
                                                                                               data_edge,
                                                                                               &subctx->new2old,
                                                                                               &subctx->sub_vertex_nlf,
                                                                                               sub_mapped_new_orbit_peers);

                        if (!sub_valid) {
                            continue;
                        }
                        phase1_before_sub_valid_check_count_ += 1;
                        LocalCandidatesSnapshot sub_candidates_snapshot =
                                subctx->lvm.extract_candidates_snapshot(subctx->sub_query.get());

                        LocalCandidatesSnapshot sub_candidates_snapshot_for_tasks = sub_candidates_snapshot;
                        subctx->lvm.restore_candidates_from(std::move(sub_candidates_snapshot));

                        std::vector<std::vector<uint32_t>> sub_embeddings_for_group;
                        bool sub_embeddings_ready = false;

                        static const std::vector<std::vector<uint32_t>> k_empty_mapped_new_by_old;

                        for (auto &pt : bucket.tasks) {
                            OrdersPerEdge* orders_for_group =
                                    std::get<1>(om_.get_automorphism_meta_special(pt.automorphism_id, group_sub_id));

                
                            const uint32_t task_u0 = orders_for_group->indexing_order_[0];
                            const uint32_t task_u1 = orders_for_group->indexing_order_[1];
                            const uint32_t task_ns = ns;
                            const uint32_t task_nt = nt;
                            const uint32_t key = (task_u0 & 0xFFu) | ((task_u1 & 0xFFu) << 8) |
                                                 ((task_ns & 0xFFu) << 16) | ((task_nt & 0xFFu) << 24);
                            auto it_mapped_new_set = subctx->mapped_new_set_by_key_and_old.find(key);

                            const std::vector<std::vector<uint32_t>> *mapped_for_task =
                                    (it_mapped_new_set != subctx->mapped_new_set_by_key_and_old.end())
                                            ? &it_mapped_new_set->second
                                            : &k_empty_mapped_new_by_old;
                          
                            bool task_pruned_valid = lvm_.create_candidates_only_with_sub_constraint(
                                    query_graph,
                                    *orders_for_group,
                                    gvm_,
                                    data_edge,
                                    sub_candidates_snapshot_for_tasks,
                                    *mapped_for_task,
                                    subctx->sub_query.get(),
                                    &subctx->new2old);


                            phase1_before_task_pruned_check_count_ += 1;
                            if (!task_pruned_valid) {
                                continue;
                            }

                            if (!sub_embeddings_ready) {
                                subctx->lvm.build_view_from_candidates_for_subgraph(
                                        subctx->sub_query.get(), *sub_orders, gvm_, &subctx->new2old);
                                sub_embeddings_for_group = subctx->sm.special_search_on_reduced_query(subctx->sub_query.get(),
                                                                                                           *sub_orders,
                                                                                                           subctx->lvm,
                                                                                                           gvm_,
                                                                                                           subctx->sm.target_number,
                                                                                                           &subctx->new2old);
                                subctx->lvm.destroy_view();
                                sub_embeddings_ready = true;
                                
                            }

                            phase1_before_sub_embeddings_check_count_ += 1;
                            if (sub_embeddings_for_group.empty()) {
                                continue;
                            }
                            phase1_after_sub_embeddings_check_count_ += 1;

                                
                            auto it_mapping_ids = subctx->mapping_ids_index.find(key);
                            const std::vector<uint32_t>& candidate_mapping_ids = it_mapping_ids->second;

                            lvm_.build_view_from_candidates(query_graph, *orders_for_group, gvm_);
                            for (const auto& emb : sub_embeddings_for_group) {
                                for (uint32_t selected_mapping_id : candidate_mapping_ids) {
                                    auto it_tpl = subctx->forced_pairs_template.find(selected_mapping_id);

                                    uint64_t local_result_count =
                                            sm_.constrained_search_with_subemb(query_graph,
                                                                                *orders_for_group,
                                                                                lvm_,
                                                                                gvm_,
                                                                                it_tpl->second,
                                                                                emb);
                                    result_count += local_result_count * pt.automorphism_size;
                                    is_searched_ = true;

                                    if (result_count >= sm_.target_number) break;
                                    if (g_exit) break;
                                }
                                if (result_count >= sm_.target_number) break;
                                if (g_exit) break;
                            }
                            lvm_.destroy_view();
                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;

                        }
                    }
                }
            }
        }
    }
    else {
        label_triple = label_triple_initial; 
        data_edge = update.edge_;
        std::vector<Edge> *mapped_query_edges_del = om_.get_mapped_edges(label_triple);
        const size_t phase1_matching_query_edge_count_del = mapped_query_edges_del ? mapped_query_edges_del->size() : 0;
        auto mapped_automorphism = om_.get_mapped_automorphism(label_triple); 
        if (mapped_automorphism != nullptr) { 
            is_relevant_ = true; 
            if (enable_local_view) { 
                if (phase1_matching_query_edge_count_del <= 1) { 
                    for (auto automorphism_id: *mapped_automorphism) { 
                        auto meta = om_.get_automorphism_meta(automorphism_id); 
                        OrdersPerEdge *orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta); 
                        bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge); 

                        if (is_valid) { 
                            search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_; 
                            search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_; 
                        } else { 
                            non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_; 
                            non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_; 

                            if (lvm_.generate_visited_neighbor_count_ == 0) 
                                direct_rejection_count_ += 1; 
                        }
                        first_indexing_vertex_ += lvm_.first_vertex_neighbor_; 
                        if (enable_search && is_valid) {
                            uint64_t local_result_count = sm_.search_on_reduced_query(query_graph, *orders, lvm_, gvm_);
                            result_count += local_result_count * automorphism_size; 
                            is_searched_ = true; 
                        }

                        lvm_.destroy_view(); 

                        if (result_count >= sm_.target_number) break; 
                        if (g_exit) break; 
                    }
                } else {
                     struct PreGroupTask {
                        OrdersPerEdge* orders = nullptr;
                        uint32_t automorphism_size = 1;
                        uint32_t automorphism_id = 0;
                    };
                    struct PreGroupBucket {
                        uint32_t group_id = 0;
                        std::vector<PreGroupTask> tasks;
                        uint32_t chosen_sub_id = std::numeric_limits<uint32_t>::max();
                    };

                    std::vector<PreGroupTask> plain_search_pre;
                    std::unordered_map<uint32_t, PreGroupBucket> grouped_pre;
                    std::unordered_map<uint32_t, PreGroupBucket> iso_grouped_pre;
                    plain_search_pre.reserve(mapped_automorphism->size());
                    grouped_pre.reserve(32);
                    iso_grouped_pre.reserve(32);

                    for (auto automorphism_id : *mapped_automorphism) {
                        auto meta = om_.get_automorphism_meta(automorphism_id);
                        OrdersPerEdge* orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta);
                        if (orders == nullptr) continue;

                        const uint32_t u0 = orders->indexing_order_[0];
                        const uint32_t u1 = orders->indexing_order_[1];

                        if (!gvm_.nlf_check(u0, data_edge.first) || !gvm_.nlf_check(u1, data_edge.second)) {
                            continue;
                        }

                        DirectedEdgeKey key{u0, u1, 0};
                        uint32_t gid = 0;
                        auto it_gid = phase1_edge_group_id_.find(key);
                        if (it_gid != phase1_edge_group_id_.end()) {
                             gid = it_gid->second;
                        }

                        if (gid == 0) {
                            uint32_t iso_gid = 0;
                            auto it_iso = phase2_edge_group_id_.find(key);
                            if (it_iso != phase2_edge_group_id_.end()) {
                                iso_gid = it_iso->second;
                            }
                            if (iso_gid >= 1) {
                                auto& ib = iso_grouped_pre[iso_gid];
                                ib.group_id = iso_gid;
                                ib.tasks.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                                ib.chosen_sub_id = phase2_group_chosen_sub_id_[iso_gid];
                                continue;
                            }
                            plain_search_pre.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                            continue;
                        }

                        auto& b = grouped_pre[gid];
                        b.group_id = gid;
                        b.tasks.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                        b.chosen_sub_id = phase1_group_chosen_sub_id_[gid]; 
                    }

                    for (auto it = grouped_pre.begin(); it != grouped_pre.end(); ) {
                        if (it->second.tasks.size() <= 1) {
                            for (auto& t : it->second.tasks)
                                plain_search_pre.push_back(t);
                            it = grouped_pre.erase(it);
                        } else {
                            ++it;
                        }
                    }
                    for (auto it = iso_grouped_pre.begin(); it != iso_grouped_pre.end(); ) {
                        if (it->second.tasks.size() <= 1) {
                            for (auto& t : it->second.tasks)
                                plain_search_pre.push_back(t);
                            it = iso_grouped_pre.erase(it);
                        } else {
                            ++it;
                        }
                    }

                    const auto t_plain_total_begin = std::chrono::high_resolution_clock::now();
                    for (const auto &pt : plain_search_pre) {
                        if (pt.orders == nullptr) continue;
                        // 这里走“单边组”回溯：直接按当前 orders 构建局部视图并搜索（不依赖 pending candidates）。
                        auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                        bool is_valid = lvm_.create_view(query_graph, *pt.orders, gvm_, data_edge);
                        auto t_lvm_end = std::chrono::high_resolution_clock::now();
                        time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                t_lvm_end - t_lvm_begin).count();
                        if (enable_search && is_valid) {
                            auto t_plain_begin = std::chrono::high_resolution_clock::now();
                            uint64_t local_result_count =
                                    sm_.search_on_reduced_query(query_graph, *pt.orders, lvm_, gvm_);
                            auto t_plain_end = std::chrono::high_resolution_clock::now();
                            time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                    t_plain_end - t_plain_begin).count();
                            result_count += local_result_count * pt.automorphism_size;
                            is_searched_ = true;
                        }
                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view();
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                t_lvm_destroy_end - t_lvm_destroy_begin).count();
                        if (result_count >= sm_.target_number) break;
                        if (g_exit) break;
                    }
                    const auto t_plain_total_end = std::chrono::high_resolution_clock::now();
                    time_phase1on_plain_total_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            t_plain_total_end - t_plain_total_begin).count();


                    for (auto &kv : iso_grouped_pre) {
                        auto &bucket = kv.second;
                        // TODO: 同构组复用搜索逻辑
                        const uint32_t group_sub_id = bucket.chosen_sub_id;
                        auto *subctx = get_isosubquery_context(group_sub_id);
                        if (subctx == nullptr) {
                            continue;
                        }

                        OrdersPerEdge* group_seed_orders =
                                std::get<1>(om_.get_automorphism_meta_special2(bucket.tasks[0].automorphism_id, group_sub_id));
                        if (group_seed_orders == nullptr) {
                            group_seed_orders = bucket.tasks[0].orders;
                        }
                        if (group_seed_orders == nullptr) {
                            continue;
                        }
                        const uint32_t u0 = group_seed_orders->indexing_order_[0];
                        const uint32_t u1 = group_seed_orders->indexing_order_[1];
                
                        const int32_t n0 = subctx->old2new[u0];
                        const int32_t n1 = subctx->old2new[u1];
                
                        const uint32_t ns = static_cast<uint32_t>(n0);
                        const uint32_t nt = static_cast<uint32_t>(n1);

                        const uint32_t sub_step0_key =
                                (u0 & 0xFFu) | ((u1 & 0xFFu) << 8) | ((ns & 0xFFu) << 16) | ((nt & 0xFFu) << 24);
                    
                        OrdersPerEdge *isosub_orders = subctx->om.get_orders({ns, nt});

                        if (isosub_orders == nullptr) {
                            continue;
                        }
                        uint32_t sub_k = subctx->sub_query->getVerticesCount();
                        bool sub_valid = subctx->lvm.create_candidates_only_for_subgraph_iso(
                                subctx->sub_query.get(),
                                query_graph,
                                *isosub_orders,
                                gvm_,
                                data_edge,
                                &subctx->old2new,
                                sub_k);
                        if (!sub_valid) {
                            continue;
                        }
                        LocalCandidatesSnapshot sub_candidates_snapshot =
                                subctx->lvm.extract_candidates_snapshot(subctx->sub_query.get());

                        LocalCandidatesSnapshot sub_candidates_snapshot_for_tasks = sub_candidates_snapshot;
                        subctx->lvm.restore_candidates_from(std::move(sub_candidates_snapshot));

                        std::vector<std::vector<uint32_t>> sub_embeddings_for_group;
                        bool sub_embeddings_ready = false;
                        static const std::vector<std::vector<uint32_t>> k_empty_mapped_new_by_old;
                        for (auto &pt : bucket.tasks) {
                            OrdersPerEdge* orders_for_group =
                                    std::get<1>(om_.get_automorphism_meta_special2(pt.automorphism_id, group_sub_id));
                            const uint32_t task_u0 = orders_for_group->indexing_order_[0];
                            const uint32_t task_u1 = orders_for_group->indexing_order_[1];

                            int32_t chosen_instance = -1;
                            for (size_t inst = 0; inst < subctx->mappings.size(); ++inst) {
                                bool has_u0 = false, has_u1 = false;
                                for (uint32_t old_u : subctx->mappings[inst]) {
                                    if (old_u == task_u0) has_u0 = true;
                                    if (old_u == task_u1) has_u1 = true;
                                }
                                if (has_u0 && has_u1) { chosen_instance = static_cast<int32_t>(inst); break; }
                            }
                            if (chosen_instance < 0) {
                                if (subctx->mappings.empty()) continue;
                                chosen_instance = 0;
                            }

                            std::vector<int32_t> local_new2old(sub_k, -1);
                            const auto& mapping = subctx->mappings[chosen_instance];
                            for (uint32_t c = 0; c < mapping.size(); ++c) {
                                local_new2old[c] = static_cast<int32_t>(mapping[c]);
                            }

                            uint32_t old_limit = static_cast<uint32_t>(subctx->old2new.size());
                            std::vector<std::vector<uint32_t>> mapped_for_task(old_limit);
                            for (uint32_t c = 0; c < mapping.size(); ++c) {
                                uint32_t old_u = mapping[c];
                                if (old_u < old_limit) {
                                    mapped_for_task[old_u].push_back(c);
                                }
                            }

                            bool task_pruned_valid = lvm_.create_candidates_only_with_sub_constraint(
                                    query_graph,
                                    *orders_for_group,
                                    gvm_,
                                    data_edge,
                                    sub_candidates_snapshot_for_tasks,
                                    mapped_for_task,
                                    subctx->sub_query.get(),
                                    &local_new2old);

                            if (!task_pruned_valid) {
                                continue;
                            }

                            if (!sub_embeddings_ready) {
                                subctx->lvm.build_view_from_candidates_for_subgraph(
                                        subctx->sub_query.get(), *isosub_orders, gvm_, &local_new2old);
                                sub_embeddings_for_group = subctx->sm.special_search_on_reduced_query(subctx->sub_query.get(),
                                                                                                           *isosub_orders,
                                                                                                           subctx->lvm,
                                                                                                           gvm_,
                                                                                                           subctx->sm.target_number,
                                                                                                           &local_new2old);
                                subctx->lvm.destroy_view();
                                sub_embeddings_ready = true;
                            }

                            if (sub_embeddings_for_group.empty()) {
                                continue;
                            }

                                
                            std::vector<std::pair<uint32_t, uint32_t>> forced_pairs;
                            forced_pairs.reserve(local_new2old.size());
                            for (size_t c = 0; c < local_new2old.size(); ++c) {
                                if (local_new2old[c] >= 0) {
                                    forced_pairs.emplace_back(static_cast<uint32_t>(local_new2old[c]), static_cast<uint32_t>(c));
                                }
                            }

                            lvm_.build_view_from_candidates(query_graph, *orders_for_group, gvm_);
                            for (const auto& emb : sub_embeddings_for_group) {
                                uint64_t local_result_count =
                                        sm_.constrained_search_with_subemb(query_graph,
                                                                            *orders_for_group,
                                                                            lvm_,
                                                                            gvm_,
                                                                            forced_pairs,
                                                                            emb);
                                result_count += local_result_count * pt.automorphism_size;
                                is_searched_ = true;

                                if (result_count >= sm_.target_number) break;
                                if (g_exit) break;
                            }
                            lvm_.destroy_view();
                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;


                        }



                    }

                    for (auto &kv : grouped_pre) {
                        auto &bucket = kv.second;
                        const uint32_t group_sub_id = bucket.chosen_sub_id;

                        auto *subctx = get_subquery_context(group_sub_id);
                        if (subctx == nullptr) {
                            continue;
                        }

                        OrdersPerEdge* group_seed_orders =
                                std::get<1>(om_.get_automorphism_meta_special(bucket.tasks[0].automorphism_id, group_sub_id));
                        if (group_seed_orders == nullptr) {
                            group_seed_orders = bucket.tasks[0].orders;
                        }
                        if (group_seed_orders == nullptr) {
                            continue;
                        }
                        const uint32_t u0 = group_seed_orders->indexing_order_[0];
                        const uint32_t u1 = group_seed_orders->indexing_order_[1];
                
                        const int32_t n0 = subctx->old2new[u0];
                        const int32_t n1 = subctx->old2new[u1];
                
                        const uint32_t ns = static_cast<uint32_t>(n0);
                        const uint32_t nt = static_cast<uint32_t>(n1);

                        const uint32_t sub_step0_key =
                                (u0 & 0xFFu) | ((u1 & 0xFFu) << 8) | ((ns & 0xFFu) << 16) | ((nt & 0xFFu) << 24);
                        const std::vector<std::vector<uint32_t>>* sub_mapped_new_orbit_peers = nullptr;
                        if (n0 >= 0 && n1 >= 0) {
                            auto it_orbit = subctx->mapped_new_orbit_peers_by_key.find(sub_step0_key);
                            if (it_orbit != subctx->mapped_new_orbit_peers_by_key.end())
                                sub_mapped_new_orbit_peers = &it_orbit->second;
                        }

                        OrdersPerEdge *sub_orders = subctx->om.get_orders({ns, nt});

                        if (sub_orders == nullptr) {
                            continue;
                        }
                        bool sub_valid = subctx->lvm.create_candidates_only_for_subgraph(subctx->sub_query.get(),
                                                                                               query_graph,
                                                                                               *sub_orders,
                                                                                               gvm_,
                                                                                               data_edge,
                                                                                               &subctx->new2old,
                                                                                               &subctx->sub_vertex_nlf,
                                                                                               sub_mapped_new_orbit_peers);

                        if (!sub_valid) {
                            continue;
                        }
                        phase1_before_sub_valid_check_count_ += 1;
                        LocalCandidatesSnapshot sub_candidates_snapshot =
                                subctx->lvm.extract_candidates_snapshot(subctx->sub_query.get());

                        LocalCandidatesSnapshot sub_candidates_snapshot_for_tasks = sub_candidates_snapshot;
                        subctx->lvm.restore_candidates_from(std::move(sub_candidates_snapshot));

                        std::vector<std::vector<uint32_t>> sub_embeddings_for_group;
                        bool sub_embeddings_ready = false;

                        static const std::vector<std::vector<uint32_t>> k_empty_mapped_new_by_old;

                        for (auto &pt : bucket.tasks) {
                            OrdersPerEdge* orders_for_group =
                                    std::get<1>(om_.get_automorphism_meta_special(pt.automorphism_id, group_sub_id));

                
                            const uint32_t task_u0 = orders_for_group->indexing_order_[0];
                            const uint32_t task_u1 = orders_for_group->indexing_order_[1];
                            const uint32_t task_ns = ns;
                            const uint32_t task_nt = nt;
                            const uint32_t key = (task_u0 & 0xFFu) | ((task_u1 & 0xFFu) << 8) |
                                                 ((task_ns & 0xFFu) << 16) | ((task_nt & 0xFFu) << 24);
                            auto it_mapped_new_set = subctx->mapped_new_set_by_key_and_old.find(key);

                            const std::vector<std::vector<uint32_t>> *mapped_for_task =
                                    (it_mapped_new_set != subctx->mapped_new_set_by_key_and_old.end())
                                            ? &it_mapped_new_set->second
                                            : &k_empty_mapped_new_by_old;
                          
                            bool task_pruned_valid = lvm_.create_candidates_only_with_sub_constraint(
                                    query_graph,
                                    *orders_for_group,
                                    gvm_,
                                    data_edge,
                                    sub_candidates_snapshot_for_tasks,
                                    *mapped_for_task,
                                    subctx->sub_query.get(),
                                    &subctx->new2old);


                            phase1_before_task_pruned_check_count_ += 1;
                            if (!task_pruned_valid) {
                                continue;
                            }

                            if (!sub_embeddings_ready) {
                                subctx->lvm.build_view_from_candidates_for_subgraph(
                                        subctx->sub_query.get(), *sub_orders, gvm_, &subctx->new2old);
                                sub_embeddings_for_group = subctx->sm.special_search_on_reduced_query(subctx->sub_query.get(),
                                                                                                           *sub_orders,
                                                                                                           subctx->lvm,
                                                                                                           gvm_,
                                                                                                           subctx->sm.target_number,
                                                                                                           &subctx->new2old);
                                subctx->lvm.destroy_view();
                                sub_embeddings_ready = true;
                                
                            }

                            phase1_before_sub_embeddings_check_count_ += 1;
                            if (sub_embeddings_for_group.empty()) {
                                continue;
                            }
                            phase1_after_sub_embeddings_check_count_ += 1;

                                
                            auto it_mapping_ids = subctx->mapping_ids_index.find(key);
                            const std::vector<uint32_t>& candidate_mapping_ids = it_mapping_ids->second;

                            lvm_.build_view_from_candidates(query_graph, *orders_for_group, gvm_);
                            for (const auto& emb : sub_embeddings_for_group) {
                                for (uint32_t selected_mapping_id : candidate_mapping_ids) {
                                    auto it_tpl = subctx->forced_pairs_template.find(selected_mapping_id);

                                    uint64_t local_result_count =
                                            sm_.constrained_search_with_subemb(query_graph,
                                                                                *orders_for_group,
                                                                                lvm_,
                                                                                gvm_,
                                                                                it_tpl->second,
                                                                                emb);
                                    result_count += local_result_count * pt.automorphism_size;
                                    is_searched_ = true;

                                    if (result_count >= sm_.target_number) break;
                                    if (g_exit) break;
                                }
                                if (result_count >= sm_.target_number) break;
                                if (g_exit) break;
                            }
                            lvm_.destroy_view();
                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;

                        }
                    }

                }
            }
        }

        auto it = gvm_.get_mapped_views(label_triple); // 局部处理完成后，再更新 global view（删除边上的匹配结构）
        if (it != nullptr) { // 查询中存在该标签模式对应的 global view
            is_relevant_ = true; // 本更新与全局视图维护相关
            relevant_update_count_ += 1; // 相关更新统计 +1
            // [无向边核对] 删除边同样在两种有向朝向下各更新一次 global view，再恢复 data_edge/label_triple。
            for (uint32_t i = 0; i < 2; ++i) { // 同样两种有向朝向各删一次
                it = gvm_.get_mapped_views(label_triple); // 按当前朝向重新取映射

                uint32_t view_id = it->second; // 当前朝向下 global view 编号
                gvm_.update_view(update.op_, view_id, data_edge, label_triple); // op 为 '-'，从对应 global view 中移除该边

                std::swap(data_edge.first, data_edge.second); // 换向以更新反向 global view
                std::swap(label_triple.src_label_, label_triple.dst_label_); // 同步交换端点标签
            }

            std::swap(data_edge.first, data_edge.second); // 恢复 data_edge 为 update 原始顺序
            std::swap(label_triple.src_label_, label_triple.dst_label_); // 恢复 label_triple
        }

        if (is_relevant_) { // 若曾命中 global view，则 NLF 也需随删除更新
            gvm_.update_nlf_view(update.op_, data_edge, label_triple); // 删除边：更新 NLF
        }
    }
    edge_process_count_ += 1; 
    if (is_relevant_) { 
        if (result_count > 0) { 
            positive_count_ += 1; 
        }
        if (is_searched_) {
            search_count_ += 1; 
        }

        result_count_ += result_count; 
        invalid_partial_result_count_ += sm_.invalid_partial_result_count_; 
        partial_result_count_ += sm_.partial_result_count_; 
        iso_conflict_count_ += sm_.iso_conflict_count_; 
        si_empty_count_ += sm_.si_empty_count_; 
        lc_empty_count_ += sm_.lc_empty_count_; 

        sm_.reset_performance_counters(); 
    }
    return result_count; 
}

uint64_t StreamingEngine::execute_phase4_on(const Graph *query_graph, const Update &update, bool enable_local_view,
                                           bool enable_search) {
   

    LabelTriple label_triple = update.labels_; 
    Edge data_edge = update.edge_; 
    const LabelTriple label_triple_initial = label_triple; 
    uint64_t result_count = 0; 
    is_relevant_ = false; 
    is_searched_ = false;

    if (update.op_ == '+') { 
      
        auto it = gvm_.get_mapped_views(label_triple); 
        if (it != nullptr) {
            is_relevant_ = true; 
            relevant_update_count_ += 1; 
            auto t_gvm_begin = std::chrono::high_resolution_clock::now();
            for (uint32_t i = 0; i < 2; ++i) { 
                it = gvm_.get_mapped_views(label_triple); 
                uint32_t view_id = it->second; 
                gvm_.update_view(update.op_, view_id, data_edge, label_triple); 
                std::swap(data_edge.first, data_edge.second); 
                std::swap(label_triple.src_label_, label_triple.dst_label_); 
            }
            std::swap(data_edge.first, data_edge.second); 
            std::swap(label_triple.src_label_, label_triple.dst_label_); 
            auto t_gvm_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_views_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_gvm_end - t_gvm_begin).count();
        }
        if (is_relevant_) { 
            auto t_nlf_begin = std::chrono::high_resolution_clock::now();
            gvm_.update_nlf_view(update.op_, data_edge, label_triple); 
            auto t_nlf_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_nlf_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_nlf_end - t_nlf_begin).count();
        }     
        std::vector<Edge> *mapped_query_edges = om_.get_mapped_edges(label_triple);
        const size_t phase1_matching_query_edge_count = mapped_query_edges ? mapped_query_edges->size() : 0;
        auto mapped_automorphism = om_.get_mapped_automorphism(label_triple); 
        if (mapped_automorphism != nullptr) { 
            is_relevant_ = true; 
            if (enable_local_view) { 
                if (phase1_matching_query_edge_count <= 1) { 
                
                    for (auto automorphism_id: *mapped_automorphism) { 
                        auto meta = om_.get_automorphism_meta(automorphism_id); 
                        OrdersPerEdge *orders = std::get<1>(meta); 
                        uint32_t automorphism_size = std::get<2>(meta); 
                        auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                        bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge); 
                        auto t_lvm_end = std::chrono::high_resolution_clock::now();
                        time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_end - t_lvm_begin).count();

                        if (is_valid) { 
                            search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_; 
                            search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_; 
                        } else { // 视图无效或提前剪枝
                            non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                            non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;

                            if (lvm_.generate_visited_neighbor_count_ == 0) 
                                direct_rejection_count_ += 1; 
                        }
                        first_indexing_vertex_ += lvm_.first_vertex_neighbor_; 
                        if (enable_search && is_valid) { 
                            auto t_plain_begin = std::chrono::high_resolution_clock::now();
                            uint64_t local_result_count = sm_.search_on_reduced_query(query_graph, *orders, lvm_, gvm_); 
                            auto t_plain_end = std::chrono::high_resolution_clock::now();
                            time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_plain_end - t_plain_begin).count();
                            result_count += local_result_count * automorphism_size;
                            is_searched_ = true; 
                        }

                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view(); 
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_destroy_end - t_lvm_destroy_begin).count();

                        if (result_count >= sm_.target_number) break; 
                        if (g_exit) break;
                    }
                }else {
                    struct PreGroupTask {
                        OrdersPerEdge* orders = nullptr;
                        uint32_t automorphism_size = 1;
                        uint32_t automorphism_id = 0;
                    };
                    struct PreGroupBucket {
                        uint32_t group_id = 0;
                        std::vector<PreGroupTask> tasks;
                        uint32_t chosen_sub_id = std::numeric_limits<uint32_t>::max();
                    };
                    std::vector<PreGroupTask> plain_search_pre;
                    std::unordered_map<uint32_t, PreGroupBucket> grouped_pre;
                    plain_search_pre.reserve(mapped_automorphism->size());
                    grouped_pre.reserve(32);

                    for (auto automorphism_id : *mapped_automorphism) {
                        auto meta = om_.get_automorphism_meta(automorphism_id);
                        OrdersPerEdge* orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta);
                        if (orders == nullptr) continue;

                        const uint32_t u0 = orders->indexing_order_[0];
                        const uint32_t u1 = orders->indexing_order_[1];

                        if (!gvm_.nlf_check(u0, data_edge.first) || !gvm_.nlf_check(u1, data_edge.second)) {
                            continue;
                        }

                        DirectedEdgeKey key{u0, u1, 0};
                        uint32_t gid = 0;
                        auto it_gid = phase1_edge_group_id_.find(key);
                        if (it_gid != phase1_edge_group_id_.end()) {
                             gid = it_gid->second;
                        }

                        if (gid == 0) {
                            plain_search_pre.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                            continue;
                        }

                        auto& b = grouped_pre[gid];
                        b.group_id = gid;
                        b.tasks.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                        b.chosen_sub_id = phase1_group_chosen_sub_id_[gid]; 
                    }
                    const auto t_plain_total_begin = std::chrono::high_resolution_clock::now();
                    if (plain_search_pre.size() <= 1){
                        for (const auto &pt : plain_search_pre) {
                            if (pt.orders == nullptr) continue;
                            auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                            bool is_valid = lvm_.create_view(query_graph, *pt.orders, gvm_, data_edge);
                            auto t_lvm_end = std::chrono::high_resolution_clock::now();
                            time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                    t_lvm_end - t_lvm_begin).count();

                            if (enable_search && is_valid) {
                                auto t_plain_begin = std::chrono::high_resolution_clock::now();
                                uint64_t local_result_count =
                                        sm_.search_on_reduced_query(query_graph, *pt.orders, lvm_, gvm_);
                                auto t_plain_end = std::chrono::high_resolution_clock::now();
                                time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                        t_plain_end - t_plain_begin).count();
                                result_count += local_result_count * pt.automorphism_size;
                                is_searched_ = true;
                            }

                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view();
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                t_lvm_destroy_end - t_lvm_destroy_begin).count();

                        if (result_count >= sm_.target_number) break;
                        if (g_exit) break;
                        }
                    }else{
                        sm_.clear_phase3_block_cache();
                        for (const auto &pt : plain_search_pre) {
                            if (pt.orders == nullptr) continue;
                            auto meta = om_.get_block_reuse_automorphism_meta(pt.automorphism_id);
                            OrdersPerEdge *orders = std::get<1>(meta);
                            uint32_t automorphism_size = std::get<2>(meta);
                            auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                            bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge);
                            auto t_lvm_end = std::chrono::high_resolution_clock::now();
                            time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_end - t_lvm_begin).count();
                            if (is_valid) {
                                search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                                search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                            } else {
                                non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                                non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;

                                if (lvm_.generate_visited_neighbor_count_ == 0)
                                    direct_rejection_count_ += 1;
                            }
                            first_indexing_vertex_ += lvm_.first_vertex_neighbor_;
                            if (enable_search && is_valid) {
                                phase3_insert_block_reuse_calls_ += 1;
                                auto t_plain_begin = std::chrono::high_resolution_clock::now();
                                uint64_t local_result_count =
                                        sm_.block_reuse_search_on_reduced_query(query_graph, *orders, lvm_, gvm_, query_bctree_);
                                auto t_plain_end = std::chrono::high_resolution_clock::now();
                                time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_plain_end - t_plain_begin).count();

                                result_count += local_result_count * automorphism_size;
                                is_searched_ = true;
                            }

                            auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                            lvm_.destroy_view();
                            auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                            time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_destroy_end - t_lvm_destroy_begin).count();

                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;
                        }
                        sm_.clear_phase3_block_cache();
                    }
                    const auto t_plain_total_end = std::chrono::high_resolution_clock::now();
                    time_phase1on_plain_total_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            t_plain_total_end - t_plain_total_begin).count();
                      
                    for (auto &kv : grouped_pre) {
                        auto &bucket = kv.second;
                        const uint32_t group_sub_id = bucket.chosen_sub_id;

                        auto *subctx = get_subquery_context(group_sub_id);
                        if (subctx == nullptr) {
                            continue;
                        }

                        OrdersPerEdge* group_seed_orders =
                                std::get<1>(om_.get_automorphism_meta_special(bucket.tasks[0].automorphism_id, group_sub_id));
                        if (group_seed_orders == nullptr) {
                            group_seed_orders = bucket.tasks[0].orders;
                        }
                        if (group_seed_orders == nullptr) {
                            continue;
                        }
                        const uint32_t u0 = group_seed_orders->indexing_order_[0];
                        const uint32_t u1 = group_seed_orders->indexing_order_[1];
                
                        const int32_t n0 = subctx->old2new[u0];
                        const int32_t n1 = subctx->old2new[u1];
                
                        const uint32_t ns = static_cast<uint32_t>(n0);
                        const uint32_t nt = static_cast<uint32_t>(n1);

                        const uint32_t sub_step0_key =
                                (u0 & 0xFFu) | ((u1 & 0xFFu) << 8) | ((ns & 0xFFu) << 16) | ((nt & 0xFFu) << 24);
                        const std::vector<std::vector<uint32_t>>* sub_mapped_new_orbit_peers = nullptr;
                        if (n0 >= 0 && n1 >= 0) {
                            auto it_orbit = subctx->mapped_new_orbit_peers_by_key.find(sub_step0_key);
                            if (it_orbit != subctx->mapped_new_orbit_peers_by_key.end())
                                sub_mapped_new_orbit_peers = &it_orbit->second;
                        }

                        OrdersPerEdge *sub_orders = subctx->om.get_orders({ns, nt});

                        if (sub_orders == nullptr) {
                            continue;
                        }
                        bool sub_valid = subctx->lvm.create_candidates_only_for_subgraph(subctx->sub_query.get(),
                                                                                               query_graph,
                                                                                               *sub_orders,
                                                                                               gvm_,
                                                                                               data_edge,
                                                                                               &subctx->new2old,
                                                                                               &subctx->sub_vertex_nlf,
                                                                                               sub_mapped_new_orbit_peers);

                        if (!sub_valid) {
                            continue;
                        }
                        phase1_before_sub_valid_check_count_ += 1;
                        LocalCandidatesSnapshot sub_candidates_snapshot =
                                subctx->lvm.extract_candidates_snapshot(subctx->sub_query.get());

                        LocalCandidatesSnapshot sub_candidates_snapshot_for_tasks = sub_candidates_snapshot;
                        subctx->lvm.restore_candidates_from(std::move(sub_candidates_snapshot));

                        std::vector<std::vector<uint32_t>> sub_embeddings_for_group;
                        bool sub_embeddings_ready = false;

                        static const std::vector<std::vector<uint32_t>> k_empty_mapped_new_by_old;

                        for (auto &pt : bucket.tasks) {
                            OrdersPerEdge* orders_for_group =
                                    std::get<1>(om_.get_automorphism_meta_special(pt.automorphism_id, group_sub_id));
                            const uint32_t task_u0 = orders_for_group->indexing_order_[0];
                            const uint32_t task_u1 = orders_for_group->indexing_order_[1];
                            const uint32_t task_ns = ns;
                            const uint32_t task_nt = nt;
                            const uint32_t key = (task_u0 & 0xFFu) | ((task_u1 & 0xFFu) << 8) |
                                                 ((task_ns & 0xFFu) << 16) | ((task_nt & 0xFFu) << 24);
                            auto it_mapped_new_set = subctx->mapped_new_set_by_key_and_old.find(key);

                            const std::vector<std::vector<uint32_t>> *mapped_for_task =
                                    (it_mapped_new_set != subctx->mapped_new_set_by_key_and_old.end())
                                            ? &it_mapped_new_set->second
                                            : &k_empty_mapped_new_by_old;
                          
                            bool task_pruned_valid = lvm_.create_candidates_only_with_sub_constraint(
                                    query_graph,
                                    *orders_for_group,
                                    gvm_,
                                    data_edge,
                                    sub_candidates_snapshot_for_tasks,
                                    *mapped_for_task,
                                    subctx->sub_query.get(),
                                    &subctx->new2old);


                            phase1_before_task_pruned_check_count_ += 1;
                            if (!task_pruned_valid) {
                                continue;
                            }

                            if (!sub_embeddings_ready) {
                                subctx->lvm.build_view_from_candidates_for_subgraph(
                                        subctx->sub_query.get(), *sub_orders, gvm_, &subctx->new2old);
                                sub_embeddings_for_group = subctx->sm.special_search_on_reduced_query(subctx->sub_query.get(),
                                                                                                           *sub_orders,
                                                                                                           subctx->lvm,
                                                                                                           gvm_,
                                                                                                           subctx->sm.target_number,
                                                                                                           &subctx->new2old);
                                subctx->lvm.destroy_view();
                                sub_embeddings_ready = true;
                                
                            }

                            phase1_before_sub_embeddings_check_count_ += 1;
                            if (sub_embeddings_for_group.empty()) {
                                continue;
                            }
                            phase1_after_sub_embeddings_check_count_ += 1;

                                
                            auto it_mapping_ids = subctx->mapping_ids_index.find(key);
                            const std::vector<uint32_t>& candidate_mapping_ids = it_mapping_ids->second;

                            lvm_.build_view_from_candidates(query_graph, *orders_for_group, gvm_);
                            for (const auto& emb : sub_embeddings_for_group) {
                                for (uint32_t selected_mapping_id : candidate_mapping_ids) {
                                    auto it_tpl = subctx->forced_pairs_template.find(selected_mapping_id);
                                    uint64_t local_result_count =
                                            sm_.constrained_search_with_subemb(query_graph,
                                                                                *orders_for_group,
                                                                                lvm_,
                                                                                gvm_,
                                                                                it_tpl->second,
                                                                                emb);
                                    result_count += local_result_count * pt.automorphism_size;
                                    is_searched_ = true;
                                    if (result_count >= sm_.target_number) break;
                                    if (g_exit) break;
                                }
                                if (result_count >= sm_.target_number) break;
                                if (g_exit) break;
                            }
                            lvm_.destroy_view();
                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;

                        }
                    }
                }
            }
        }
    }
    else {

        label_triple = label_triple_initial; 
        data_edge = update.edge_;
        std::vector<Edge> *mapped_query_edges_del = om_.get_mapped_edges(label_triple);
        const size_t phase1_matching_query_edge_count_del = mapped_query_edges_del ? mapped_query_edges_del->size() : 0;

        auto mapped_automorphism = om_.get_mapped_automorphism(label_triple); // 与插入分支相同：按标签三元组取订单桶
        if (mapped_automorphism != nullptr) { // 存在匹配顺序时才可能建局部视图
            is_relevant_ = true; // 标记与查询订单相关
            if (enable_local_view) { // 与插入分支相同：由调用方决定是否物化 LocalView

                if (phase1_matching_query_edge_count_del <= 1) { // 单条匹配查询边：原逻辑
                    for (auto automorphism_id: *mapped_automorphism) { // 逐个 automorphism 桶处理
                        auto meta = om_.get_automorphism_meta(automorphism_id); // 元组：订单指针、轨道大小等
                        OrdersPerEdge *orders = std::get<1>(meta); // 当前桶的匹配/索引顺序
                        uint32_t automorphism_size = std::get<2>(meta); // 轨道大小，用于结果计数放大
                        bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge); // 注意：此时全局视图尚未删边，与 RapidFlow 删除语义一致

                        if (is_valid) { // 删除分支下局部视图仍可能成功（在删边前的数据快照上）
                            search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_; // 搜索路径邻居生成统计
                            search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_; // 搜索路径构建统计
                        } else { // 删除分支下视图无效或剪枝
                            non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_; // 非搜索路径生成统计
                            non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_; // 非搜索路径构建统计

                            if (lvm_.generate_visited_neighbor_count_ == 0) // 无邻居可扩展则直接拒绝
                                direct_rejection_count_ += 1; // 直接拒绝计数
                        }
                        first_indexing_vertex_ += lvm_.first_vertex_neighbor_; // 首个索引顶点邻居指标
                        if (enable_search && is_valid) { // 需要搜索且视图有效
                            uint64_t local_result_count = sm_.search_on_reduced_query(query_graph, *orders, lvm_, gvm_); // 枚举当前状态下的匹配（删边前）
                            result_count += local_result_count * automorphism_size; // 乘轨道大小累计
                            is_searched_ = true; // 标记已执行搜索
                        }

                        lvm_.destroy_view(); // 释放局部视图

                        if (result_count >= sm_.target_number) break; // 达目标嵌入数则停
                        if (g_exit) break; // 全局退出
                    }
                } else {


                    struct PreGroupTask {
                        OrdersPerEdge* orders = nullptr;
                        uint32_t automorphism_size = 1;
                        uint32_t automorphism_id = 0;
                    };
                    struct PreGroupBucket {
                        uint32_t group_id = 0;
                        std::vector<PreGroupTask> tasks;
                        uint32_t chosen_sub_id = std::numeric_limits<uint32_t>::max();
                    };

                    std::vector<PreGroupTask> plain_search_pre;
                    std::unordered_map<uint32_t, PreGroupBucket> grouped_pre;
                    plain_search_pre.reserve(mapped_automorphism->size());
                    grouped_pre.reserve(32);

                    for (auto automorphism_id : *mapped_automorphism) {
                        auto meta = om_.get_automorphism_meta(automorphism_id);
                        OrdersPerEdge* orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta);
                        if (orders == nullptr) continue;

                        const uint32_t u0 = orders->indexing_order_[0];
                        const uint32_t u1 = orders->indexing_order_[1];

                        if (!gvm_.nlf_check(u0, data_edge.first) || !gvm_.nlf_check(u1, data_edge.second)) {
                            continue;
                        }

                        DirectedEdgeKey key{u0, u1, 0};
                        uint32_t gid = 0;
                        auto it_gid = phase1_edge_group_id_.find(key);
                        if (it_gid != phase1_edge_group_id_.end()) {
                             gid = it_gid->second;
                        }

                        if (gid == 0) {
                            plain_search_pre.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                            continue;
                        }

                        auto& b = grouped_pre[gid];
                        b.group_id = gid;
                        b.tasks.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                        b.chosen_sub_id = phase1_group_chosen_sub_id_[gid]; 
                    }

//////
                    // Part 1: plain_search_pre
                    const auto t_plain_total_begin = std::chrono::high_resolution_clock::now();
                    if (plain_search_pre.size() <= 1){
                        for (const auto &pt : plain_search_pre) {
                            if (pt.orders == nullptr) continue;
                            // 这里走“单边组”回溯：直接按当前 orders 构建局部视图并搜索（不依赖 pending candidates）。
                            auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                            bool is_valid = lvm_.create_view(query_graph, *pt.orders, gvm_, data_edge);
                            auto t_lvm_end = std::chrono::high_resolution_clock::now();
                            time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                    t_lvm_end - t_lvm_begin).count();

                            if (enable_search && is_valid) {
                                auto t_plain_begin = std::chrono::high_resolution_clock::now();
                                uint64_t local_result_count =
                                        sm_.search_on_reduced_query(query_graph, *pt.orders, lvm_, gvm_);
                                auto t_plain_end = std::chrono::high_resolution_clock::now();
                                time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                        t_plain_end - t_plain_begin).count();
                                result_count += local_result_count * pt.automorphism_size;
                                is_searched_ = true;
                            }

                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view();
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                t_lvm_destroy_end - t_lvm_destroy_begin).count();

                        if (result_count >= sm_.target_number) break;
                        if (g_exit) break;
                        }
                    }else{
                        sm_.clear_phase3_block_cache();
                        for (const auto &pt : plain_search_pre) {
                            if (pt.orders == nullptr) continue;
                            auto meta = om_.get_block_reuse_automorphism_meta(pt.automorphism_id);
                            OrdersPerEdge *orders = std::get<1>(meta);
                            uint32_t automorphism_size = std::get<2>(meta);
                            auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                            bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge);
                            auto t_lvm_end = std::chrono::high_resolution_clock::now();
                            time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_end - t_lvm_begin).count();
                            if (is_valid) {
                                search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                                search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                            } else {
                                non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                                non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;

                                if (lvm_.generate_visited_neighbor_count_ == 0)
                                    direct_rejection_count_ += 1;
                            }
                            first_indexing_vertex_ += lvm_.first_vertex_neighbor_;
                            if (enable_search && is_valid) {
                                phase3_insert_block_reuse_calls_ += 1;
                                auto t_plain_begin = std::chrono::high_resolution_clock::now();
                                uint64_t local_result_count =
                                        sm_.block_reuse_search_on_reduced_query(query_graph, *orders, lvm_, gvm_, query_bctree_);
                                auto t_plain_end = std::chrono::high_resolution_clock::now();
                                time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_plain_end - t_plain_begin).count();

                                result_count += local_result_count * automorphism_size;
                                is_searched_ = true;
                            }

                            auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                            lvm_.destroy_view();
                            auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                            time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_destroy_end - t_lvm_destroy_begin).count();

                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;
                        }
                        sm_.clear_phase3_block_cache();
                    }
                    const auto t_plain_total_end = std::chrono::high_resolution_clock::now();
                    time_phase1on_plain_total_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            t_plain_total_end - t_plain_total_begin).count();





                            
/////////////////多边组
                      
                    for (auto &kv : grouped_pre) {
                        auto &bucket = kv.second;
                        const uint32_t group_sub_id = bucket.chosen_sub_id;

                        auto *subctx = get_subquery_context(group_sub_id);
                        if (subctx == nullptr) {
                            continue;
                        }

                        OrdersPerEdge* group_seed_orders =
                                std::get<1>(om_.get_automorphism_meta_special(bucket.tasks[0].automorphism_id, group_sub_id));
                        if (group_seed_orders == nullptr) {
                            group_seed_orders = bucket.tasks[0].orders;
                        }
                        if (group_seed_orders == nullptr) {
                            continue;
                        }
                        const uint32_t u0 = group_seed_orders->indexing_order_[0];
                        const uint32_t u1 = group_seed_orders->indexing_order_[1];
                
                        const int32_t n0 = subctx->old2new[u0];
                        const int32_t n1 = subctx->old2new[u1];
                
                        const uint32_t ns = static_cast<uint32_t>(n0);
                        const uint32_t nt = static_cast<uint32_t>(n1);

                        const uint32_t sub_step0_key =
                                (u0 & 0xFFu) | ((u1 & 0xFFu) << 8) | ((ns & 0xFFu) << 16) | ((nt & 0xFFu) << 24);
                        const std::vector<std::vector<uint32_t>>* sub_mapped_new_orbit_peers = nullptr;
                        if (n0 >= 0 && n1 >= 0) {
                            auto it_orbit = subctx->mapped_new_orbit_peers_by_key.find(sub_step0_key);
                            if (it_orbit != subctx->mapped_new_orbit_peers_by_key.end())
                                sub_mapped_new_orbit_peers = &it_orbit->second;
                        }

                        OrdersPerEdge *sub_orders = subctx->om.get_orders({ns, nt});

                        if (sub_orders == nullptr) {
                            continue;
                        }
                        bool sub_valid = subctx->lvm.create_candidates_only_for_subgraph(subctx->sub_query.get(),
                                                                                               query_graph,
                                                                                               *sub_orders,
                                                                                               gvm_,
                                                                                               data_edge,
                                                                                               &subctx->new2old,
                                                                                               &subctx->sub_vertex_nlf,
                                                                                               sub_mapped_new_orbit_peers);

                        if (!sub_valid) {
                            continue;
                        }
                        phase1_before_sub_valid_check_count_ += 1;
                        LocalCandidatesSnapshot sub_candidates_snapshot =
                                subctx->lvm.extract_candidates_snapshot(subctx->sub_query.get());

                        LocalCandidatesSnapshot sub_candidates_snapshot_for_tasks = sub_candidates_snapshot;
                        subctx->lvm.restore_candidates_from(std::move(sub_candidates_snapshot));

                        std::vector<std::vector<uint32_t>> sub_embeddings_for_group;
                        bool sub_embeddings_ready = false;

                        static const std::vector<std::vector<uint32_t>> k_empty_mapped_new_by_old;

                        for (auto &pt : bucket.tasks) {
                            OrdersPerEdge* orders_for_group =
                                    std::get<1>(om_.get_automorphism_meta_special(pt.automorphism_id, group_sub_id));

                
                            const uint32_t task_u0 = orders_for_group->indexing_order_[0];
                            const uint32_t task_u1 = orders_for_group->indexing_order_[1];
                            const uint32_t task_ns = ns;
                            const uint32_t task_nt = nt;
                            const uint32_t key = (task_u0 & 0xFFu) | ((task_u1 & 0xFFu) << 8) |
                                                 ((task_ns & 0xFFu) << 16) | ((task_nt & 0xFFu) << 24);
                            auto it_mapped_new_set = subctx->mapped_new_set_by_key_and_old.find(key);

                            const std::vector<std::vector<uint32_t>> *mapped_for_task =
                                    (it_mapped_new_set != subctx->mapped_new_set_by_key_and_old.end())
                                            ? &it_mapped_new_set->second
                                            : &k_empty_mapped_new_by_old;
                          
                            bool task_pruned_valid = lvm_.create_candidates_only_with_sub_constraint(
                                    query_graph,
                                    *orders_for_group,
                                    gvm_,
                                    data_edge,
                                    sub_candidates_snapshot_for_tasks,
                                    *mapped_for_task,
                                    subctx->sub_query.get(),
                                    &subctx->new2old);


                            phase1_before_task_pruned_check_count_ += 1;
                            if (!task_pruned_valid) {
                                continue;
                            }

                            if (!sub_embeddings_ready) {
                                subctx->lvm.build_view_from_candidates_for_subgraph(
                                        subctx->sub_query.get(), *sub_orders, gvm_, &subctx->new2old);
                                sub_embeddings_for_group = subctx->sm.special_search_on_reduced_query(subctx->sub_query.get(),
                                                                                                           *sub_orders,
                                                                                                           subctx->lvm,
                                                                                                           gvm_,
                                                                                                           subctx->sm.target_number,
                                                                                                           &subctx->new2old);
                                subctx->lvm.destroy_view();
                                sub_embeddings_ready = true;
                                
                            }

                            phase1_before_sub_embeddings_check_count_ += 1;
                            if (sub_embeddings_for_group.empty()) {
                                continue;
                            }
                            phase1_after_sub_embeddings_check_count_ += 1;

                                
                            auto it_mapping_ids = subctx->mapping_ids_index.find(key);
                            const std::vector<uint32_t>& candidate_mapping_ids = it_mapping_ids->second;

                            lvm_.build_view_from_candidates(query_graph, *orders_for_group, gvm_);
                            for (const auto& emb : sub_embeddings_for_group) {
                                for (uint32_t selected_mapping_id : candidate_mapping_ids) {
                                    auto it_tpl = subctx->forced_pairs_template.find(selected_mapping_id);
                                    uint64_t local_result_count =
                                            sm_.constrained_search_with_subemb(query_graph,
                                                                                *orders_for_group,
                                                                                lvm_,
                                                                                gvm_,
                                                                                it_tpl->second,
                                                                                emb);
                                    result_count += local_result_count * pt.automorphism_size;
                                    is_searched_ = true;
                                    if (result_count >= sm_.target_number) break;
                                    if (g_exit) break;
                                }
                                if (result_count >= sm_.target_number) break;
                                if (g_exit) break;
                            }
                            lvm_.destroy_view();
                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;

                        }
                    }
                }
            }
        }

        auto it = gvm_.get_mapped_views(label_triple); 
        if (it != nullptr) { 
            is_relevant_ = true; 
            relevant_update_count_ += 1; 
            for (uint32_t i = 0; i < 2; ++i) { 
                it = gvm_.get_mapped_views(label_triple); 
                uint32_t view_id = it->second; 
                gvm_.update_view(update.op_, view_id, data_edge, label_triple); 
                std::swap(data_edge.first, data_edge.second);
                std::swap(label_triple.src_label_, label_triple.dst_label_);
            }
            std::swap(data_edge.first, data_edge.second);
            std::swap(label_triple.src_label_, label_triple.dst_label_);
        }
        if (is_relevant_) { 
            gvm_.update_nlf_view(update.op_, data_edge, label_triple); 
        }
    }
    edge_process_count_ += 1; 
    if (is_relevant_) { 
        if (result_count > 0) { 
            positive_count_ += 1; 
        }
        if (is_searched_) { 
            search_count_ += 1; 
        }
        result_count_ += result_count; 
        invalid_partial_result_count_ += sm_.invalid_partial_result_count_; 
        partial_result_count_ += sm_.partial_result_count_; 
        iso_conflict_count_ += sm_.iso_conflict_count_; 
        si_empty_count_ += sm_.si_empty_count_; 
        lc_empty_count_ += sm_.lc_empty_count_; 
        sm_.reset_performance_counters(); 
    }
    return result_count;
}

uint64_t StreamingEngine::execute_phase5_on(const Graph *query_graph, const Update &update, bool enable_local_view,
                                           bool enable_search) {
   
    LabelTriple label_triple = update.labels_; 
    Edge data_edge = update.edge_; 
    
    const LabelTriple label_triple_initial = label_triple;
    uint64_t result_count = 0; 
    is_relevant_ = false; 
    is_searched_ = false;
    if (update.op_ == '+') {       
        auto it = gvm_.get_mapped_views(label_triple); 
        if (it != nullptr) { 
            is_relevant_ = true; 
            relevant_update_count_ += 1; 
            auto t_gvm_begin = std::chrono::high_resolution_clock::now();
            for (uint32_t i = 0; i < 2; ++i) { 
                it = gvm_.get_mapped_views(label_triple); 
                uint32_t view_id = it->second; 
                gvm_.update_view(update.op_, view_id, data_edge, label_triple); 
                // 处理无向边的反向键：(v,u) 与交换后的标签三元组
                std::swap(data_edge.first, data_edge.second);
                std::swap(label_triple.src_label_, label_triple.dst_label_); 
            }
            std::swap(data_edge.first, data_edge.second); 
            std::swap(label_triple.src_label_, label_triple.dst_label_); 
            auto t_gvm_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_views_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_gvm_end - t_gvm_begin).count();
        }
        if (is_relevant_) { 
            auto t_nlf_begin = std::chrono::high_resolution_clock::now();
            gvm_.update_nlf_view(update.op_, data_edge, label_triple); 
            auto t_nlf_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_nlf_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_nlf_end - t_nlf_begin).count();
        }
   
        std::vector<Edge> *mapped_query_edges = om_.get_mapped_edges(label_triple);
        const size_t phase1_matching_query_edge_count = mapped_query_edges ? mapped_query_edges->size() : 0;
        auto mapped_automorphism = om_.get_mapped_automorphism(label_triple);
        if (mapped_automorphism != nullptr) { 
            is_relevant_ = true; 
            if (enable_local_view) { 
                if (phase1_matching_query_edge_count <= 1) { 
                    for (auto automorphism_id: *mapped_automorphism) { 
                        auto meta = om_.get_automorphism_meta(automorphism_id); 
                        OrdersPerEdge *orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta); 
                        auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                        bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge); 
                        auto t_lvm_end = std::chrono::high_resolution_clock::now();
                        time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_end - t_lvm_begin).count();
                        if (is_valid) { 
                            search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_; 
                            search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                        } else { 
                            non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_; 
                            non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_; 
                            if (lvm_.generate_visited_neighbor_count_ == 0) 
                                direct_rejection_count_ += 1; 
                        }
                        first_indexing_vertex_ += lvm_.first_vertex_neighbor_; 
                        if (enable_search && is_valid) { 
                            auto t_plain_begin = std::chrono::high_resolution_clock::now();
                            uint64_t local_result_count = sm_.search_on_reduced_query(query_graph, *orders, lvm_, gvm_); 
                            auto t_plain_end = std::chrono::high_resolution_clock::now();
                            time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_plain_end - t_plain_begin).count();
                            result_count += local_result_count * automorphism_size;
                            is_searched_ = true; 
                        }
                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view(); 
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_destroy_end - t_lvm_destroy_begin).count();

                        if (result_count >= sm_.target_number) break;
                        if (g_exit) break; 
                    }
                }else {
                    struct PreGroupTask {
                        OrdersPerEdge* orders = nullptr;
                        uint32_t automorphism_size = 1;
                        uint32_t automorphism_id = 0;
                    };
                    struct PreGroupBucket {
                        uint32_t group_id = 0;
                        std::vector<PreGroupTask> tasks;
                        uint32_t chosen_sub_id = std::numeric_limits<uint32_t>::max();
                    };

                    std::vector<PreGroupTask> plain_search_pre;
                    std::unordered_map<uint32_t, PreGroupBucket> grouped_pre;
                    std::unordered_map<uint32_t, PreGroupBucket> iso_grouped_pre;
                    plain_search_pre.reserve(mapped_automorphism->size());
                    grouped_pre.reserve(32);
                    iso_grouped_pre.reserve(32);

                    for (auto automorphism_id : *mapped_automorphism) {
                        auto meta = om_.get_automorphism_meta(automorphism_id);
                        OrdersPerEdge* orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta);
                        if (orders == nullptr) continue;

                        const uint32_t u0 = orders->indexing_order_[0];
                        const uint32_t u1 = orders->indexing_order_[1];

                        if (!gvm_.nlf_check(u0, data_edge.first) || !gvm_.nlf_check(u1, data_edge.second)) {
                            continue;
                        }

                        DirectedEdgeKey key{u0, u1, 0};
                        uint32_t gid = 0;
                        auto it_gid = phase1_edge_group_id_.find(key);
                        if (it_gid != phase1_edge_group_id_.end()) {
                             gid = it_gid->second;
                        }

                        if (gid == 0) {
                            uint32_t iso_gid = 0;
                            auto it_iso = phase2_edge_group_id_.find(key);
                            if (it_iso != phase2_edge_group_id_.end()) {
                                iso_gid = it_iso->second;
                            }
                            if (iso_gid >= 1) {
                                auto& ib = iso_grouped_pre[iso_gid];
                                ib.group_id = iso_gid;
                                ib.tasks.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                                ib.chosen_sub_id = phase2_group_chosen_sub_id_[iso_gid];
                                continue;
                            }
                            plain_search_pre.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                            continue;
                        }

                        auto& b = grouped_pre[gid];
                        b.group_id = gid;
                        b.tasks.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                        b.chosen_sub_id = phase1_group_chosen_sub_id_[gid]; 
                    }

                    for (auto it = grouped_pre.begin(); it != grouped_pre.end(); ) {
                        if (it->second.tasks.size() <= 1) {
                            for (auto& t : it->second.tasks)
                                plain_search_pre.push_back(t);
                            it = grouped_pre.erase(it);
                        } else {
                            ++it;
                        }
                    }
                    for (auto it = iso_grouped_pre.begin(); it != iso_grouped_pre.end(); ) {
                        if (it->second.tasks.size() <= 1) {
                            for (auto& t : it->second.tasks)
                                plain_search_pre.push_back(t);
                            it = iso_grouped_pre.erase(it);
                        } else {
                            ++it;
                        }
                    }

                    const auto t_plain_total_begin = std::chrono::high_resolution_clock::now();
                    if (plain_search_pre.size() <= 1){
                        for (const auto &pt : plain_search_pre) {
                            if (pt.orders == nullptr) continue;
                            auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                            bool is_valid = lvm_.create_view(query_graph, *pt.orders, gvm_, data_edge);
                            auto t_lvm_end = std::chrono::high_resolution_clock::now();
                            time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                    t_lvm_end - t_lvm_begin).count();

                            if (enable_search && is_valid) {
                                auto t_plain_begin = std::chrono::high_resolution_clock::now();
                                uint64_t local_result_count =
                                        sm_.search_on_reduced_query(query_graph, *pt.orders, lvm_, gvm_);
                                auto t_plain_end = std::chrono::high_resolution_clock::now();
                                time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                        t_plain_end - t_plain_begin).count();
                                result_count += local_result_count * pt.automorphism_size;
                                is_searched_ = true;
                            }

                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view();
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                t_lvm_destroy_end - t_lvm_destroy_begin).count();

                        if (result_count >= sm_.target_number) break;
                        if (g_exit) break;
                        }
                    }else{
                        sm_.clear_phase3_block_cache();
                        for (const auto &pt : plain_search_pre) {
                            if (pt.orders == nullptr) continue;
                            auto meta = om_.get_block_reuse_automorphism_meta(pt.automorphism_id);
                            OrdersPerEdge *orders = std::get<1>(meta);
                            uint32_t automorphism_size = std::get<2>(meta);
                            auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                            bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge);
                            auto t_lvm_end = std::chrono::high_resolution_clock::now();
                            time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_end - t_lvm_begin).count();
                            if (is_valid) {
                                search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                                search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                            } else {
                                non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                                non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;

                                if (lvm_.generate_visited_neighbor_count_ == 0)
                                    direct_rejection_count_ += 1;
                            }
                            first_indexing_vertex_ += lvm_.first_vertex_neighbor_;
                            if (enable_search && is_valid) {
                                phase3_insert_block_reuse_calls_ += 1;
                                auto t_plain_begin = std::chrono::high_resolution_clock::now();
                                uint64_t local_result_count =
                                        sm_.block_reuse_search_on_reduced_query(query_graph, *orders, lvm_, gvm_, query_bctree_);
                                auto t_plain_end = std::chrono::high_resolution_clock::now();
                                time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_plain_end - t_plain_begin).count();

                                result_count += local_result_count * automorphism_size;
                                is_searched_ = true;
                            }

                            auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                            lvm_.destroy_view();
                            auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                            time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_destroy_end - t_lvm_destroy_begin).count();

                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;
                        }
                        sm_.clear_phase3_block_cache();
                    }
                    const auto t_plain_total_end = std::chrono::high_resolution_clock::now();
                    time_phase1on_plain_total_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            t_plain_total_end - t_plain_total_begin).count();


                    for (auto &kv : iso_grouped_pre) {
                        auto &bucket = kv.second;
                        // TODO: 同构组复用搜索逻辑
                        const uint32_t group_sub_id = bucket.chosen_sub_id;
                        auto *subctx = get_isosubquery_context(group_sub_id);
                        if (subctx == nullptr) {
                            continue;
                        }

                        OrdersPerEdge* group_seed_orders =
                                std::get<1>(om_.get_automorphism_meta_special2(bucket.tasks[0].automorphism_id, group_sub_id));
                        if (group_seed_orders == nullptr) {
                            group_seed_orders = bucket.tasks[0].orders;
                        }
                        if (group_seed_orders == nullptr) {
                            continue;
                        }
                        const uint32_t u0 = group_seed_orders->indexing_order_[0];
                        const uint32_t u1 = group_seed_orders->indexing_order_[1];
                
                        const int32_t n0 = subctx->old2new[u0];
                        const int32_t n1 = subctx->old2new[u1];
                
                        const uint32_t ns = static_cast<uint32_t>(n0);
                        const uint32_t nt = static_cast<uint32_t>(n1);

                        const uint32_t sub_step0_key =
                                (u0 & 0xFFu) | ((u1 & 0xFFu) << 8) | ((ns & 0xFFu) << 16) | ((nt & 0xFFu) << 24);
                    
                        OrdersPerEdge *isosub_orders = subctx->om.get_orders({ns, nt});

                        if (isosub_orders == nullptr) {
                            continue;
                        }
                        uint32_t sub_k = subctx->sub_query->getVerticesCount();
                        bool sub_valid = subctx->lvm.create_candidates_only_for_subgraph_iso(
                                subctx->sub_query.get(),
                                query_graph,
                                *isosub_orders,
                                gvm_,
                                data_edge,
                                &subctx->old2new,
                                sub_k);
                        if (!sub_valid) {
                            continue;
                        }
                        LocalCandidatesSnapshot sub_candidates_snapshot =
                                subctx->lvm.extract_candidates_snapshot(subctx->sub_query.get());

                        LocalCandidatesSnapshot sub_candidates_snapshot_for_tasks = sub_candidates_snapshot;
                        subctx->lvm.restore_candidates_from(std::move(sub_candidates_snapshot));

                        std::vector<std::vector<uint32_t>> sub_embeddings_for_group;
                        bool sub_embeddings_ready = false;
                        static const std::vector<std::vector<uint32_t>> k_empty_mapped_new_by_old;
                        for (auto &pt : bucket.tasks) {
                            OrdersPerEdge* orders_for_group =
                                    std::get<1>(om_.get_automorphism_meta_special2(pt.automorphism_id, group_sub_id));
                            const uint32_t task_u0 = orders_for_group->indexing_order_[0];
                            const uint32_t task_u1 = orders_for_group->indexing_order_[1];

                            int32_t chosen_instance = -1;
                            for (size_t inst = 0; inst < subctx->mappings.size(); ++inst) {
                                bool has_u0 = false, has_u1 = false;
                                for (uint32_t old_u : subctx->mappings[inst]) {
                                    if (old_u == task_u0) has_u0 = true;
                                    if (old_u == task_u1) has_u1 = true;
                                }
                                if (has_u0 && has_u1) { chosen_instance = static_cast<int32_t>(inst); break; }
                            }
                            if (chosen_instance < 0) {
                                if (subctx->mappings.empty()) continue;
                                chosen_instance = 0;
                            }

                            std::vector<int32_t> local_new2old(sub_k, -1);
                            const auto& mapping = subctx->mappings[chosen_instance];
                            for (uint32_t c = 0; c < mapping.size(); ++c) {
                                local_new2old[c] = static_cast<int32_t>(mapping[c]);
                            }

                            uint32_t old_limit = static_cast<uint32_t>(subctx->old2new.size());
                            std::vector<std::vector<uint32_t>> mapped_for_task(old_limit);
                            for (uint32_t c = 0; c < mapping.size(); ++c) {
                                uint32_t old_u = mapping[c];
                                if (old_u < old_limit) {
                                    mapped_for_task[old_u].push_back(c);
                                }
                            }

                            bool task_pruned_valid = lvm_.create_candidates_only_with_sub_constraint(
                                    query_graph,
                                    *orders_for_group,
                                    gvm_,
                                    data_edge,
                                    sub_candidates_snapshot_for_tasks,
                                    mapped_for_task,
                                    subctx->sub_query.get(),
                                    &local_new2old);

                            if (!task_pruned_valid) {
                                continue;
                            }

                            if (!sub_embeddings_ready) {
                                subctx->lvm.build_view_from_candidates_for_subgraph(
                                        subctx->sub_query.get(), *isosub_orders, gvm_, &local_new2old);
                                sub_embeddings_for_group = subctx->sm.special_search_on_reduced_query(subctx->sub_query.get(),
                                                                                                           *isosub_orders,
                                                                                                           subctx->lvm,
                                                                                                           gvm_,
                                                                                                           subctx->sm.target_number,
                                                                                                           &local_new2old);
                                subctx->lvm.destroy_view();
                                sub_embeddings_ready = true;
                            }

                            if (sub_embeddings_for_group.empty()) {
                                continue;
                            }

                                
                            std::vector<std::pair<uint32_t, uint32_t>> forced_pairs;
                            forced_pairs.reserve(local_new2old.size());
                            for (size_t c = 0; c < local_new2old.size(); ++c) {
                                if (local_new2old[c] >= 0) {
                                    forced_pairs.emplace_back(static_cast<uint32_t>(local_new2old[c]), static_cast<uint32_t>(c));
                                }
                            }

                            lvm_.build_view_from_candidates(query_graph, *orders_for_group, gvm_);
                            for (const auto& emb : sub_embeddings_for_group) {
                                uint64_t local_result_count =
                                        sm_.constrained_search_with_subemb(query_graph,
                                                                            *orders_for_group,
                                                                            lvm_,
                                                                            gvm_,
                                                                            forced_pairs,
                                                                            emb);
                                result_count += local_result_count * pt.automorphism_size;
                                is_searched_ = true;

                                if (result_count >= sm_.target_number) break;
                                if (g_exit) break;
                            }
                            lvm_.destroy_view();
                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;


                        }



                    }

                    for (auto &kv : grouped_pre) {
                        auto &bucket = kv.second;
                        const uint32_t group_sub_id = bucket.chosen_sub_id;

                        auto *subctx = get_subquery_context(group_sub_id);
                        if (subctx == nullptr) {
                            continue;
                        }

                        OrdersPerEdge* group_seed_orders =
                                std::get<1>(om_.get_automorphism_meta_special(bucket.tasks[0].automorphism_id, group_sub_id));
                        if (group_seed_orders == nullptr) {
                            group_seed_orders = bucket.tasks[0].orders;
                        }
                        if (group_seed_orders == nullptr) {
                            continue;
                        }
                        const uint32_t u0 = group_seed_orders->indexing_order_[0];
                        const uint32_t u1 = group_seed_orders->indexing_order_[1];
                
                        const int32_t n0 = subctx->old2new[u0];
                        const int32_t n1 = subctx->old2new[u1];
                
                        const uint32_t ns = static_cast<uint32_t>(n0);
                        const uint32_t nt = static_cast<uint32_t>(n1);

                        const uint32_t sub_step0_key =
                                (u0 & 0xFFu) | ((u1 & 0xFFu) << 8) | ((ns & 0xFFu) << 16) | ((nt & 0xFFu) << 24);
                        const std::vector<std::vector<uint32_t>>* sub_mapped_new_orbit_peers = nullptr;
                        if (n0 >= 0 && n1 >= 0) {
                            auto it_orbit = subctx->mapped_new_orbit_peers_by_key.find(sub_step0_key);
                            if (it_orbit != subctx->mapped_new_orbit_peers_by_key.end())
                                sub_mapped_new_orbit_peers = &it_orbit->second;
                        }

                        OrdersPerEdge *sub_orders = subctx->om.get_orders({ns, nt});

                        if (sub_orders == nullptr) {
                            continue;
                        }
                        bool sub_valid = subctx->lvm.create_candidates_only_for_subgraph(subctx->sub_query.get(),
                                                                                               query_graph,
                                                                                               *sub_orders,
                                                                                               gvm_,
                                                                                               data_edge,
                                                                                               &subctx->new2old,
                                                                                               &subctx->sub_vertex_nlf,
                                                                                               sub_mapped_new_orbit_peers);

                        if (!sub_valid) {
                            continue;
                        }
                        phase1_before_sub_valid_check_count_ += 1;
                        LocalCandidatesSnapshot sub_candidates_snapshot =
                                subctx->lvm.extract_candidates_snapshot(subctx->sub_query.get());

                        LocalCandidatesSnapshot sub_candidates_snapshot_for_tasks = sub_candidates_snapshot;
                        subctx->lvm.restore_candidates_from(std::move(sub_candidates_snapshot));

                        std::vector<std::vector<uint32_t>> sub_embeddings_for_group;
                        bool sub_embeddings_ready = false;

                        static const std::vector<std::vector<uint32_t>> k_empty_mapped_new_by_old;

                        for (auto &pt : bucket.tasks) {
                            OrdersPerEdge* orders_for_group =
                                    std::get<1>(om_.get_automorphism_meta_special(pt.automorphism_id, group_sub_id));

                
                            const uint32_t task_u0 = orders_for_group->indexing_order_[0];
                            const uint32_t task_u1 = orders_for_group->indexing_order_[1];
                            const uint32_t task_ns = ns;
                            const uint32_t task_nt = nt;
                            const uint32_t key = (task_u0 & 0xFFu) | ((task_u1 & 0xFFu) << 8) |
                                                 ((task_ns & 0xFFu) << 16) | ((task_nt & 0xFFu) << 24);
                            auto it_mapped_new_set = subctx->mapped_new_set_by_key_and_old.find(key);

                            const std::vector<std::vector<uint32_t>> *mapped_for_task =
                                    (it_mapped_new_set != subctx->mapped_new_set_by_key_and_old.end())
                                            ? &it_mapped_new_set->second
                                            : &k_empty_mapped_new_by_old;
                          
                            bool task_pruned_valid = lvm_.create_candidates_only_with_sub_constraint(
                                    query_graph,
                                    *orders_for_group,
                                    gvm_,
                                    data_edge,
                                    sub_candidates_snapshot_for_tasks,
                                    *mapped_for_task,
                                    subctx->sub_query.get(),
                                    &subctx->new2old);


                            phase1_before_task_pruned_check_count_ += 1;
                            if (!task_pruned_valid) {
                                continue;
                            }

                            if (!sub_embeddings_ready) {
                                subctx->lvm.build_view_from_candidates_for_subgraph(
                                        subctx->sub_query.get(), *sub_orders, gvm_, &subctx->new2old);
                                sub_embeddings_for_group = subctx->sm.special_search_on_reduced_query(subctx->sub_query.get(),
                                                                                                           *sub_orders,
                                                                                                           subctx->lvm,
                                                                                                           gvm_,
                                                                                                           subctx->sm.target_number,
                                                                                                           &subctx->new2old);
                                subctx->lvm.destroy_view();
                                sub_embeddings_ready = true;
                                
                            }

                            phase1_before_sub_embeddings_check_count_ += 1;
                            if (sub_embeddings_for_group.empty()) {
                                continue;
                            }
                            phase1_after_sub_embeddings_check_count_ += 1;

                                
                            auto it_mapping_ids = subctx->mapping_ids_index.find(key);
                            const std::vector<uint32_t>& candidate_mapping_ids = it_mapping_ids->second;

                            lvm_.build_view_from_candidates(query_graph, *orders_for_group, gvm_);
                            for (const auto& emb : sub_embeddings_for_group) {
                                for (uint32_t selected_mapping_id : candidate_mapping_ids) {
                                    auto it_tpl = subctx->forced_pairs_template.find(selected_mapping_id);

                                    uint64_t local_result_count =
                                            sm_.constrained_search_with_subemb(query_graph,
                                                                                *orders_for_group,
                                                                                lvm_,
                                                                                gvm_,
                                                                                it_tpl->second,
                                                                                emb);
                                    result_count += local_result_count * pt.automorphism_size;
                                    is_searched_ = true;

                                    if (result_count >= sm_.target_number) break;
                                    if (g_exit) break;
                                }
                                if (result_count >= sm_.target_number) break;
                                if (g_exit) break;
                            }
                            lvm_.destroy_view();
                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;

                        }
                    }
                }
            }
        }
    }
    else {
        label_triple = label_triple_initial; 
        data_edge = update.edge_;
        std::vector<Edge> *mapped_query_edges_del = om_.get_mapped_edges(label_triple);
        const size_t phase1_matching_query_edge_count_del = mapped_query_edges_del ? mapped_query_edges_del->size() : 0;
        auto mapped_automorphism = om_.get_mapped_automorphism(label_triple); 
        if (mapped_automorphism != nullptr) { 
            is_relevant_ = true; 
            if (enable_local_view) { 
                if (phase1_matching_query_edge_count_del <= 1) { 
                    for (auto automorphism_id: *mapped_automorphism) { 
                        auto meta = om_.get_automorphism_meta(automorphism_id); 
                        OrdersPerEdge *orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta); 
                        bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge); 

                        if (is_valid) { 
                            search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_; 
                            search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_; 
                        } else { 
                            non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_; 
                            non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_; 

                            if (lvm_.generate_visited_neighbor_count_ == 0) 
                                direct_rejection_count_ += 1; 
                        }
                        first_indexing_vertex_ += lvm_.first_vertex_neighbor_; 
                        if (enable_search && is_valid) {
                            uint64_t local_result_count = sm_.search_on_reduced_query(query_graph, *orders, lvm_, gvm_);
                            result_count += local_result_count * automorphism_size; 
                            is_searched_ = true; 
                        }

                        lvm_.destroy_view(); 

                        if (result_count >= sm_.target_number) break; 
                        if (g_exit) break; 
                    }
                } else {
                     struct PreGroupTask {
                        OrdersPerEdge* orders = nullptr;
                        uint32_t automorphism_size = 1;
                        uint32_t automorphism_id = 0;
                    };
                    struct PreGroupBucket {
                        uint32_t group_id = 0;
                        std::vector<PreGroupTask> tasks;
                        uint32_t chosen_sub_id = std::numeric_limits<uint32_t>::max();
                    };

                    std::vector<PreGroupTask> plain_search_pre;
                    std::unordered_map<uint32_t, PreGroupBucket> grouped_pre;
                    std::unordered_map<uint32_t, PreGroupBucket> iso_grouped_pre;
                    plain_search_pre.reserve(mapped_automorphism->size());
                    grouped_pre.reserve(32);
                    iso_grouped_pre.reserve(32);

                    for (auto automorphism_id : *mapped_automorphism) {
                        auto meta = om_.get_automorphism_meta(automorphism_id);
                        OrdersPerEdge* orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta);
                        if (orders == nullptr) continue;

                        const uint32_t u0 = orders->indexing_order_[0];
                        const uint32_t u1 = orders->indexing_order_[1];

                        if (!gvm_.nlf_check(u0, data_edge.first) || !gvm_.nlf_check(u1, data_edge.second)) {
                            continue;
                        }

                        DirectedEdgeKey key{u0, u1, 0};
                        uint32_t gid = 0;
                        auto it_gid = phase1_edge_group_id_.find(key);
                        if (it_gid != phase1_edge_group_id_.end()) {
                             gid = it_gid->second;
                        }

                        if (gid == 0) {
                            uint32_t iso_gid = 0;
                            auto it_iso = phase2_edge_group_id_.find(key);
                            if (it_iso != phase2_edge_group_id_.end()) {
                                iso_gid = it_iso->second;
                            }
                            if (iso_gid >= 1) {
                                auto& ib = iso_grouped_pre[iso_gid];
                                ib.group_id = iso_gid;
                                ib.tasks.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                                ib.chosen_sub_id = phase2_group_chosen_sub_id_[iso_gid];
                                continue;
                            }
                            plain_search_pre.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                            continue;
                        }

                        auto& b = grouped_pre[gid];
                        b.group_id = gid;
                        b.tasks.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                        b.chosen_sub_id = phase1_group_chosen_sub_id_[gid]; 
                    }

                    for (auto it = grouped_pre.begin(); it != grouped_pre.end(); ) {
                        if (it->second.tasks.size() <= 1) {
                            for (auto& t : it->second.tasks)
                                plain_search_pre.push_back(t);
                            it = grouped_pre.erase(it);
                        } else {
                            ++it;
                        }
                    }
                    for (auto it = iso_grouped_pre.begin(); it != iso_grouped_pre.end(); ) {
                        if (it->second.tasks.size() <= 1) {
                            for (auto& t : it->second.tasks)
                                plain_search_pre.push_back(t);
                            it = iso_grouped_pre.erase(it);
                        } else {
                            ++it;
                        }
                    }

                  const auto t_plain_total_begin = std::chrono::high_resolution_clock::now();
                    if (plain_search_pre.size() <= 1){
                        for (const auto &pt : plain_search_pre) {
                            if (pt.orders == nullptr) continue;
                            auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                            bool is_valid = lvm_.create_view(query_graph, *pt.orders, gvm_, data_edge);
                            auto t_lvm_end = std::chrono::high_resolution_clock::now();
                            time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                    t_lvm_end - t_lvm_begin).count();

                            if (enable_search && is_valid) {
                                auto t_plain_begin = std::chrono::high_resolution_clock::now();
                                uint64_t local_result_count =
                                        sm_.search_on_reduced_query(query_graph, *pt.orders, lvm_, gvm_);
                                auto t_plain_end = std::chrono::high_resolution_clock::now();
                                time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                        t_plain_end - t_plain_begin).count();
                                result_count += local_result_count * pt.automorphism_size;
                                is_searched_ = true;
                            }

                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view();
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                t_lvm_destroy_end - t_lvm_destroy_begin).count();

                        if (result_count >= sm_.target_number) break;
                        if (g_exit) break;
                        }
                    }else{
                        sm_.clear_phase3_block_cache();
                        for (const auto &pt : plain_search_pre) {
                            if (pt.orders == nullptr) continue;
                            auto meta = om_.get_block_reuse_automorphism_meta(pt.automorphism_id);
                            OrdersPerEdge *orders = std::get<1>(meta);
                            uint32_t automorphism_size = std::get<2>(meta);
                            auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                            bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge);
                            auto t_lvm_end = std::chrono::high_resolution_clock::now();
                            time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_end - t_lvm_begin).count();
                            if (is_valid) {
                                search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                                search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                            } else {
                                non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                                non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;

                                if (lvm_.generate_visited_neighbor_count_ == 0)
                                    direct_rejection_count_ += 1;
                            }
                            first_indexing_vertex_ += lvm_.first_vertex_neighbor_;
                            if (enable_search && is_valid) {
                                phase3_insert_block_reuse_calls_ += 1;
                                auto t_plain_begin = std::chrono::high_resolution_clock::now();
                                uint64_t local_result_count =
                                        sm_.block_reuse_search_on_reduced_query(query_graph, *orders, lvm_, gvm_, query_bctree_);
                                auto t_plain_end = std::chrono::high_resolution_clock::now();
                                time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_plain_end - t_plain_begin).count();

                                result_count += local_result_count * automorphism_size;
                                is_searched_ = true;
                            }

                            auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                            lvm_.destroy_view();
                            auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                            time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_destroy_end - t_lvm_destroy_begin).count();

                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;
                        }
                        sm_.clear_phase3_block_cache();
                    }
                    const auto t_plain_total_end = std::chrono::high_resolution_clock::now();
                    time_phase1on_plain_total_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            t_plain_total_end - t_plain_total_begin).count();


                    for (auto &kv : iso_grouped_pre) {
                        auto &bucket = kv.second;
                        // TODO: 同构组复用搜索逻辑
                        const uint32_t group_sub_id = bucket.chosen_sub_id;
                        auto *subctx = get_isosubquery_context(group_sub_id);
                        if (subctx == nullptr) {
                            continue;
                        }

                        OrdersPerEdge* group_seed_orders =
                                std::get<1>(om_.get_automorphism_meta_special2(bucket.tasks[0].automorphism_id, group_sub_id));
                        if (group_seed_orders == nullptr) {
                            group_seed_orders = bucket.tasks[0].orders;
                        }
                        if (group_seed_orders == nullptr) {
                            continue;
                        }
                        const uint32_t u0 = group_seed_orders->indexing_order_[0];
                        const uint32_t u1 = group_seed_orders->indexing_order_[1];
                
                        const int32_t n0 = subctx->old2new[u0];
                        const int32_t n1 = subctx->old2new[u1];
                
                        const uint32_t ns = static_cast<uint32_t>(n0);
                        const uint32_t nt = static_cast<uint32_t>(n1);

                        const uint32_t sub_step0_key =
                                (u0 & 0xFFu) | ((u1 & 0xFFu) << 8) | ((ns & 0xFFu) << 16) | ((nt & 0xFFu) << 24);
                    
                        OrdersPerEdge *isosub_orders = subctx->om.get_orders({ns, nt});

                        if (isosub_orders == nullptr) {
                            continue;
                        }
                        uint32_t sub_k = subctx->sub_query->getVerticesCount();
                        bool sub_valid = subctx->lvm.create_candidates_only_for_subgraph_iso(
                                subctx->sub_query.get(),
                                query_graph,
                                *isosub_orders,
                                gvm_,
                                data_edge,
                                &subctx->old2new,
                                sub_k);
                        if (!sub_valid) {
                            continue;
                        }
                        LocalCandidatesSnapshot sub_candidates_snapshot =
                                subctx->lvm.extract_candidates_snapshot(subctx->sub_query.get());

                        LocalCandidatesSnapshot sub_candidates_snapshot_for_tasks = sub_candidates_snapshot;
                        subctx->lvm.restore_candidates_from(std::move(sub_candidates_snapshot));

                        std::vector<std::vector<uint32_t>> sub_embeddings_for_group;
                        bool sub_embeddings_ready = false;
                        static const std::vector<std::vector<uint32_t>> k_empty_mapped_new_by_old;
                        for (auto &pt : bucket.tasks) {
                            OrdersPerEdge* orders_for_group =
                                    std::get<1>(om_.get_automorphism_meta_special2(pt.automorphism_id, group_sub_id));
                            const uint32_t task_u0 = orders_for_group->indexing_order_[0];
                            const uint32_t task_u1 = orders_for_group->indexing_order_[1];

                            int32_t chosen_instance = -1;
                            for (size_t inst = 0; inst < subctx->mappings.size(); ++inst) {
                                bool has_u0 = false, has_u1 = false;
                                for (uint32_t old_u : subctx->mappings[inst]) {
                                    if (old_u == task_u0) has_u0 = true;
                                    if (old_u == task_u1) has_u1 = true;
                                }
                                if (has_u0 && has_u1) { chosen_instance = static_cast<int32_t>(inst); break; }
                            }
                            if (chosen_instance < 0) {
                                if (subctx->mappings.empty()) continue;
                                chosen_instance = 0;
                            }

                            std::vector<int32_t> local_new2old(sub_k, -1);
                            const auto& mapping = subctx->mappings[chosen_instance];
                            for (uint32_t c = 0; c < mapping.size(); ++c) {
                                local_new2old[c] = static_cast<int32_t>(mapping[c]);
                            }

                            uint32_t old_limit = static_cast<uint32_t>(subctx->old2new.size());
                            std::vector<std::vector<uint32_t>> mapped_for_task(old_limit);
                            for (uint32_t c = 0; c < mapping.size(); ++c) {
                                uint32_t old_u = mapping[c];
                                if (old_u < old_limit) {
                                    mapped_for_task[old_u].push_back(c);
                                }
                            }

                            bool task_pruned_valid = lvm_.create_candidates_only_with_sub_constraint(
                                    query_graph,
                                    *orders_for_group,
                                    gvm_,
                                    data_edge,
                                    sub_candidates_snapshot_for_tasks,
                                    mapped_for_task,
                                    subctx->sub_query.get(),
                                    &local_new2old);

                            if (!task_pruned_valid) {
                                continue;
                            }

                            if (!sub_embeddings_ready) {
                                subctx->lvm.build_view_from_candidates_for_subgraph(
                                        subctx->sub_query.get(), *isosub_orders, gvm_, &local_new2old);
                                sub_embeddings_for_group = subctx->sm.special_search_on_reduced_query(subctx->sub_query.get(),
                                                                                                           *isosub_orders,
                                                                                                           subctx->lvm,
                                                                                                           gvm_,
                                                                                                           subctx->sm.target_number,
                                                                                                           &local_new2old);
                                subctx->lvm.destroy_view();
                                sub_embeddings_ready = true;
                            }

                            if (sub_embeddings_for_group.empty()) {
                                continue;
                            }

                                
                            std::vector<std::pair<uint32_t, uint32_t>> forced_pairs;
                            forced_pairs.reserve(local_new2old.size());
                            for (size_t c = 0; c < local_new2old.size(); ++c) {
                                if (local_new2old[c] >= 0) {
                                    forced_pairs.emplace_back(static_cast<uint32_t>(local_new2old[c]), static_cast<uint32_t>(c));
                                }
                            }

                            lvm_.build_view_from_candidates(query_graph, *orders_for_group, gvm_);
                            for (const auto& emb : sub_embeddings_for_group) {
                                uint64_t local_result_count =
                                        sm_.constrained_search_with_subemb(query_graph,
                                                                            *orders_for_group,
                                                                            lvm_,
                                                                            gvm_,
                                                                            forced_pairs,
                                                                            emb);
                                result_count += local_result_count * pt.automorphism_size;
                                is_searched_ = true;

                                if (result_count >= sm_.target_number) break;
                                if (g_exit) break;
                            }
                            lvm_.destroy_view();
                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;


                        }



                    }

                    for (auto &kv : grouped_pre) {
                        auto &bucket = kv.second;
                        const uint32_t group_sub_id = bucket.chosen_sub_id;

                        auto *subctx = get_subquery_context(group_sub_id);
                        if (subctx == nullptr) {
                            continue;
                        }

                        OrdersPerEdge* group_seed_orders =
                                std::get<1>(om_.get_automorphism_meta_special(bucket.tasks[0].automorphism_id, group_sub_id));
                        if (group_seed_orders == nullptr) {
                            group_seed_orders = bucket.tasks[0].orders;
                        }
                        if (group_seed_orders == nullptr) {
                            continue;
                        }
                        const uint32_t u0 = group_seed_orders->indexing_order_[0];
                        const uint32_t u1 = group_seed_orders->indexing_order_[1];
                
                        const int32_t n0 = subctx->old2new[u0];
                        const int32_t n1 = subctx->old2new[u1];
                
                        const uint32_t ns = static_cast<uint32_t>(n0);
                        const uint32_t nt = static_cast<uint32_t>(n1);

                        const uint32_t sub_step0_key =
                                (u0 & 0xFFu) | ((u1 & 0xFFu) << 8) | ((ns & 0xFFu) << 16) | ((nt & 0xFFu) << 24);
                        const std::vector<std::vector<uint32_t>>* sub_mapped_new_orbit_peers = nullptr;
                        if (n0 >= 0 && n1 >= 0) {
                            auto it_orbit = subctx->mapped_new_orbit_peers_by_key.find(sub_step0_key);
                            if (it_orbit != subctx->mapped_new_orbit_peers_by_key.end())
                                sub_mapped_new_orbit_peers = &it_orbit->second;
                        }

                        OrdersPerEdge *sub_orders = subctx->om.get_orders({ns, nt});

                        if (sub_orders == nullptr) {
                            continue;
                        }
                        bool sub_valid = subctx->lvm.create_candidates_only_for_subgraph(subctx->sub_query.get(),
                                                                                               query_graph,
                                                                                               *sub_orders,
                                                                                               gvm_,
                                                                                               data_edge,
                                                                                               &subctx->new2old,
                                                                                               &subctx->sub_vertex_nlf,
                                                                                               sub_mapped_new_orbit_peers);

                        if (!sub_valid) {
                            continue;
                        }
                        phase1_before_sub_valid_check_count_ += 1;
                        LocalCandidatesSnapshot sub_candidates_snapshot =
                                subctx->lvm.extract_candidates_snapshot(subctx->sub_query.get());

                        LocalCandidatesSnapshot sub_candidates_snapshot_for_tasks = sub_candidates_snapshot;
                        subctx->lvm.restore_candidates_from(std::move(sub_candidates_snapshot));

                        std::vector<std::vector<uint32_t>> sub_embeddings_for_group;
                        bool sub_embeddings_ready = false;

                        static const std::vector<std::vector<uint32_t>> k_empty_mapped_new_by_old;

                        for (auto &pt : bucket.tasks) {
                            OrdersPerEdge* orders_for_group =
                                    std::get<1>(om_.get_automorphism_meta_special(pt.automorphism_id, group_sub_id));

                
                            const uint32_t task_u0 = orders_for_group->indexing_order_[0];
                            const uint32_t task_u1 = orders_for_group->indexing_order_[1];
                            const uint32_t task_ns = ns;
                            const uint32_t task_nt = nt;
                            const uint32_t key = (task_u0 & 0xFFu) | ((task_u1 & 0xFFu) << 8) |
                                                 ((task_ns & 0xFFu) << 16) | ((task_nt & 0xFFu) << 24);
                            auto it_mapped_new_set = subctx->mapped_new_set_by_key_and_old.find(key);

                            const std::vector<std::vector<uint32_t>> *mapped_for_task =
                                    (it_mapped_new_set != subctx->mapped_new_set_by_key_and_old.end())
                                            ? &it_mapped_new_set->second
                                            : &k_empty_mapped_new_by_old;
                          
                            bool task_pruned_valid = lvm_.create_candidates_only_with_sub_constraint(
                                    query_graph,
                                    *orders_for_group,
                                    gvm_,
                                    data_edge,
                                    sub_candidates_snapshot_for_tasks,
                                    *mapped_for_task,
                                    subctx->sub_query.get(),
                                    &subctx->new2old);


                            phase1_before_task_pruned_check_count_ += 1;
                            if (!task_pruned_valid) {
                                continue;
                            }

                            if (!sub_embeddings_ready) {
                                subctx->lvm.build_view_from_candidates_for_subgraph(
                                        subctx->sub_query.get(), *sub_orders, gvm_, &subctx->new2old);
                                sub_embeddings_for_group = subctx->sm.special_search_on_reduced_query(subctx->sub_query.get(),
                                                                                                           *sub_orders,
                                                                                                           subctx->lvm,
                                                                                                           gvm_,
                                                                                                           subctx->sm.target_number,
                                                                                                           &subctx->new2old);
                                subctx->lvm.destroy_view();
                                sub_embeddings_ready = true;
                                
                            }

                            phase1_before_sub_embeddings_check_count_ += 1;
                            if (sub_embeddings_for_group.empty()) {
                                continue;
                            }
                            phase1_after_sub_embeddings_check_count_ += 1;

                                
                            auto it_mapping_ids = subctx->mapping_ids_index.find(key);
                            const std::vector<uint32_t>& candidate_mapping_ids = it_mapping_ids->second;

                            lvm_.build_view_from_candidates(query_graph, *orders_for_group, gvm_);
                            for (const auto& emb : sub_embeddings_for_group) {
                                for (uint32_t selected_mapping_id : candidate_mapping_ids) {
                                    auto it_tpl = subctx->forced_pairs_template.find(selected_mapping_id);

                                    uint64_t local_result_count =
                                            sm_.constrained_search_with_subemb(query_graph,
                                                                                *orders_for_group,
                                                                                lvm_,
                                                                                gvm_,
                                                                                it_tpl->second,
                                                                                emb);
                                    result_count += local_result_count * pt.automorphism_size;
                                    is_searched_ = true;

                                    if (result_count >= sm_.target_number) break;
                                    if (g_exit) break;
                                }
                                if (result_count >= sm_.target_number) break;
                                if (g_exit) break;
                            }
                            lvm_.destroy_view();
                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;

                        }
                    }

                }
            }
        }

        auto it = gvm_.get_mapped_views(label_triple); // 局部处理完成后，再更新 global view（删除边上的匹配结构）
        if (it != nullptr) { // 查询中存在该标签模式对应的 global view
            is_relevant_ = true; // 本更新与全局视图维护相关
            relevant_update_count_ += 1; // 相关更新统计 +1
            // [无向边核对] 删除边同样在两种有向朝向下各更新一次 global view，再恢复 data_edge/label_triple。
            for (uint32_t i = 0; i < 2; ++i) { // 同样两种有向朝向各删一次
                it = gvm_.get_mapped_views(label_triple); // 按当前朝向重新取映射

                uint32_t view_id = it->second; // 当前朝向下 global view 编号
                gvm_.update_view(update.op_, view_id, data_edge, label_triple); // op 为 '-'，从对应 global view 中移除该边

                std::swap(data_edge.first, data_edge.second); // 换向以更新反向 global view
                std::swap(label_triple.src_label_, label_triple.dst_label_); // 同步交换端点标签
            }

            std::swap(data_edge.first, data_edge.second); // 恢复 data_edge 为 update 原始顺序
            std::swap(label_triple.src_label_, label_triple.dst_label_); // 恢复 label_triple
        }

        if (is_relevant_) { // 若曾命中 global view，则 NLF 也需随删除更新
            gvm_.update_nlf_view(update.op_, data_edge, label_triple); // 删除边：更新 NLF
        }
    }
    edge_process_count_ += 1; 
    if (is_relevant_) { 
        if (result_count > 0) { 
            positive_count_ += 1; 
        }
        if (is_searched_) {
            search_count_ += 1; 
        }

        result_count_ += result_count; 
        invalid_partial_result_count_ += sm_.invalid_partial_result_count_; 
        partial_result_count_ += sm_.partial_result_count_; 
        iso_conflict_count_ += sm_.iso_conflict_count_; 
        si_empty_count_ += sm_.si_empty_count_; 
        lc_empty_count_ += sm_.lc_empty_count_; 

        sm_.reset_performance_counters(); 
    }
    return result_count; 
}


uint64_t StreamingEngine::execute_phase1_off(const Graph *query_graph, const Update &update, bool enable_local_view,
                                             bool enable_search) {
    LabelTriple label_triple = update.labels_; 
    Edge data_edge = update.edge_;           

    uint64_t result_count = 0;  
    is_relevant_ = false;      
    is_searched_ = false;       

    if (update.op_ == '+') { 
        auto it = gvm_.get_mapped_views(label_triple); 
        if (it != nullptr) {
            is_relevant_ = true;
            relevant_update_count_ += 1;
            auto t_gvm_begin = std::chrono::high_resolution_clock::now();
            for (uint32_t i = 0; i < 2; ++i) {
                it = gvm_.get_mapped_views(label_triple);
                uint32_t view_id = it->second;
                gvm_.update_view(update.op_, view_id, data_edge, label_triple);
                std::swap(data_edge.first, data_edge.second);
                std::swap(label_triple.src_label_, label_triple.dst_label_);
            }
            std::swap(data_edge.first, data_edge.second);
            std::swap(label_triple.src_label_, label_triple.dst_label_); 
            auto t_gvm_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_views_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_gvm_end - t_gvm_begin).count();
        }
        if (is_relevant_) {
            auto t_nlf_begin = std::chrono::high_resolution_clock::now();
            gvm_.update_nlf_view(update.op_, data_edge, label_triple);
            auto t_nlf_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_nlf_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_nlf_end - t_nlf_begin).count();
        }
        auto mapped_automorphism = om_.get_mapped_automorphism(label_triple);
        if (mapped_automorphism != nullptr) {
            is_relevant_ = true;
            if (enable_local_view) {
                for (auto automorphism_id: *mapped_automorphism) {
                    auto meta = om_.get_automorphism_meta(automorphism_id);
                    OrdersPerEdge *orders = std::get<1>(meta);
                    uint32_t automorphism_size = std::get<2>(meta);
                    auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                    bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge);
                    auto t_lvm_end = std::chrono::high_resolution_clock::now();
                    time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_end - t_lvm_begin).count();
                    if (is_valid) {
                        search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                        search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                    } else {
                        non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                        non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                        if (lvm_.generate_visited_neighbor_count_ == 0)
                            direct_rejection_count_ += 1;
                    }
                    first_indexing_vertex_ += lvm_.first_vertex_neighbor_;
                    if (enable_search && is_valid) {
                        auto t_plain_begin = std::chrono::high_resolution_clock::now();
                        uint64_t local_result_count = sm_.search_on_reduced_query(query_graph, *orders, lvm_, gvm_);
                        auto t_plain_end = std::chrono::high_resolution_clock::now();
                        const auto dt_plain = std::chrono::duration_cast<std::chrono::nanoseconds>(t_plain_end - t_plain_begin).count();
                        time_plain_search_ns_ += dt_plain;
                        time_phase1off_backtracking_search_ns_ += dt_plain;
                        result_count += local_result_count * automorphism_size;
                        is_searched_ = true;
                    }
                    auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                    lvm_.destroy_view();
                    auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                    time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_destroy_end - t_lvm_destroy_begin).count();

                    if (result_count >= sm_.target_number) break;
                    if (g_exit) break;
                }
            }
        }
    }
    else {
        auto mapped_automorphism = om_.get_mapped_automorphism(label_triple);
        if (mapped_automorphism != nullptr) {
            is_relevant_ = true;
            if (enable_local_view) {
                for (auto automorphism_id: *mapped_automorphism) {
                    auto meta = om_.get_automorphism_meta(automorphism_id);
                    OrdersPerEdge *orders = std::get<1>(meta);
                    uint32_t automorphism_size = std::get<2>(meta);
                    auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                    bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge);
                    auto t_lvm_end = std::chrono::high_resolution_clock::now();
                    time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_end - t_lvm_begin).count();
                    if (is_valid) {
                        search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                        search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                    } else {
                        non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                        non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;

                        if (lvm_.generate_visited_neighbor_count_ == 0)
                            direct_rejection_count_ += 1;
                    }
                    first_indexing_vertex_ += lvm_.first_vertex_neighbor_;
                    if (enable_search && is_valid) {
                        auto t_plain_begin = std::chrono::high_resolution_clock::now();
                        uint64_t local_result_count = sm_.search_on_reduced_query(query_graph, *orders, lvm_, gvm_);
                        auto t_plain_end = std::chrono::high_resolution_clock::now();
                        const auto dt_plain = std::chrono::duration_cast<std::chrono::nanoseconds>(t_plain_end - t_plain_begin).count();
                        time_plain_search_ns_ += dt_plain;
                        time_phase1off_backtracking_search_ns_ += dt_plain;
                        result_count += local_result_count * automorphism_size;
                        is_searched_ = true;
                    }

                    auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                    lvm_.destroy_view();
                    auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                    time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_destroy_end - t_lvm_destroy_begin).count();

                    if (result_count >= sm_.target_number) break;
                    if (g_exit) break;
                }
            }
        }

        auto it = gvm_.get_mapped_views(label_triple);
        if (it != nullptr) {
            is_relevant_ = true;
            relevant_update_count_ += 1;
            auto t_gvm_begin = std::chrono::high_resolution_clock::now();
            for (uint32_t i = 0; i < 2; ++i) {
                it = gvm_.get_mapped_views(label_triple);

                uint32_t view_id = it->second;
                gvm_.update_view(update.op_, view_id, data_edge, label_triple);

                std::swap(data_edge.first, data_edge.second);
                std::swap(label_triple.src_label_, label_triple.dst_label_);
            }

            std::swap(data_edge.first, data_edge.second);
            std::swap(label_triple.src_label_, label_triple.dst_label_);
            auto t_gvm_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_views_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_gvm_end - t_gvm_begin).count();
        }

        if (is_relevant_) {
            auto t_nlf_begin = std::chrono::high_resolution_clock::now();
            gvm_.update_nlf_view(update.op_, data_edge, label_triple);
            auto t_nlf_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_nlf_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_nlf_end - t_nlf_begin).count();
        }
    }
    edge_process_count_ += 1;
    if (is_relevant_) {
        if (result_count > 0) {
            positive_count_ += 1;
        }
        if (is_searched_) {
            search_count_ += 1;
        }
        result_count_ += result_count;
        invalid_partial_result_count_ += sm_.invalid_partial_result_count_;
        partial_result_count_ += sm_.partial_result_count_;
        iso_conflict_count_ += sm_.iso_conflict_count_;
        si_empty_count_ += sm_.si_empty_count_;
        lc_empty_count_ += sm_.lc_empty_count_;
        sm_.reset_performance_counters();
    }
    return result_count;
}

uint64_t StreamingEngine::execute_phase3_on(const Graph *query_graph, const Update &update, bool enable_local_view,
                                            bool enable_search) {
    LabelTriple label_triple = update.labels_; 
    Edge data_edge = update.edge_;            

    uint64_t result_count = 0;  
    is_relevant_ = false;        
    is_searched_ = false;        
    if (update.op_ == '+') { 
        auto it = gvm_.get_mapped_views(label_triple); 
        if (it != nullptr) {
            is_relevant_ = true;
            relevant_update_count_ += 1;
            auto t_gvm_begin = std::chrono::high_resolution_clock::now();
            for (uint32_t i = 0; i < 2; ++i) {
                it = gvm_.get_mapped_views(label_triple);
                uint32_t view_id = it->second;
                gvm_.update_view(update.op_, view_id, data_edge, label_triple);
                std::swap(data_edge.first, data_edge.second);
                std::swap(label_triple.src_label_, label_triple.dst_label_);
            }
            std::swap(data_edge.first, data_edge.second);
            std::swap(label_triple.src_label_, label_triple.dst_label_); 
            auto t_gvm_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_views_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_gvm_end - t_gvm_begin).count();
        }
        if (is_relevant_) {
            auto t_nlf_begin = std::chrono::high_resolution_clock::now();
            gvm_.update_nlf_view(update.op_, data_edge, label_triple);
            auto t_nlf_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_nlf_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_nlf_end - t_nlf_begin).count();
        }
        std::vector<Edge> *mapped_query_edges = om_.get_mapped_edges(label_triple);
        const size_t phase3_matching_query_edge_count = mapped_query_edges ? mapped_query_edges->size() : 0;
        auto mapped_automorphism = om_.get_mapped_automorphism(label_triple);
        if (mapped_automorphism != nullptr) {
            is_relevant_ = true;
            if (enable_local_view) {
                if (phase3_matching_query_edge_count <= 1) {
                    for (auto automorphism_id: *mapped_automorphism) {
                        auto meta = om_.get_automorphism_meta(automorphism_id);
                        OrdersPerEdge *orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta);
                        auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                        bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge);
                        auto t_lvm_end = std::chrono::high_resolution_clock::now();
                        time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_end - t_lvm_begin).count();


                        if (is_valid) {
                            search_generate_neighbor_count_ += +lvm_.generate_visited_neighbor_count_;
                            search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                        } else {
                            non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                            non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;

                            if (lvm_.generate_visited_neighbor_count_ == 0)
                                direct_rejection_count_ += 1;
                        }
                        first_indexing_vertex_ += lvm_.first_vertex_neighbor_;
                        if (enable_search && is_valid) {
                            phase3_insert_plain_search_le1_calls_ += 1;
                            auto t_plain_begin = std::chrono::high_resolution_clock::now();
                            uint64_t local_result_count = sm_.search_on_reduced_query(query_graph, *orders, lvm_, gvm_);
                            auto t_plain_end = std::chrono::high_resolution_clock::now();
                            time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_plain_end - t_plain_begin).count();
                            result_count += local_result_count * automorphism_size;
                            is_searched_ = true;
                        }

                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view();
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_destroy_end - t_lvm_destroy_begin).count();

                        if (result_count >= sm_.target_number) break;
                        if (g_exit) break;
                    }
                } 
                else {
                    sm_.clear_phase3_block_cache();
                    for (auto automorphism_id: *mapped_automorphism) {
                        auto meta = om_.get_block_reuse_automorphism_meta(automorphism_id);
                        OrdersPerEdge *orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta);

                        auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                        bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge);
                        auto t_lvm_end = std::chrono::high_resolution_clock::now();
                        time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_end - t_lvm_begin).count();

                        if (is_valid) {
                            search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                            search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                        } else {
                            non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                            non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;

                            if (lvm_.generate_visited_neighbor_count_ == 0)
                                direct_rejection_count_ += 1;
                        }
                        first_indexing_vertex_ += lvm_.first_vertex_neighbor_;
                        if (enable_search && is_valid) {
                            phase3_insert_block_reuse_calls_ += 1;
                            auto t_plain_begin = std::chrono::high_resolution_clock::now();
                            uint64_t local_result_count =
                                    sm_.block_reuse_search_on_reduced_query(query_graph, *orders, lvm_, gvm_, query_bctree_);
                            auto t_plain_end = std::chrono::high_resolution_clock::now();
                            time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_plain_end - t_plain_begin).count();

                            result_count += local_result_count * automorphism_size;
                            is_searched_ = true;
                        }
                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view();
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_destroy_end - t_lvm_destroy_begin).count();

                        if (result_count >= sm_.target_number) break;
                        if (g_exit) break;
                    }
                    sm_.clear_phase3_block_cache();
                }
            }
        }
    }

    else {
        std::vector<Edge> *mapped_query_edges = om_.get_mapped_edges(label_triple);
        const size_t phase3_matching_query_edge_count = mapped_query_edges ? mapped_query_edges->size() : 0;
        auto mapped_automorphism = om_.get_mapped_automorphism(label_triple);
        if (mapped_automorphism != nullptr) {
            is_relevant_ = true;
            if (enable_local_view) {
                if (phase3_matching_query_edge_count <= 1) {
                    for (auto automorphism_id: *mapped_automorphism) {
                        auto meta = om_.get_automorphism_meta(automorphism_id);
                        OrdersPerEdge *orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta);
                        auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                        bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge);
                        auto t_lvm_end = std::chrono::high_resolution_clock::now();
                        time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_end - t_lvm_begin).count();


                        if (is_valid) {
                            search_generate_neighbor_count_ += +lvm_.generate_visited_neighbor_count_;
                            search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                        } else {
                            non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                            non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;

                            if (lvm_.generate_visited_neighbor_count_ == 0)
                                direct_rejection_count_ += 1;
                        }
                        first_indexing_vertex_ += lvm_.first_vertex_neighbor_;
                        if (enable_search && is_valid) {
                            phase3_insert_plain_search_le1_calls_ += 1;
                            auto t_plain_begin = std::chrono::high_resolution_clock::now();
                            uint64_t local_result_count = sm_.search_on_reduced_query(query_graph, *orders, lvm_, gvm_);
                            auto t_plain_end = std::chrono::high_resolution_clock::now();
                            time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_plain_end - t_plain_begin).count();
                            result_count += local_result_count * automorphism_size;
                            is_searched_ = true;
                        }

                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view();
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_destroy_end - t_lvm_destroy_begin).count();

                        if (result_count >= sm_.target_number) break;
                        if (g_exit) break;
                    }
                } 
                else {
                    sm_.clear_phase3_block_cache();
                    for (auto automorphism_id: *mapped_automorphism) {
                        auto meta = om_.get_block_reuse_automorphism_meta(automorphism_id);
                        OrdersPerEdge *orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta);

                        auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                        bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge);
                        auto t_lvm_end = std::chrono::high_resolution_clock::now();
                        time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_end - t_lvm_begin).count();

                        if (is_valid) {
                            search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                            search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;
                        } else {
                            non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                            non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_;

                            if (lvm_.generate_visited_neighbor_count_ == 0)
                                direct_rejection_count_ += 1;
                        }
                        first_indexing_vertex_ += lvm_.first_vertex_neighbor_;
                        if (enable_search && is_valid) {
                            phase3_insert_block_reuse_calls_ += 1;
                            auto t_plain_begin = std::chrono::high_resolution_clock::now();
                            uint64_t local_result_count =
                                    sm_.block_reuse_search_on_reduced_query(query_graph, *orders, lvm_, gvm_, query_bctree_);
                            auto t_plain_end = std::chrono::high_resolution_clock::now();
                            time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_plain_end - t_plain_begin).count();

                            result_count += local_result_count * automorphism_size;
                            is_searched_ = true;
                        }
                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view();
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_destroy_end - t_lvm_destroy_begin).count();

                        if (result_count >= sm_.target_number) break;
                        if (g_exit) break;
                    }
                    sm_.clear_phase3_block_cache();
                }
            }
        }

        auto it = gvm_.get_mapped_views(label_triple);
        if (it != nullptr) {
            is_relevant_ = true;
            relevant_update_count_ += 1;
            // [无向边核对] 同插入分支：无向更新双向各维护一次 view，再 swap 回与 NLF/后续一致。
            auto t_gvm_begin = std::chrono::high_resolution_clock::now();
            for (uint32_t i = 0; i < 2; ++i) {
                it = gvm_.get_mapped_views(label_triple);

                uint32_t view_id = it->second;
                gvm_.update_view(update.op_, view_id, data_edge, label_triple);

                std::swap(data_edge.first, data_edge.second);
                std::swap(label_triple.src_label_, label_triple.dst_label_);
            }

            std::swap(data_edge.first, data_edge.second);
            std::swap(label_triple.src_label_, label_triple.dst_label_);
            auto t_gvm_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_views_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_gvm_end - t_gvm_begin).count();
        }

        // Update nlf view.
        if (is_relevant_) {
            auto t_nlf_begin = std::chrono::high_resolution_clock::now();
            gvm_.update_nlf_view(update.op_, data_edge, label_triple);
            auto t_nlf_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_nlf_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_nlf_end - t_nlf_begin).count();
        }
    }
    // Update performance counters
    edge_process_count_ += 1;
    if (is_relevant_) {

        if (result_count > 0) {
            positive_count_ += 1;
        }

        if (is_searched_) {
            search_count_ += 1;
        }

        result_count_ += result_count;
        invalid_partial_result_count_ += sm_.invalid_partial_result_count_;
        partial_result_count_ += sm_.partial_result_count_;
        iso_conflict_count_ += sm_.iso_conflict_count_;
        si_empty_count_ += sm_.si_empty_count_;
        lc_empty_count_ += sm_.lc_empty_count_;

        phase3_sum_block_reuse_entries_ += sm_.phase3_block_reuse_entry_count_;
        phase3_sum_block_reuse_fallback_plain_ += sm_.phase3_block_reuse_fallback_plain_count_;
        phase3_sum_block_reuse_skeleton_runs_ += sm_.phase3_block_reuse_skeleton_run_count_;
        phase3_sum_skeleton_invalid_depth_breaks_ += sm_.phase3_skeleton_invalid_depth_break_count_;
        phase3_sum_cut_vertex_branches_ += sm_.phase3_cut_vertex_branch_count_;
        phase3_sum_cut_skip_no_phase3_block_ += sm_.phase3_cut_skip_no_phase3_block_count_;
        phase3_sum_block_inner_dfs_runs_ += sm_.phase3_block_inner_dfs_run_count_;
        phase3_sum_block_cache_hits_ += sm_.phase3_block_cache_hit_count_;
        phase3_sum_block_cache_misses_ += sm_.phase3_block_cache_miss_count_;
        phase3_sum_block_cache_stores_ += sm_.phase3_block_cache_store_count_;

        sm_.reset_performance_counters();
    }
    return result_count;
}

uint64_t StreamingEngine::execute_phase1_on(const Graph *query_graph, const Update &update, bool enable_local_view, // Phase1：处理单条流更新；本参数为是否物化局部视图（LocalView）
                                            bool enable_search) { 
    LabelTriple label_triple = update.labels_; 
    Edge data_edge = update.edge_;
    const LabelTriple label_triple_initial = label_triple; 
    uint64_t result_count = 0; 
    is_relevant_ = false;
    is_searched_ = false; 
    if (update.op_ == '+') { 
        auto it = gvm_.get_mapped_views(label_triple); 
        if (it != nullptr) { 
            is_relevant_ = true; 
            relevant_update_count_ += 1;
            auto t_gvm_begin = std::chrono::high_resolution_clock::now();
            for (uint32_t i = 0; i < 2; ++i) { 
                it = gvm_.get_mapped_views(label_triple); 
                uint32_t view_id = it->second;
                gvm_.update_view(update.op_, view_id, data_edge, label_triple); 
                std::swap(data_edge.first, data_edge.second);
                std::swap(label_triple.src_label_, label_triple.dst_label_); 
            }
            std::swap(data_edge.first, data_edge.second); 
            std::swap(label_triple.src_label_, label_triple.dst_label_);
            auto t_gvm_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_views_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_gvm_end - t_gvm_begin).count();
        }
        if (is_relevant_) { 
            auto t_nlf_begin = std::chrono::high_resolution_clock::now();
            gvm_.update_nlf_view(update.op_, data_edge, label_triple); 
            auto t_nlf_end = std::chrono::high_resolution_clock::now();
            time_gvm_update_nlf_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_nlf_end - t_nlf_begin).count();
        }
   
        std::vector<Edge> *mapped_query_edges = om_.get_mapped_edges(label_triple);
        const size_t phase1_matching_query_edge_count = mapped_query_edges ? mapped_query_edges->size() : 0;
        auto mapped_automorphism = om_.get_mapped_automorphism(label_triple); 
        if (mapped_automorphism != nullptr) {
            is_relevant_ = true;
            if (enable_local_view) { 
                if (phase1_matching_query_edge_count <= 1) { 
                    
                    for (auto automorphism_id: *mapped_automorphism) { 
                        auto meta = om_.get_automorphism_meta(automorphism_id); 
                        OrdersPerEdge *orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta); 
                        auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                        bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge);
                        auto t_lvm_end = std::chrono::high_resolution_clock::now();
                        time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_end - t_lvm_begin).count();
                        if (is_valid) { 
                            search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_; 
                            search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_; 
                        } else { // 视图无效或提前剪枝
                            non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_; 
                            non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_; 

                            if (lvm_.generate_visited_neighbor_count_ == 0) 
                                direct_rejection_count_ += 1; 
                        }
                        first_indexing_vertex_ += lvm_.first_vertex_neighbor_; 
                        if (enable_search && is_valid) {
                            auto t_plain_begin = std::chrono::high_resolution_clock::now();
                            uint64_t local_result_count = sm_.search_on_reduced_query(query_graph, *orders, lvm_, gvm_); 
                            auto t_plain_end = std::chrono::high_resolution_clock::now();
                            time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_plain_end - t_plain_begin).count();
                            result_count += local_result_count * automorphism_size; 
                            is_searched_ = true;
                        }

                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view(); 
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(t_lvm_destroy_end - t_lvm_destroy_begin).count();

                        if (result_count >= sm_.target_number) break; 
                        if (g_exit) break; 
                    }
                }else {
                    struct PreGroupTask {
                        OrdersPerEdge* orders = nullptr;
                        uint32_t automorphism_size = 1;
                        uint32_t automorphism_id = 0;
                    };
                    struct PreGroupBucket {
                        uint32_t group_id = 0;
                        std::vector<PreGroupTask> tasks;
                        uint32_t chosen_sub_id = std::numeric_limits<uint32_t>::max();
                    };
                    std::vector<PreGroupTask> plain_search_pre;
                    std::unordered_map<uint32_t, PreGroupBucket> grouped_pre;
                    plain_search_pre.reserve(mapped_automorphism->size());
                    grouped_pre.reserve(32);
                    for (auto automorphism_id : *mapped_automorphism) {
                        auto meta = om_.get_automorphism_meta(automorphism_id);
                        OrdersPerEdge* orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta);
                        if (orders == nullptr) continue;
                        const uint32_t u0 = orders->indexing_order_[0];
                        const uint32_t u1 = orders->indexing_order_[1];
                        if (!gvm_.nlf_check(u0, data_edge.first) || !gvm_.nlf_check(u1, data_edge.second)) {
                            continue;
                        }
                        DirectedEdgeKey key{u0, u1, 0};
                        uint32_t gid = 0;
                        auto it_gid = phase1_edge_group_id_.find(key);
                        if (it_gid != phase1_edge_group_id_.end()) {
                             gid = it_gid->second;
                        }
                        if (gid == 0) {
                            plain_search_pre.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                            continue;
                        }
                        auto& b = grouped_pre[gid];
                        b.group_id = gid;
                        b.tasks.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                        b.chosen_sub_id = phase1_group_chosen_sub_id_[gid]; 
                    }

                    const auto t_plain_total_begin = std::chrono::high_resolution_clock::now();
                    for (const auto &pt : plain_search_pre) {
                        if (pt.orders == nullptr) continue;
                        auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                        bool is_valid = lvm_.create_view(query_graph, *pt.orders, gvm_, data_edge);
                        auto t_lvm_end = std::chrono::high_resolution_clock::now();
                        time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                t_lvm_end - t_lvm_begin).count();
                        if (enable_search && is_valid) {
                            auto t_plain_begin = std::chrono::high_resolution_clock::now();
                            uint64_t local_result_count =
                                    sm_.search_on_reduced_query(query_graph, *pt.orders, lvm_, gvm_);
                            auto t_plain_end = std::chrono::high_resolution_clock::now();
                            time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                    t_plain_end - t_plain_begin).count();
                            result_count += local_result_count * pt.automorphism_size;
                            is_searched_ = true;
                        }
                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view();
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                t_lvm_destroy_end - t_lvm_destroy_begin).count();
                        if (result_count >= sm_.target_number) break;
                        if (g_exit) break;
                    }
                    const auto t_plain_total_end = std::chrono::high_resolution_clock::now();
                    time_phase1on_plain_total_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            t_plain_total_end - t_plain_total_begin).count();
                      
                    for (auto &kv : grouped_pre) {
                        auto &bucket = kv.second;
                        const uint32_t group_sub_id = bucket.chosen_sub_id;
                        auto *subctx = get_subquery_context(group_sub_id);
                        if (subctx == nullptr) {
                            continue;
                        }
                        OrdersPerEdge* group_seed_orders =
                                std::get<1>(om_.get_automorphism_meta_special(bucket.tasks[0].automorphism_id, group_sub_id));
                        if (group_seed_orders == nullptr) {
                            group_seed_orders = bucket.tasks[0].orders;
                        }
                        if (group_seed_orders == nullptr) {
                            continue;
                        }
                        const uint32_t u0 = group_seed_orders->indexing_order_[0];
                        const uint32_t u1 = group_seed_orders->indexing_order_[1];
                        const int32_t n0 = subctx->old2new[u0];
                        const int32_t n1 = subctx->old2new[u1];
                        const uint32_t ns = static_cast<uint32_t>(n0);
                        const uint32_t nt = static_cast<uint32_t>(n1);

                        const uint32_t sub_step0_key =
                                (u0 & 0xFFu) | ((u1 & 0xFFu) << 8) | ((ns & 0xFFu) << 16) | ((nt & 0xFFu) << 24);
                        const std::vector<std::vector<uint32_t>>* sub_mapped_new_orbit_peers = nullptr;
                        if (n0 >= 0 && n1 >= 0) {
                            auto it_orbit = subctx->mapped_new_orbit_peers_by_key.find(sub_step0_key);
                            if (it_orbit != subctx->mapped_new_orbit_peers_by_key.end())
                                sub_mapped_new_orbit_peers = &it_orbit->second;
                        }
                        OrdersPerEdge *sub_orders = subctx->om.get_orders({ns, nt});
                        if (sub_orders == nullptr) {
                            continue;
                        }
                        bool sub_valid = subctx->lvm.create_candidates_only_for_subgraph(subctx->sub_query.get(),
                                                                                               query_graph,
                                                                                               *sub_orders,
                                                                                               gvm_,
                                                                                               data_edge,
                                                                                               &subctx->new2old,
                                                                                               &subctx->sub_vertex_nlf,
                                                                                               sub_mapped_new_orbit_peers);

                        if (!sub_valid) {
                            continue;
                        }
                        phase1_before_sub_valid_check_count_ += 1;
                        LocalCandidatesSnapshot sub_candidates_snapshot =
                                subctx->lvm.extract_candidates_snapshot(subctx->sub_query.get());

                        LocalCandidatesSnapshot sub_candidates_snapshot_for_tasks = sub_candidates_snapshot;
                        subctx->lvm.restore_candidates_from(std::move(sub_candidates_snapshot));
                        std::vector<std::vector<uint32_t>> sub_embeddings_for_group;
                        bool sub_embeddings_ready = false;
                        static const std::vector<std::vector<uint32_t>> k_empty_mapped_new_by_old;
                        for (auto &pt : bucket.tasks) {
                            OrdersPerEdge* orders_for_group =
                                    std::get<1>(om_.get_automorphism_meta_special(pt.automorphism_id, group_sub_id));

                
                            const uint32_t task_u0 = orders_for_group->indexing_order_[0];
                            const uint32_t task_u1 = orders_for_group->indexing_order_[1];
                            const uint32_t task_ns = ns;
                            const uint32_t task_nt = nt;
                            const uint32_t key = (task_u0 & 0xFFu) | ((task_u1 & 0xFFu) << 8) |
                                                 ((task_ns & 0xFFu) << 16) | ((task_nt & 0xFFu) << 24);
                            auto it_mapped_new_set = subctx->mapped_new_set_by_key_and_old.find(key);

                            const std::vector<std::vector<uint32_t>> *mapped_for_task =
                                    (it_mapped_new_set != subctx->mapped_new_set_by_key_and_old.end())
                                            ? &it_mapped_new_set->second
                                            : &k_empty_mapped_new_by_old;
                          
                            bool task_pruned_valid = lvm_.create_candidates_only_with_sub_constraint(
                                    query_graph,
                                    *orders_for_group,
                                    gvm_,
                                    data_edge,
                                    sub_candidates_snapshot_for_tasks,
                                    *mapped_for_task,
                                    subctx->sub_query.get(),
                                    &subctx->new2old);


                            phase1_before_task_pruned_check_count_ += 1;
                            if (!task_pruned_valid) {
                                continue;
                            }

                            if (!sub_embeddings_ready) {
                                subctx->lvm.build_view_from_candidates_for_subgraph(
                                        subctx->sub_query.get(), *sub_orders, gvm_, &subctx->new2old);
                                sub_embeddings_for_group = subctx->sm.special_search_on_reduced_query(subctx->sub_query.get(),
                                                                                                           *sub_orders,
                                                                                                           subctx->lvm,
                                                                                                           gvm_,
                                                                                                           subctx->sm.target_number,
                                                                                                           &subctx->new2old);
                                subctx->lvm.destroy_view();
                                sub_embeddings_ready = true;
                                
                            }

                            phase1_before_sub_embeddings_check_count_ += 1;
                            if (sub_embeddings_for_group.empty()) {
                                continue;
                            }
                            phase1_after_sub_embeddings_check_count_ += 1;

                                
                            auto it_mapping_ids = subctx->mapping_ids_index.find(key);
                            const std::vector<uint32_t>& candidate_mapping_ids = it_mapping_ids->second;

                            lvm_.build_view_from_candidates(query_graph, *orders_for_group, gvm_);
                            for (const auto& emb : sub_embeddings_for_group) {
                                for (uint32_t selected_mapping_id : candidate_mapping_ids) {
                                    auto it_tpl = subctx->forced_pairs_template.find(selected_mapping_id);

                                    uint64_t local_result_count =
                                            sm_.constrained_search_with_subemb(query_graph,
                                                                                *orders_for_group,
                                                                                lvm_,
                                                                                gvm_,
                                                                                it_tpl->second,
                                                                                emb);
                                    result_count += local_result_count * pt.automorphism_size;
                                    is_searched_ = true;
                                    if (result_count >= sm_.target_number) break;
                                    if (g_exit) break;
                                }
                                if (result_count >= sm_.target_number) break;
                                if (g_exit) break;
                            }
                            lvm_.destroy_view();
                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;

                        }
                    }
                }
            }
        }
    }
    else { 

        label_triple = label_triple_initial; 
        data_edge = update.edge_;

        std::vector<Edge> *mapped_query_edges_del = om_.get_mapped_edges(label_triple);
        const size_t phase1_matching_query_edge_count_del = mapped_query_edges_del ? mapped_query_edges_del->size() : 0;





        auto mapped_automorphism = om_.get_mapped_automorphism(label_triple); 
        if (mapped_automorphism != nullptr) { 
            is_relevant_ = true; 
            if (enable_local_view) { 
                if (phase1_matching_query_edge_count_del <= 1) { 
                    for (auto automorphism_id: *mapped_automorphism) { 
                        auto meta = om_.get_automorphism_meta(automorphism_id); 
                        OrdersPerEdge *orders = std::get<1>(meta); 
                        uint32_t automorphism_size = std::get<2>(meta);
                        bool is_valid = lvm_.create_view(query_graph, *orders, gvm_, data_edge);
                        if (is_valid) { 
                            search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_;
                            search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_; 
                        } else { 
                            non_search_generate_neighbor_count_ += lvm_.generate_visited_neighbor_count_; 
                            non_search_build_neighbor_count_ += lvm_.build_visited_neighbor_count_; 

                            if (lvm_.generate_visited_neighbor_count_ == 0) 
                                direct_rejection_count_ += 1; 
                        }
                        first_indexing_vertex_ += lvm_.first_vertex_neighbor_; 
                        if (enable_search && is_valid) { 
                            uint64_t local_result_count = sm_.search_on_reduced_query(query_graph, *orders, lvm_, gvm_); // 枚举当前状态下的匹配（删边前）
                            result_count += local_result_count * automorphism_size; 
                            is_searched_ = true;
                        }
                        lvm_.destroy_view(); 

                        if (result_count >= sm_.target_number) break; 
                        if (g_exit) break; 
                    }
                } else {
                     struct PreGroupTask {
                        OrdersPerEdge* orders = nullptr;
                        uint32_t automorphism_size = 1;
                        uint32_t automorphism_id = 0;
                    };
                    struct PreGroupBucket {
                        uint32_t group_id = 0;
                        std::vector<PreGroupTask> tasks;
                        uint32_t chosen_sub_id = std::numeric_limits<uint32_t>::max();
                    };
                    std::vector<PreGroupTask> plain_search_pre;
                    std::unordered_map<uint32_t, PreGroupBucket> grouped_pre;
                    plain_search_pre.reserve(mapped_automorphism->size());
                    grouped_pre.reserve(32);
                    for (auto automorphism_id : *mapped_automorphism) {
                        auto meta = om_.get_automorphism_meta(automorphism_id);
                        OrdersPerEdge* orders = std::get<1>(meta);
                        uint32_t automorphism_size = std::get<2>(meta);
                        if (orders == nullptr) continue;
                        const uint32_t u0 = orders->indexing_order_[0];
                        const uint32_t u1 = orders->indexing_order_[1];
                        if (!gvm_.nlf_check(u0, data_edge.first) || !gvm_.nlf_check(u1, data_edge.second)) {
                            continue;
                        }
                        DirectedEdgeKey key{u0, u1, 0};
                        uint32_t gid = 0;
                        auto it_gid = phase1_edge_group_id_.find(key);
                        if (it_gid != phase1_edge_group_id_.end()) {
                             gid = it_gid->second;
                        }
                        if (gid == 0) {
                            plain_search_pre.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                            continue;
                        }
                        auto& b = grouped_pre[gid];
                        b.group_id = gid;
                        b.tasks.push_back(PreGroupTask{orders, automorphism_size, automorphism_id});
                        b.chosen_sub_id = phase1_group_chosen_sub_id_[gid]; 
                    }

                    const auto t_plain_total_begin = std::chrono::high_resolution_clock::now();
                    for (const auto &pt : plain_search_pre) {
                        if (pt.orders == nullptr) continue;
                        auto t_lvm_begin = std::chrono::high_resolution_clock::now();
                        bool is_valid = lvm_.create_view(query_graph, *pt.orders, gvm_, data_edge);
                        auto t_lvm_end = std::chrono::high_resolution_clock::now();
                        time_lvm_create_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                t_lvm_end - t_lvm_begin).count();
                        if (enable_search && is_valid) {
                            auto t_plain_begin = std::chrono::high_resolution_clock::now();
                            uint64_t local_result_count =
                                    sm_.search_on_reduced_query(query_graph, *pt.orders, lvm_, gvm_);
                            auto t_plain_end = std::chrono::high_resolution_clock::now();
                            time_plain_search_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                    t_plain_end - t_plain_begin).count();
                            result_count += local_result_count * pt.automorphism_size;
                            is_searched_ = true;
                        }
                        auto t_lvm_destroy_begin = std::chrono::high_resolution_clock::now();
                        lvm_.destroy_view();
                        auto t_lvm_destroy_end = std::chrono::high_resolution_clock::now();
                        time_lvm_destroy_view_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                t_lvm_destroy_end - t_lvm_destroy_begin).count();
                        if (result_count >= sm_.target_number) break;
                        if (g_exit) break;
                    }
                    const auto t_plain_total_end = std::chrono::high_resolution_clock::now();
                    time_phase1on_plain_total_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            t_plain_total_end - t_plain_total_begin).count();
                      
                    for (auto &kv : grouped_pre) {
                        auto &bucket = kv.second;
                        const uint32_t group_sub_id = bucket.chosen_sub_id;
                        auto *subctx = get_subquery_context(group_sub_id);
                        if (subctx == nullptr) {
                            continue;
                        }
                        OrdersPerEdge* group_seed_orders =
                                std::get<1>(om_.get_automorphism_meta_special(bucket.tasks[0].automorphism_id, group_sub_id));
                        if (group_seed_orders == nullptr) {
                            group_seed_orders = bucket.tasks[0].orders;
                        }
                        if (group_seed_orders == nullptr) {
                            continue;
                        }
                        const uint32_t u0 = group_seed_orders->indexing_order_[0];
                        const uint32_t u1 = group_seed_orders->indexing_order_[1];
                        const int32_t n0 = subctx->old2new[u0];
                        const int32_t n1 = subctx->old2new[u1];
                        const uint32_t ns = static_cast<uint32_t>(n0);
                        const uint32_t nt = static_cast<uint32_t>(n1);

                        const uint32_t sub_step0_key =
                                (u0 & 0xFFu) | ((u1 & 0xFFu) << 8) | ((ns & 0xFFu) << 16) | ((nt & 0xFFu) << 24);
                        const std::vector<std::vector<uint32_t>>* sub_mapped_new_orbit_peers = nullptr;
                        if (n0 >= 0 && n1 >= 0) {
                            auto it_orbit = subctx->mapped_new_orbit_peers_by_key.find(sub_step0_key);
                            if (it_orbit != subctx->mapped_new_orbit_peers_by_key.end())
                                sub_mapped_new_orbit_peers = &it_orbit->second;
                        }
                        OrdersPerEdge *sub_orders = subctx->om.get_orders({ns, nt});
                        if (sub_orders == nullptr) {
                            continue;
                        }
                        bool sub_valid = subctx->lvm.create_candidates_only_for_subgraph(subctx->sub_query.get(),
                                                                                               query_graph,
                                                                                               *sub_orders,
                                                                                               gvm_,
                                                                                               data_edge,
                                                                                               &subctx->new2old,
                                                                                               &subctx->sub_vertex_nlf,
                                                                                               sub_mapped_new_orbit_peers);

                        if (!sub_valid) {
                            continue;
                        }
                        phase1_before_sub_valid_check_count_ += 1;
                        LocalCandidatesSnapshot sub_candidates_snapshot =
                                subctx->lvm.extract_candidates_snapshot(subctx->sub_query.get());

                        LocalCandidatesSnapshot sub_candidates_snapshot_for_tasks = sub_candidates_snapshot;
                        subctx->lvm.restore_candidates_from(std::move(sub_candidates_snapshot));
                        std::vector<std::vector<uint32_t>> sub_embeddings_for_group;
                        bool sub_embeddings_ready = false;
                        static const std::vector<std::vector<uint32_t>> k_empty_mapped_new_by_old;
                        for (auto &pt : bucket.tasks) {
                            OrdersPerEdge* orders_for_group =
                                    std::get<1>(om_.get_automorphism_meta_special(pt.automorphism_id, group_sub_id));

                
                            const uint32_t task_u0 = orders_for_group->indexing_order_[0];
                            const uint32_t task_u1 = orders_for_group->indexing_order_[1];
                            const uint32_t task_ns = ns;
                            const uint32_t task_nt = nt;
                            const uint32_t key = (task_u0 & 0xFFu) | ((task_u1 & 0xFFu) << 8) |
                                                 ((task_ns & 0xFFu) << 16) | ((task_nt & 0xFFu) << 24);
                            auto it_mapped_new_set = subctx->mapped_new_set_by_key_and_old.find(key);

                            const std::vector<std::vector<uint32_t>> *mapped_for_task =
                                    (it_mapped_new_set != subctx->mapped_new_set_by_key_and_old.end())
                                            ? &it_mapped_new_set->second
                                            : &k_empty_mapped_new_by_old;
                          
                            bool task_pruned_valid = lvm_.create_candidates_only_with_sub_constraint(
                                    query_graph,
                                    *orders_for_group,
                                    gvm_,
                                    data_edge,
                                    sub_candidates_snapshot_for_tasks,
                                    *mapped_for_task,
                                    subctx->sub_query.get(),
                                    &subctx->new2old);


                            phase1_before_task_pruned_check_count_ += 1;
                            if (!task_pruned_valid) {
                                continue;
                            }

                            if (!sub_embeddings_ready) {
                                subctx->lvm.build_view_from_candidates_for_subgraph(
                                        subctx->sub_query.get(), *sub_orders, gvm_, &subctx->new2old);
                                sub_embeddings_for_group = subctx->sm.special_search_on_reduced_query(subctx->sub_query.get(),
                                                                                                           *sub_orders,
                                                                                                           subctx->lvm,
                                                                                                           gvm_,
                                                                                                           subctx->sm.target_number,
                                                                                                           &subctx->new2old);
                                subctx->lvm.destroy_view();
                                sub_embeddings_ready = true;
                                
                            }

                            phase1_before_sub_embeddings_check_count_ += 1;
                            if (sub_embeddings_for_group.empty()) {
                                continue;
                            }
                            phase1_after_sub_embeddings_check_count_ += 1;

                                
                            auto it_mapping_ids = subctx->mapping_ids_index.find(key);
                            const std::vector<uint32_t>& candidate_mapping_ids = it_mapping_ids->second;

                            lvm_.build_view_from_candidates(query_graph, *orders_for_group, gvm_);
                            for (const auto& emb : sub_embeddings_for_group) {
                                for (uint32_t selected_mapping_id : candidate_mapping_ids) {
                                    auto it_tpl = subctx->forced_pairs_template.find(selected_mapping_id);

                                    uint64_t local_result_count =
                                            sm_.constrained_search_with_subemb(query_graph,
                                                                                *orders_for_group,
                                                                                lvm_,
                                                                                gvm_,
                                                                                it_tpl->second,
                                                                                emb);
                                    result_count += local_result_count * pt.automorphism_size;
                                    is_searched_ = true;
                                    if (result_count >= sm_.target_number) break;
                                    if (g_exit) break;
                                }
                                if (result_count >= sm_.target_number) break;
                                if (g_exit) break;
                            }
                            lvm_.destroy_view();
                            if (result_count >= sm_.target_number) break;
                            if (g_exit) break;

                        }
                    }
                }
            }
        }

        auto it = gvm_.get_mapped_views(label_triple); 
        if (it != nullptr) { 
            is_relevant_ = true;
            relevant_update_count_ += 1; 
            for (uint32_t i = 0; i < 2; ++i) { 
                it = gvm_.get_mapped_views(label_triple);
                uint32_t view_id = it->second;
                gvm_.update_view(update.op_, view_id, data_edge, label_triple); 
                std::swap(data_edge.first, data_edge.second);
                std::swap(label_triple.src_label_, label_triple.dst_label_); 
            }

            std::swap(data_edge.first, data_edge.second); 
            std::swap(label_triple.src_label_, label_triple.dst_label_); 
        }

        if (is_relevant_) { 
            gvm_.update_nlf_view(update.op_, data_edge, label_triple); 
        }
    }
    edge_process_count_ += 1;
    if (is_relevant_) { 
        if (result_count > 0) { 
            positive_count_ += 1; 
        }
        if (is_searched_) { 
            search_count_ += 1; 
        }

        result_count_ += result_count; 
        invalid_partial_result_count_ += sm_.invalid_partial_result_count_; 
        partial_result_count_ += sm_.partial_result_count_; 
        iso_conflict_count_ += sm_.iso_conflict_count_; 
        si_empty_count_ += sm_.si_empty_count_;
        lc_empty_count_ += sm_.lc_empty_count_; 
        sm_.reset_performance_counters(); 
    }
    return result_count;
}

void StreamingEngine::evaluate_view_update(const Graph *query_graph, const Graph *data_graph,
                                           const std::vector<Update> &stream) {
    const uint32_t repetition_time = 1;
    {
        global_view_update_time_ = 0;
        gvm_.release();
    
        for (uint32_t i = 0; i < repetition_time; ++i) {
            gvm_.initialize(query_graph);
            gvm_.create_views(query_graph, data_graph);

            auto start = std::chrono::high_resolution_clock::now();

#ifdef MEASURE_INDEXING_COST
            uint32_t local_update_count = 0;
#endif

            for (auto& update : stream) {
                execute(query_graph, update, false, false);

#ifdef MEASURE_INDEXING_COST
                local_update_count += 1;
                if (local_update_count >= g_update_count) break;
#endif
            }

            auto end = std::chrono::high_resolution_clock::now();
            global_view_update_time_ += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

            gvm_.release();
        }

        global_view_update_time_ /= repetition_time;
        printf("Global view update time (seconds): %.6f\n", NANOSECTOSEC( global_view_update_time_));
    }
    {
        local_view_update_time_ = 0;

        search_build_neighbor_count_ = 0;
        non_search_build_neighbor_count_ = 0;
        search_generate_neighbor_count_ = 0;
        non_search_generate_neighbor_count_ = 0;
        direct_rejection_count_ = 0;
        first_indexing_vertex_ = 0;
        for (uint32_t i = 0; i < repetition_time; ++i) {
            gvm_.initialize(query_graph);
            gvm_.create_views(query_graph, data_graph);

#ifdef MEASURE_INDEXING_COST
            uint32_t local_update_count = 0;
#endif

            auto start = std::chrono::high_resolution_clock::now();
            for (auto& update : stream) {
                execute(query_graph, update, true, false);

#ifdef MEASURE_INDEXING_COST
                local_update_count += 1;
                if (local_update_count >= g_update_count) break;
#endif
            }

            auto end = std::chrono::high_resolution_clock::now();
            local_view_update_time_ += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

            gvm_.release();
        }

        local_view_update_time_ /= repetition_time;
        local_view_update_time_ = local_view_update_time_ > global_view_update_time_ ?
                local_view_update_time_ - global_view_update_time_ : 0;
        printf("Local view update time (seconds): %.6f\n", NANOSECTOSEC(local_view_update_time_));
    }

    printf("View update time (seconds): %.6f\n", NANOSECTOSEC(global_view_update_time_ + local_view_update_time_));
}

void
StreamingEngine::evaluate_search(const Graph *query_graph, const Graph *data_graph, const std::vector<Update> &stream) {
    gvm_.release();
    gvm_.initialize(query_graph);
    gvm_.create_views(query_graph, data_graph);

    reset_performance_counters();
    auto start = std::chrono::high_resolution_clock::now();

#ifdef MEASURE_UPDATE_COST
    uint32_t count = 0;
    auto measure_begin = std::chrono::high_resolution_clock::now();
    auto measure_end = std::chrono::high_resolution_clock::now();
    bool measure_flip = true;
#endif
#ifdef MEASURE_INDEXING_COST
    g_update_count = 0;
#endif


    for (auto& update : stream) {
        execute(query_graph, update, true, true);

#ifdef MEASURE_INDEXING_COST
        g_update_count += 1;
#endif
#ifdef MEASURE_UPDATE_COST
        count += 1;
        if (count % MEASURE_BATCH_SIZE == 0) {
            uint64_t update_cost = 0;
            if (measure_flip) {
                measure_end = std::chrono::high_resolution_clock::now();
                update_cost = std::chrono::duration_cast<std::chrono::nanoseconds>(measure_end - measure_begin).count();
            }
            else {
                measure_begin = std::chrono::high_resolution_clock::now();
                update_cost = std::chrono::duration_cast<std::chrono::nanoseconds>(measure_begin - measure_end).count();
            }
            measure_flip = !measure_flip;
            g_measure_time_cost_per_update.emplace_back(update_cost);
        }
#endif

        if (g_exit) {

#ifdef MEASURE_UPDATE_COST
            if (count % MEASURE_BATCH_SIZE != 0) {
                uint64_t update_cost = 0;
                if (measure_flip) {
                    measure_end = std::chrono::high_resolution_clock::now();
                    update_cost = std::chrono::duration_cast<std::chrono::nanoseconds>(measure_end - measure_begin).count();
                }
                else {
                    measure_begin = std::chrono::high_resolution_clock::now();
                    update_cost = std::chrono::duration_cast<std::chrono::nanoseconds>(measure_begin - measure_end).count();
                }
                g_measure_time_cost_per_update.emplace_back(update_cost);
            }
#endif
            break;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    query_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
}

std::string getFileNameWithoutExtension(const std::string& s);

void StreamingEngine::print_metrics() {
    uint64_t view_update_time = global_view_update_time_ + local_view_update_time_;
    uint64_t search_time = query_time_ > view_update_time ? query_time_ - view_update_time : 0;
    printf("%.6f\n", NANOSECTOSEC(search_time));
}

void load_stream(const std::string& file_path, std::vector<Update>& stream, const Graph* data_graph) {
    auto start = std::chrono::high_resolution_clock::now();

    uint32_t vertex_num = data_graph->getVerticesCount();
    spp::sparse_hash_map<uint32_t, uint32_t> new_vertex_label;
    std::vector<uint32_t> vertex_label;
    vertex_label.reserve(vertex_num);
    for (uint32_t u = 0; u < vertex_num; ++u) {
        vertex_label.push_back(data_graph->getVertexLabel(u));
    }

    std::ifstream ifs(file_path);

    if (!ifs.is_open()) {
        std::cout << "Can not open the stream file " << file_path << " ." << std::endl;
        exit(-1);
    }

    Update update;
    std::string tmp_str;
    // 用 getline 作循环条件，避免 while(good)+getline 在 EOF/空行时多跑一次；空行上 tmp_str[0] 为 UB，可表现为间歇性崩溃。
    while (std::getline(ifs, tmp_str)) {
        if (tmp_str.empty() || tmp_str[0] == '#')
            continue;

        std::stringstream ss(tmp_str);
        std::string op_str;
        ss >> op_str;
        //以上改过

        if (op_str == "v") {
            uint32_t id;
            uint32_t label;
            ss >> id >> label;
            if (id < vertex_num) {
                vertex_label[id] = label;
            }
            else {
                new_vertex_label[id] = label;
            }
        }
        else if (op_str == "e" || op_str == "-e") {
            update.op_ = op_str == "e" ? '+' : '-';
            ss >> update.edge_.first >> update.edge_.second >> update.labels_.edge_label_;

            if (update.edge_.first < vertex_num) {
                update.labels_.src_label_ = vertex_label[update.edge_.first];
            }
            else {
                update.labels_.src_label_ = new_vertex_label[update.edge_.first];
            }

            if (update.edge_.second < vertex_num) {
                update.labels_.dst_label_ = vertex_label[update.edge_.second];
            }
            else {
                update.labels_.dst_label_ = new_vertex_label[update.edge_.second];
            }

            stream.emplace_back(update);
        }
    }

    ifs.close();

    auto end = std::chrono::high_resolution_clock::now();
    auto load_stream_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
}

bool g_exit = false;

void execute_within_time_limit(const Graph *query_graph, const Graph *data_graph, std::vector<Update> &stream,
                               StreamingEngine *sm, uint64_t time_limit) {
    g_exit = false;
    auto future = std::async(std::launch::async, [query_graph, data_graph, &stream, sm](){
        sm->evaluate_search(query_graph, data_graph, stream);
    });

    std::future_status status;
    do {
        status = future.wait_for(std::chrono::seconds(time_limit));
        if (status == std::future_status::deferred) {
            std::cout << "Deferred\n";
            exit(-1);
        } else if (status == std::future_status::timeout) {
            g_exit = true;
        }
    } while (status != std::future_status::ready);
}

std::string getFileNameWithoutExtension(const std::string& s) {
    char sep = '/';
#ifdef _WIN32
    sep = '\\';
#endif
    size_t i = s.rfind(sep, s.length());
    std::string filename = s;
    if (i != std::string::npos) {
        filename = s.substr(i + 1, s.length() - i);
    }
    size_t lastindex = filename.find_last_of('.');
    if (lastindex != std::string::npos) {
        return filename.substr(0, lastindex);
    }
    return filename;
}

int main(int argc, char** argv) {
#ifdef _WIN32
    static bool crash_handler_installed = false;
    if (!crash_handler_installed) {
        crash_handler_installed = true;
        SetUnhandledExceptionFilter([](EXCEPTION_POINTERS* ep) -> LONG {
            (void)ep;
            HANDLE process = GetCurrentProcess();
            HMODULE self = GetModuleHandleW(nullptr);
            const DWORD64 base = (DWORD64)self;
            SymSetOptions(SYMOPT_DEFERRED_LOADS | SYMOPT_UNDNAME);
            SymInitialize(process, nullptr, TRUE);

            void* stack[64];
            USHORT frames = CaptureStackBackTrace(0, 64, stack, nullptr);

            fprintf(stderr, "\n=== Unhandled exception: stack trace (%u frames) ===\n", (unsigned)frames);
            fprintf(stderr, "Module base: 0x%llx\n", (unsigned long long)base);
            SYMBOL_INFO* sym = (SYMBOL_INFO*)calloc(sizeof(SYMBOL_INFO) + 256, 1);
            if (sym) {
                sym->MaxNameLen = 255;
                sym->SizeOfStruct = sizeof(SYMBOL_INFO);
            }
            for (USHORT i = 0; i < frames; ++i) {
                DWORD64 addr = (DWORD64)(stack[i]);
                DWORD64 disp = 0;
                const unsigned long long off = (addr >= base) ? (unsigned long long)(addr - base) : 0ull;
                if (sym && SymFromAddr(process, addr, &disp, sym)) {
                    fprintf(stderr, "#%u  %s + 0x%llx  [0x%llx, +0x%llx]\n",
                            (unsigned)i, sym->Name,
                            (unsigned long long)disp, (unsigned long long)addr, off);
                } else {
                    fprintf(stderr, "#%u  [0x%llx, +0x%llx]\n",
                            (unsigned)i, (unsigned long long)addr, off);
                }
            }
            if (sym) free(sym);
            fprintf(stderr, "=== End stack trace ===\n");
            fflush(stderr);
            return EXCEPTION_EXECUTE_HANDLER;
        });
    }
#endif
    InputParser cmd_parser(argc, argv);
    std::string input_data_graph_file = cmd_parser.get_cmd_option("-d");
    std::string input_data_graph_update_file = cmd_parser.get_cmd_option("-u");
    std::string input_query_graph_file = cmd_parser.get_cmd_option("-q");
    std::string input_target_embedding_number = cmd_parser.get_cmd_option("-num");
    std::string input_step_length = cmd_parser.get_cmd_option("-s");
    std::string input_time_limit = cmd_parser.get_cmd_option("-time_limit");
    std::string input_phase1 = cmd_parser.get_cmd_option("-phase1");
    std::string input_phase2 = cmd_parser.get_cmd_option("-phase2");
    std::string input_phase3 = cmd_parser.get_cmd_option("-phase3");
    std::string input_phase4 = cmd_parser.get_cmd_option("-phase4");
    std::string input_phase5 = cmd_parser.get_cmd_option("-phase5");

    uint64_t target_embedding_number = std::numeric_limits<uint64_t>::max();
    if (!input_target_embedding_number.empty()) {
        if (input_target_embedding_number != "MAX") {
            target_embedding_number = std::stoll(input_target_embedding_number);
        }
    }

    uint32_t step_length = 10000;
    if (!input_step_length.empty()) {
        step_length = std::stoul(input_step_length);
    }

    uint64_t time_limit = 3600;
    if (!input_time_limit.empty()) {
        time_limit = std::stoll(input_time_limit);
    }

    bool phase1_on = false;
    if (!input_phase1.empty()) {
        std::string p = input_phase1;
        for (auto& c : p) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (p == "off" || p == "0" || p == "false" || p == "no") {
            phase1_on = false;
        } else if (p == "on" || p == "1" || p == "true" || p == "yes") {
            phase1_on = true;
        } else {
            std::cerr << "Invalid -phase1 value (use on or off). Got: " << input_phase1 << std::endl;
            return 1;
        }
    }

    bool phase3_on = false;
    if (!input_phase3.empty()) {
        std::string p = input_phase3;
        for (auto& c : p) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (p == "off" || p == "0" || p == "false" || p == "no") {
            phase3_on = false;
        } else if (p == "on" || p == "1" || p == "true" || p == "yes") {
            phase3_on = true;
        } else {
            std::cerr << "Invalid -phase3 value (use on or off). Got: " << input_phase3 << std::endl;
            return 1;
        }
    }

    bool phase2_on = false;
    if (!input_phase2.empty()) {
        std::string p = input_phase2;
        for (auto& c : p) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (p == "off" || p == "0" || p == "false" || p == "no") {
            phase2_on = false;
        } else if (p == "on" || p == "1" || p == "true" || p == "yes") {
            phase2_on = true;
        } else {
            std::cerr << "Invalid -phase2 value (use on or off). Got: " << input_phase2 << std::endl;
            return 1;
        }
    }

    bool phase4_on = false;
    if (!input_phase4.empty()) {
        std::string p = input_phase4;
        for (auto& c : p) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (p == "off" || p == "0" || p == "false" || p == "no") {
            phase4_on = false;
        } else if (p == "on" || p == "1" || p == "true" || p == "yes") {
            phase4_on = true;
        } else {
            std::cerr << "Invalid -phase4 value (use on or off). Got: " << input_phase4 << std::endl;
            return 1;
        }
    }

    bool phase5_on = false;
    if (!input_phase5.empty()) {
        std::string p = input_phase5;
        for (auto& c : p) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (p == "off" || p == "0" || p == "false" || p == "no") {
            phase5_on = false;
        } else if (p == "on" || p == "1" || p == "true" || p == "yes") {
            phase5_on = true;
        } else {
            std::cerr << "Invalid -phase5 value (use on or off). Got: " << input_phase5 << std::endl;
            return 1;
        }
    }


#ifdef MEASURE_UPDATE_COST
    g_measure_time_cost_per_update.reserve(1000000);
#endif

    Graph* data_graph = new Graph(false);
    data_graph->is_edge_labeled = true;
    data_graph->loadGraphFromFileWithoutMeta(input_data_graph_file);

    Graph* query_graph = new Graph(false);
    query_graph->is_edge_labeled = true;
    query_graph->loadGraphFromFileWithoutMeta(input_query_graph_file);

    std::vector<Update> stream;
    load_stream(input_data_graph_update_file, stream, data_graph);
    StreamingEngine streaming_engine;
    streaming_engine.set_phase1_enabled(phase1_on);
    streaming_engine.set_phase2_enabled(phase2_on);
    streaming_engine.set_phase3_enabled(phase3_on);
    streaming_engine.set_phase4_enabled(phase4_on);
    streaming_engine.set_phase4_enabled(phase5_on);
    streaming_engine.set_query_graph_file(input_query_graph_file);
    streaming_engine.set_data_graph_file(input_data_graph_file);
    streaming_engine.initialize(query_graph, data_graph, target_embedding_number);
    streaming_engine.preprocess(query_graph, data_graph);

    execute_within_time_limit(query_graph, data_graph, stream, &streaming_engine, time_limit);

#ifdef MEASURE_INDEXING_COST
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Evaluate the performance of view update...\n";
    streaming_engine.evaluate_view_update(query_graph, data_graph, stream);
#endif
    streaming_engine.print_metrics();
#ifdef MEASURE_UPDATE_COST
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Dump measure update results..." << std::endl;

    {
        auto file_name = getFileNameWithoutExtension(input_query_graph_file);
        std::string file_path = file_name + "_time_cost_update.bin";
        std::ofstream ofs(file_path, std::ios::binary);

        // Format: The first 8 byte records the number of elements, the second 4 byte records the batch size, then append the vector content.
        uint64_t num_element = g_measure_time_cost_per_update.size();
        uint32_t batch_size = MEASURE_BATCH_SIZE;
        int64_t vec_size = sizeof(uint64_t) * num_element;

        ofs.write(reinterpret_cast<const char *>(&num_element), 8);
        ofs.write(reinterpret_cast<const char *>(&batch_size), 4);
        ofs.write(reinterpret_cast<const char *>(&g_measure_time_cost_per_update.front()), vec_size);
    }

    std::cout << "--------------------------------------------------------------------" << std::endl;
#endif

    if (g_exit)
        std::cout << "Time out..." << std::endl;

    delete query_graph;
    delete data_graph;

    return 0;
}