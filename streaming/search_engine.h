#ifndef RAPIDMATCH_SEARCH_ENGINE_H
#define RAPIDMATCH_SEARCH_ENGINE_H
#include "local_view_manager.h"
#include <functional>
#include <unordered_map>
#include <cstdint>
#include <vector>
#include "query_bctree.h"
using spp::sparse_hash_set;
class SearchEngine {
public:
    uint64_t target_number = 1000;
public:
    uint64_t invalid_partial_result_count_;
    uint64_t partial_result_count_;
    uint64_t iso_conflict_count_;
    uint64_t si_empty_count_;
    uint64_t lc_empty_count_;
    uint64_t phase3_block_reuse_entry_count_ = 0;
    uint64_t phase3_block_reuse_fallback_plain_count_ = 0;
    uint64_t phase3_block_reuse_skeleton_run_count_ = 0;
    uint64_t phase3_skeleton_invalid_depth_break_count_ = 0;
    uint64_t phase3_cut_vertex_branch_count_ = 0;
    uint64_t phase3_cut_skip_no_phase3_block_count_ = 0;
    uint64_t phase3_block_inner_dfs_run_count_ = 0;
    uint64_t phase3_block_cache_hit_count_ = 0;
    uint64_t phase3_block_cache_miss_count_ = 0;
    uint64_t phase3_block_cache_store_count_ = 0;
    uint64_t phase3_block_cache_skip_no_entry_map_count_ = 0;
    void reset_performance_counters() {
        invalid_partial_result_count_ = 0;
        partial_result_count_ = 0;
        iso_conflict_count_ = 0;
        si_empty_count_ = 0;
        lc_empty_count_ = 0;
        phase3_block_reuse_entry_count_ = 0;
        phase3_block_reuse_fallback_plain_count_ = 0;
        phase3_block_reuse_skeleton_run_count_ = 0;
        phase3_skeleton_invalid_depth_break_count_ = 0;
        phase3_cut_vertex_branch_count_ = 0;
        phase3_cut_skip_no_phase3_block_count_ = 0;
        phase3_block_inner_dfs_run_count_ = 0;
        phase3_block_cache_hit_count_ = 0;
        phase3_block_cache_miss_count_ = 0;
        phase3_block_cache_store_count_ = 0;
        phase3_block_cache_skip_no_entry_map_count_ = 0;
    }
private:
    std::vector<uint32_t*> local_candidates_store_;
    std::vector<uint32_t*> encoded_local_candidates_store_;

    std::vector<std::vector<uint32_t>> local_candidates_buffer1_;
    std::vector<std::vector<uint32_t>> local_candidates_buffer2_;
    std::vector<std::pair<uint32_t, uint32_t>> local_idx_;
    bool* visited_;
    std::vector<uint32_t> embedding_;
    std::vector<uint32_t> encoded_embedding_;
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint64_t>> phase3_block_cache_;

private:
    uint32_t compute_local_candidates_for_reduced_query(const Graph *query_graph, uint32_t depth,
                                                        std::vector<uint32_t> &order,
                                                        std::vector<uint32_t> &bn_offset,
                                                        std::vector<uint32_t> &bn,
                                                        std::vector<uint32_t> &view_mapping,
                                                        LocalViewManager &lvm, GlobalViewManager &gvm);
    uint32_t compute_local_candidates_for_reduced_query_with_new2old(const Graph *query_graph, uint32_t depth,
                                                                    std::vector<uint32_t> &order,
                                                                    std::vector<uint32_t> &bn_offset,
                                                                    std::vector<uint32_t> &bn,
                                                                    std::vector<uint32_t> &view_mapping,
                                                                    LocalViewManager &lvm, GlobalViewManager &gvm,
                                                                    const std::vector<int32_t> &new2old);
public:
    SearchEngine() {}
    ~SearchEngine() {}

    void initialize(const Graph *query_graph, const Graph *data_graph);
    void release();

    uint64_t search_on_reduced_query(const Graph *query_graph, OrdersPerEdge &orders, LocalViewManager &lvm,
                                     GlobalViewManager &gvm);

    uint64_t block_reuse_search_on_reduced_query(const Graph *query_graph,
                                                 OrdersPerEdge &orders,
                                                 LocalViewManager &lvm,
                                                 GlobalViewManager &gvm,
                                                 const QueryBCTreeDecomposition& query_bctree);

    
    void clear_phase3_block_cache() { phase3_block_cache_.clear(); }
    std::vector<std::vector<uint32_t>> special_search_on_reduced_query(const Graph *query_graph,
                                                                       OrdersPerEdge &orders,
                                                                       LocalViewManager &lvm,
                                                                       GlobalViewManager &gvm,
                                                                       uint64_t max_results,
                                                                       const std::vector<int32_t>* new2old = nullptr);

    uint64_t constrained_search_on_reduced_query(const Graph *query_graph,
                                                 OrdersPerEdge &orders,
                                                 LocalViewManager &lvm,
                                                 GlobalViewManager &gvm,
                                                 const std::vector<std::pair<uint32_t, uint32_t>> &forced_assignments);
    uint64_t constrained_search_with_subemb(const Graph *query_graph,
                                            OrdersPerEdge &orders,
                                            LocalViewManager &lvm,
                                            GlobalViewManager &gvm,
                                            const std::vector<std::pair<uint32_t, uint32_t>> &forced_pairs_template,
                                            const std::vector<uint32_t> &emb);
};
#endif 
