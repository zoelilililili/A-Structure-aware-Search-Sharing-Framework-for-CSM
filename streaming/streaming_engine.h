#ifndef RAPIDMATCH_STREAMING_ENGINE_H
#define RAPIDMATCH_STREAMING_ENGINE_H
#include "order_manager.h"
#include "query_bctree.h"
#include "search_engine.h"
#include <unordered_map>
#include <utility>
#include <vector>
#include <string>
class StreamingEngine {
public:
    bool is_relevant_;
    bool is_searched_;
    uint64_t query_time_ = 0;
    uint64_t global_view_initialize_time_ = 0;
    uint64_t order_generation_time_ = 0;
    uint64_t global_view_update_time_ = 0;
    uint64_t local_view_update_time_ = 0;
    uint64_t edge_process_count_ = 0;
    uint64_t relevant_update_count_ = 0;
    uint64_t phase1off_mapped_automorphism_gt1_update_count_ = 0;
    uint64_t search_count_ = 0;
    uint64_t positive_count_ = 0;
    uint64_t result_count_ = 0;
    uint64_t invalid_partial_result_count_ = 0;
    uint64_t partial_result_count_ = 0;
    uint64_t iso_conflict_count_ = 0;
    uint64_t si_empty_count_ = 0;
    uint64_t lc_empty_count_ = 0;
    uint64_t non_search_generate_neighbor_count_;
    uint64_t search_generate_neighbor_count_;
    uint64_t non_search_build_neighbor_count_;
    uint64_t search_build_neighbor_count_;
    uint64_t direct_rejection_count_;
    uint64_t first_indexing_vertex_;
    uint64_t phase1_before_sub_valid_check_count_ = 0;
    uint64_t phase1_before_task_pruned_check_count_ = 0;
    uint64_t phase1_after_sub_embeddings_check_count_ = 0;
    uint64_t phase1_before_sub_embeddings_check_count_ = 0;
    uint64_t phase1_compare_normal_search_result_count_ = 0;
    uint64_t time_gvm_update_views_ns_ = 0;
    uint64_t time_gvm_update_nlf_ns_ = 0;
    uint64_t time_lvm_create_view_ns_ = 0;
    uint64_t time_lvm_destroy_view_ns_ = 0;
    uint64_t time_plain_search_ns_ = 0;
    uint64_t time_phase1off_backtracking_search_ns_ = 0;
    uint64_t time_phase1on_grouping_ns_ = 0;
    uint64_t time_phase1on_pending_build_ns_ = 0;
    uint64_t time_phase1on_plain_total_ns_ = 0;
    uint64_t time_phase1on_subemb_total_ns_ = 0;
    uint64_t time_phase1on_step2_total_ns_ = 0;
    uint64_t time_phase1on_build_from_candidates_ns_ = 0;
    uint64_t time_phase1on_mapping_ids_lookup_ns_ = 0;
    uint64_t phase3_insert_plain_search_le1_calls_ = 0;
    uint64_t phase3_insert_block_reuse_calls_ = 0;
    uint64_t phase3_delete_plain_search_calls_ = 0;
    uint64_t phase3_sum_block_reuse_entries_ = 0;
    uint64_t phase3_sum_block_reuse_fallback_plain_ = 0;
    uint64_t phase3_sum_block_reuse_skeleton_runs_ = 0;
    uint64_t phase3_sum_skeleton_invalid_depth_breaks_ = 0;
    uint64_t phase3_sum_cut_vertex_branches_ = 0;
    uint64_t phase3_sum_cut_skip_no_phase3_block_ = 0;
    uint64_t phase3_sum_block_inner_dfs_runs_ = 0;
    uint64_t phase3_sum_block_cache_hits_ = 0;
    uint64_t phase3_sum_block_cache_misses_ = 0;
    uint64_t phase3_sum_block_cache_stores_ = 0;
    struct Phase1SubqueryHitLogEntry {
        uint64_t stream_update_index = 0;
        char op = 0;
        uint32_t data_u = 0;
        uint32_t data_v = 0;
        uint32_t edge_label = 0;
        uint32_t group_id = 0;
        uint32_t group_sub_id = 0;
        size_t sub_embedding_count = 0;
    };
    std::vector<Phase1SubqueryHitLogEntry> phase1_subquery_hit_log_;
    void reset_performance_counters() {
        edge_process_count_ = 0;
        relevant_update_count_ = 0;
        phase1off_mapped_automorphism_gt1_update_count_ = 0;
        search_count_ = 0;
        positive_count_ = 0;
        result_count_ = 0;
        invalid_partial_result_count_ = 0;
        partial_result_count_ = 0;
        iso_conflict_count_ = 0;
        si_empty_count_ = 0;
        lc_empty_count_ = 0;
        phase1_before_sub_valid_check_count_ = 0;
        phase1_before_task_pruned_check_count_ = 0;
        phase1_before_sub_embeddings_check_count_ = 0;
        phase1_after_sub_embeddings_check_count_ = 0;
        phase1_compare_normal_search_result_count_ = 0;
        time_gvm_update_views_ns_ = 0;
        time_gvm_update_nlf_ns_ = 0;
        time_lvm_create_view_ns_ = 0;
        time_lvm_destroy_view_ns_ = 0;
        time_plain_search_ns_ = 0;
        time_phase1off_backtracking_search_ns_ = 0;
        time_phase1on_grouping_ns_ = 0;
        time_phase1on_pending_build_ns_ = 0;
        time_phase1on_plain_total_ns_ = 0;
        time_phase1on_subemb_total_ns_ = 0;
        time_phase1on_step2_total_ns_ = 0;
        time_phase1on_build_from_candidates_ns_ = 0;
        time_phase1on_mapping_ids_lookup_ns_ = 0;
        phase1_subquery_hit_log_.clear();
        phase3_insert_plain_search_le1_calls_ = 0;
        phase3_insert_block_reuse_calls_ = 0;
        phase3_delete_plain_search_calls_ = 0;
        phase3_sum_block_reuse_entries_ = 0;
        phase3_sum_block_reuse_fallback_plain_ = 0;
        phase3_sum_block_reuse_skeleton_runs_ = 0;
        phase3_sum_skeleton_invalid_depth_breaks_ = 0;
        phase3_sum_cut_vertex_branches_ = 0;
        phase3_sum_cut_skip_no_phase3_block_ = 0;
        phase3_sum_block_inner_dfs_runs_ = 0;
        phase3_sum_block_cache_hits_ = 0;
        phase3_sum_block_cache_misses_ = 0;
        phase3_sum_block_cache_stores_ = 0;
    }
private:

    struct QueryEdgeKey {
        uint32_t u = 0;   
        uint32_t v = 0;   
        uint32_t label = 0;
        bool operator==(const QueryEdgeKey& o) const {
            return u == o.u && v == o.v && label == o.label;
        }
    };

    struct DirectedEdgeKey {
        uint32_t src = 0;
        uint32_t dst = 0;
        uint32_t label = 0;
        bool operator==(const DirectedEdgeKey& o) const {
            return src == o.src && dst == o.dst && label == o.label;
        }
    };

    struct QueryEdgeKeyHash {
        size_t operator()(const QueryEdgeKey& k) const noexcept {
            size_t seed = 0;                    
            spp::hash_combine(seed, k.u);      
            spp::hash_combine(seed, k.v);      
            spp::hash_combine(seed, k.label);   
            return seed;                        
        }
    };

    struct DirectedEdgeKeyHash {
        size_t operator()(const DirectedEdgeKey& k) const noexcept {
            size_t seed = 0;
            spp::hash_combine(seed, k.src);
            spp::hash_combine(seed, k.dst);
            spp::hash_combine(seed, k.label);
            return seed;
        }
    };
  
    struct Phase1EdgeReuseTag {
        bool shareable = false;                
        uint32_t tag_set_id = 0;             
      
    };

    struct Phase1TagAtom {
        uint32_t sub_id = 0;
        uint32_t x = 0;
        uint32_t y = 0;
        bool operator==(const Phase1TagAtom& o) const {
            return sub_id == o.sub_id && x == o.x && y == o.y;
        }
    };

    struct Phase1TagAtomLess {
        bool operator()(const Phase1TagAtom& a, const Phase1TagAtom& b) const {
            if (a.sub_id != b.sub_id) return a.sub_id < b.sub_id;
            if (a.x != b.x) return a.x < b.x;
            return a.y < b.y;
        }
    };

    std::unordered_map<DirectedEdgeKey, Phase1EdgeReuseTag, DirectedEdgeKeyHash> phase1_edge_reuse_tag_;
    std::unordered_map<DirectedEdgeKey, std::vector<Phase1TagAtom>, DirectedEdgeKeyHash> phase1_edge_tag_atoms_;
    std::unordered_map<DirectedEdgeKey, uint32_t, DirectedEdgeKeyHash> phase1_edge_group_id_;
    std::vector<uint32_t> phase1_group_chosen_sub_id_;
    std::unordered_map<DirectedEdgeKey, uint32_t, DirectedEdgeKeyHash> phase2_edge_group_id_;
    std::vector<uint32_t> phase2_group_chosen_sub_id_;
    OrderManager om_;
    GlobalViewManager gvm_;
    LocalViewManager lvm_;
    SearchEngine sm_;

    QueryBCTreeDecomposition query_bctree_;

    struct SubqueryContext {
        std::unique_ptr<Graph> sub_query;
        std::vector<int32_t> old2new;
        std::vector<int32_t> old_vertex_rep_new;
        std::vector<int32_t> new2old;
        std::unordered_map<uint32_t, std::vector<uint32_t>> vertex_permutations;
        std::unordered_map<uint32_t, std::vector<uint32_t>> mapping_ids_index;
        std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>> mapped_new_set_by_key_and_old;
        std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>> mapped_new_orbit_peers_by_key;
        std::vector<uint32_t> sub_old_vertices;
        std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>> forced_pairs_template;
        OrderManager om;
        LocalViewManager lvm;
        SearchEngine sm;
        std::vector<std::vector<std::pair<uint32_t, uint32_t>>> sub_vertex_nlf;
        bool initialized = false;
    };


    struct isoSubqueryContext {
        std::unique_ptr<Graph> sub_query;
        std::vector<int32_t> old2new;
        std::vector<std::vector<uint32_t>> mappings;
        OrderManager om;
        LocalViewManager lvm;
        SearchEngine sm;
        bool initialized = false;
    };

    isoSubqueryContext* get_isosubquery_context(uint32_t group_id);

    std::unordered_map<uint32_t, SubqueryContext> subquery_cache_;
    std::unordered_map<uint32_t, isoSubqueryContext> isosubquery_cache_;

    const Graph* data_graph_ = nullptr;
    std::string query_graph_file_;
    std::string data_graph_file_;
    /**
     * 结合 data_graph 与 query_graph_file，批量生成所有 sub_id 对应的 SubqueryContext 并缓存到 subquery_cache_：
     * - 针对每个 order/subgraph 或自同构子目录，生成相关的子查询图、索引、局部/全局视图等辅助对象
     * - 用于后续每次流式更新能直接命中 cache，避免重复构建和 IO 读写消耗
     * @param data_graph 数据图指针
     * @param query_graph_file 查询图文件名或所属目录
     */

    void build_subquery_cache_in_memory(const Graph* query_graph, const Graph* data_graph, uint64_t max_autos = 200000);
    void build_isosubquery_cache_in_memory(const Graph* query_graph, const Graph* data_graph, uint64_t max_isos = 200000);
    void preprocess_phase1_grouping_tags(const Graph* query_graph);
    void preprocess_phase2_grouping_tags(const Graph* query_graph);
    /**
     * 根据 sub_id 查询（已提前构建好的）SubqueryContext
     * 用例：流处理/更新时，直接命中对应子图及其索引、视图等全部辅助结构
     * @param sub_id 目标子结构编号
     * @return 指向对应 SubqueryContext 的指针，若不存在时返回 nullptr
     */
    SubqueryContext* get_subquery_context(uint32_t sub_id);

    void preprocess_phase1_off(const Graph* query_graph, const Graph* data_graph);
    void preprocess_phase2_on(const Graph* query_graph, const Graph* data_graph);
    void preprocess_phase4_on(const Graph* query_graph, const Graph* data_graph);
    void preprocess_phase5_on(const Graph* query_graph, const Graph* data_graph);
    void preprocess_phase1_on(const Graph* query_graph, const Graph* data_graph);
    void preprocess_phase3_on(const Graph* query_graph, const Graph* data_graph);
    uint64_t execute_phase1_off(const Graph* query_graph, const Update& update,
                                bool enable_local_view, bool enable_search);
    uint64_t execute_phase2_on(const Graph* query_graph, const Update& update,
                               bool enable_local_view, bool enable_search);
    uint64_t execute_phase4_on(const Graph* query_graph, const Update& update,
                               bool enable_local_view, bool enable_search);
    uint64_t execute_phase1_on(const Graph* query_graph, const Update& update,
                               bool enable_local_view, bool enable_search);
    uint64_t execute_phase3_on(const Graph* query_graph, const Update& update,
                               bool enable_local_view, bool enable_search);
    uint64_t execute_phase5_on(const Graph* query_graph, const Update& update,
                               bool enable_local_view, bool enable_search);

    bool phase4_enabled_ = false;
    bool phase2_enabled_ = false;
    bool phase1_enabled_ = false;
    bool phase3_enabled_ = false;
    bool phase5_enabled_ = false;

public:
    void set_query_graph_file(std::string query_graph_file) { query_graph_file_ = std::move(query_graph_file); }
    void set_data_graph_file(std::string data_graph_file) { data_graph_file_ = std::move(data_graph_file); }
    void set_phase4_enabled(bool enabled) { phase4_enabled_ = enabled; }
    bool phase4_enabled() const { return phase4_enabled_; }
    void set_phase2_enabled(bool enabled) { phase2_enabled_ = enabled; }
    bool phase2_enabled() const { return phase2_enabled_; }
    void set_phase1_enabled(bool enabled) { phase1_enabled_ = enabled; }
    bool phase1_enabled() const { return phase1_enabled_; }
    void set_phase3_enabled(bool enabled) { phase3_enabled_ = enabled; }
    bool phase3_enabled() const { return phase3_enabled_; }
    void set_phase5_enabled(bool enabled) { phase5_enabled_ = enabled; }
    bool phase5_enabled() const { return phase5_enabled_; }
    
    void initialize(const Graph *query_graph, const Graph *data_graph, uint64_t target_embedding_num);
    void preprocess(const Graph* query_graph, const Graph* data_graph);
    uint64_t execute(const Graph *query_graph, const Update &update, bool enable_local_view, bool enable_search);
    void release();
    void evaluate_view_update(const Graph *query_graph, const Graph *data_graph,
                              const std::vector<Update> &stream);
    void evaluate_search(const Graph *query_graph, const Graph *data_graph, const std::vector<Update> &stream);
    void print_metrics();
};

#endif 
