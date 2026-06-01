#ifndef RAPIDMATCH_LOCAL_VIEW_MANAGER_H
#define RAPIDMATCH_LOCAL_VIEW_MANAGER_H
#include "streaming_type.h"
#include "global_view_manager.h"
#include "relation/local_edge_view.h"
struct LocalViewWorkspace {
    std::vector<uint32_t> flag_array_;
    std::vector<uint32_t> updated_flag_;
    std::vector<uint32_t> si_buffer_;
    std::vector<uint32_t*> candidate_set_pointer_;
    std::vector<uint32_t> candidate_set_size_;
    void initialize(const Graph *query_graph, const Graph *data_graph) {
        flag_array_.resize(data_graph->getVerticesCount(), 0);
        updated_flag_.resize(data_graph->getVerticesCount());
        si_buffer_.resize(data_graph->getVerticesCount());
        candidate_set_size_.resize(query_graph->getVerticesCount());
        candidate_set_pointer_.resize(query_graph->getVerticesCount());
    }
    void release() {
        flag_array_.clear();
        updated_flag_.clear();
        si_buffer_.clear();
        candidate_set_size_.clear();
        candidate_set_pointer_.clear();
    }
};
struct LocalViewStorage {
    std::vector<std::vector<uint32_t>> candidates_store_;
    std::vector<LocalEdgeView> views_;
    std::vector<uint32_t> buffer_pool_;
    spp::sparse_hash_map<Edge, uint32_t> edge_view_mapping_;
    void initialize(const Graph *query_graph) {
        candidates_store_.resize(query_graph->getVerticesCount());
        for (auto& candidates_set : candidates_store_) {
            candidates_set.reserve(1024);
        }
        views_.resize(query_graph->getEdgesCount() * 2);
        buffer_pool_.reserve(1024 * 1024);
        edge_view_mapping_.reserve(256);
    }

    void initialize_candidates_only(const Graph *query_graph) {
        candidates_store_.resize(query_graph->getVerticesCount());
    }
    void clear_candidates_only() {
        for (auto& candidate_set : candidates_store_) {
            candidate_set.clear();
        }
    }
    void clear_views_only() {
        for (auto& view : views_) {
            view.clear();
        }
        buffer_pool_.clear();
        edge_view_mapping_.clear();
    }
    void clear_view_contents() {
        for (auto& view : views_) {
            view.clear();
        }
        for (auto& candidate_set : candidates_store_) {
            candidate_set.clear();
        }
        buffer_pool_.clear();
        edge_view_mapping_.clear();
    }

    void release() {
        candidates_store_.clear();
        views_.clear();
        buffer_pool_.clear();
        edge_view_mapping_.clear();
    }
};

class LocalViewSnapshot {
    friend class LocalViewManager;

public:
    uint64_t build_visited_neighbor_count_ = 0;
    uint64_t generate_visited_neighbor_count_ = 0;
    uint64_t first_vertex_neighbor_ = 0;

private:
    LocalViewStorage storage_;

public:
    LocalViewSnapshot() = default;
    explicit LocalViewSnapshot(LocalViewStorage&& storage) : storage_(std::move(storage)) {}

    LocalEdgeView* get_view(uint32_t view_id);
    LocalEdgeView* get_view(Edge query_edge);
    uint32_t get_view_id(Edge query_edge);
    uint32_t* get_candidate_set(uint32_t u);
    uint32_t get_candidate_set_size(uint32_t u);
};
class LocalCandidatesSnapshot {
    friend class LocalViewManager;
public:
    uint64_t generate_visited_neighbor_count_ = 0;
    uint64_t first_vertex_neighbor_ = 0;
    size_t get_query_vertex_count() const { return candidates_store_.size(); }
    const std::vector<uint32_t>& get_candidate_set(uint32_t u) const { return candidates_store_[u]; }
private:
    std::vector<std::vector<uint32_t>> candidates_store_;
};

class LocalViewManager {
public:
    uint64_t build_visited_neighbor_count_;
    uint64_t generate_visited_neighbor_count_;
    uint64_t first_vertex_neighbor_;
private:
    LocalViewWorkspace* workspace_ = nullptr;
    std::unique_ptr<LocalViewWorkspace> owned_workspace_;
    LocalViewStorage storage_;
    uint32_t updated_count_;
private:
    bool optimized_generate_local_candidates(const Graph *query_graph, OrdersPerEdge &orders,
                                             GlobalViewManager &gvm, Edge exclude_data_edge);

    bool optimized_generate_local_candidates_v2(const Graph *query_graph, OrdersPerEdge &orders,
                                                GlobalViewManager &gvm, Edge exclude_data_edge);

    bool optimized_generate_local_candidates_v2_for_subgraph(const Graph *query_graph, const Graph *full_query_graph,
                                                             OrdersPerEdge &orders,
                                                             GlobalViewManager &gvm, Edge exclude_data_edge,
                                                             const std::vector<int32_t>* new2old_sub_to_full,
                                                             const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>*
                                                                 sub_vertex_nlf,
                                                             const std::vector<std::vector<uint32_t>>*
                                                                 mapped_new_orbit_peers);
    bool optimized_generate_local_candidates_v2_with_subgraph(
            const Graph *query_graph,
            OrdersPerEdge &orders,
            GlobalViewManager &gvm,
            Edge exclude_data_edge,
            const LocalCandidatesSnapshot &sub_snapshot,
            const std::vector<std::vector<uint32_t>> &mapped_new_set_by_old,
            const Graph *sub_query,
            const std::vector<int32_t> *new2old_sub_to_full);

    void optimized_build_local_view(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm);

    void optimized_build_local_view_for_subgraph(const Graph *query_graph, OrdersPerEdge &orders,
                                                 GlobalViewManager &gvm,
                                                 const std::vector<int32_t>* new2old_sub_to_full);

    bool prune_local_candidates(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn);

    void set_adjacent_update_candidates_flag(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                             uint32_t encoded_v0, uint32_t encoded_v1);

    void set_update_candidates_flag(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                    uint32_t encoded_v0, uint32_t encoded_v1);

    uint32_t select_bn_with_minimum_degree_sum(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn);

    uint32_t select_bn_with_minimum_degree_sum_for_subgraph(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                                            const std::vector<int32_t>* new2old_sub_to_full);
    void set_adjacent_update_candidates_flag_for_subgraph(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                                          uint32_t encoded_v0, uint32_t encoded_v1,
                                                          const std::vector<int32_t>* new2old_sub_to_full);
    void set_update_candidates_flag_for_subgraph(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                                 uint32_t encoded_v0, uint32_t encoded_v1,
                                                 const std::vector<int32_t>* new2old_sub_to_full);
    bool prune_local_candidates_for_subgraph(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                             const std::vector<int32_t>* new2old_sub_to_full);

public:
    LocalViewManager() {}
    ~LocalViewManager() { release(); }

    void initialize(const Graph *query_graph, const Graph *data_graph);
    void initialize(const Graph *query_graph, const Graph *data_graph, LocalViewWorkspace* shared_workspace);
    void release();

    bool create_view(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm, Edge data_edge);

    bool create_candidates_only(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm, Edge data_edge);
    bool create_candidates_only_for_subgraph(const Graph *query_graph, const Graph *full_query_graph,
                                             OrdersPerEdge &orders, GlobalViewManager &gvm,
                                             Edge data_edge, const std::vector<int32_t>* new2old_sub_to_full,
                                             const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>* sub_vertex_nlf =
                                                 nullptr,
                                             const std::vector<std::vector<uint32_t>>* mapped_new_orbit_peers =
                                                 nullptr);
    bool create_candidates_only_for_subgraph_iso(const Graph *query_graph, const Graph *full_query_graph,
                                                  OrdersPerEdge &orders, GlobalViewManager &gvm,
                                                  Edge data_edge, const std::vector<int32_t>* old2new_full_to_sub,
                                                  uint32_t sub_vertex_count);
    bool create_candidates_only_with_sub_constraint(
            const Graph *query_graph,
            OrdersPerEdge &orders,
            GlobalViewManager &gvm,
            Edge data_edge,
            const LocalCandidatesSnapshot &sub_snapshot,
            const std::vector<std::vector<uint32_t>> &mapped_new_set_by_old,
            const Graph *sub_query,
            const std::vector<int32_t> *new2old_sub_to_full);

    LocalCandidatesSnapshot extract_candidates_snapshot(const Graph *query_graph);

    void restore_candidates_from(LocalCandidatesSnapshot&& snap);
  
    void build_view_from_candidates(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm);
 
    void build_view_from_candidates_for_subgraph(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm,
                                                 const std::vector<int32_t>* new2old_sub_to_full);

 
    bool create_view_for_subgraph(const Graph *query_graph, const Graph *full_query_graph,
                                  OrdersPerEdge &orders, GlobalViewManager &gvm, Edge data_edge,
                                  const std::vector<int32_t>* new2old_sub_to_full,
                                  const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>* sub_vertex_nlf =
                                      nullptr,
                                  const std::vector<std::vector<uint32_t>>* mapped_new_orbit_peers = nullptr);

  
    void destroy_view();
    LocalViewSnapshot extract_snapshot(const Graph *query_graph);
    void restore_from(LocalViewSnapshot&& snap);
    LocalEdgeView* get_view(uint32_t view_id);
    LocalEdgeView* get_view(Edge query_edge);
    uint32_t get_view_id(Edge query_edge);
    uint32_t* get_candidate_set(uint32_t u);
    uint32_t get_candidate_set_size(uint32_t u);
};


#endif 