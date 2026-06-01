#include <chrono>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <algorithm>
#include "local_view_manager.h"
namespace {
inline uint32_t map_sub_query_vertex_for_full_gvm(const std::vector<int32_t>* new2old, uint32_t sub_u) {
    if (!new2old) return sub_u;
    assert(sub_u < new2old->size());
    const int32_t fu = (*new2old)[sub_u];
    assert(fu >= 0);
    return static_cast<uint32_t>(fu);
}
} 
#include "streaming_config.h"
#include "computesetintersection.h"
void LocalViewManager::initialize(const Graph *query_graph, const Graph *data_graph) {
    owned_workspace_ = std::make_unique<LocalViewWorkspace>();
    workspace_ = owned_workspace_.get();
    workspace_->initialize(query_graph, data_graph);
    storage_.initialize(query_graph);
}

void LocalViewManager::initialize(const Graph *query_graph, const Graph *data_graph, LocalViewWorkspace* shared_workspace) {
    owned_workspace_.reset();
    workspace_ = shared_workspace;
    workspace_->initialize(query_graph, data_graph);
    storage_.initialize(query_graph);
}

void LocalViewManager::release() {
    storage_.release();
    if (owned_workspace_) {
        owned_workspace_->release();
        owned_workspace_.reset();
    }
    workspace_ = nullptr;
}

void LocalViewManager::destroy_view() {
    storage_.clear_view_contents();
}

LocalEdgeView *LocalViewManager::get_view(uint32_t view_id) {
    assert(view_id < storage_.views_.size());
    return &storage_.views_[view_id];
}

LocalEdgeView *LocalViewManager::get_view(Edge query_edge) {
    auto it = storage_.edge_view_mapping_.find(query_edge);
    if (it != storage_.edge_view_mapping_.end()) {
        return &storage_.views_[it->second];
    }
    return nullptr;
}

uint32_t LocalViewManager::get_view_id(Edge query_edge) {
    auto it = storage_.edge_view_mapping_.find(query_edge);
    if (it != storage_.edge_view_mapping_.end()) {
        return it->second;
    }
    // Local view does not exist.
    return storage_.edge_view_mapping_.size();
}

uint32_t *LocalViewManager::get_candidate_set(uint32_t u) {
    return storage_.candidates_store_[u].data();
}

uint32_t LocalViewManager::get_candidate_set_size(uint32_t u) {
    return storage_.candidates_store_[u].size();
}

LocalViewSnapshot LocalViewManager::extract_snapshot(const Graph *query_graph) {
    LocalViewSnapshot snapshot(std::move(storage_));
    snapshot.build_visited_neighbor_count_ = build_visited_neighbor_count_;
    snapshot.generate_visited_neighbor_count_ = generate_visited_neighbor_count_;
    snapshot.first_vertex_neighbor_ = first_vertex_neighbor_;

    storage_ = LocalViewStorage();
    if (query_graph != nullptr) {
        storage_.initialize(query_graph);
    }
    return snapshot;
}

void LocalViewManager::restore_from(LocalViewSnapshot&& snap) {
    destroy_view();
    storage_ = std::move(snap.storage_);
    generate_visited_neighbor_count_ = snap.generate_visited_neighbor_count_;
    build_visited_neighbor_count_ = snap.build_visited_neighbor_count_;
    first_vertex_neighbor_ = snap.first_vertex_neighbor_;
}

LocalEdgeView *LocalViewSnapshot::get_view(uint32_t view_id) {
    assert(view_id < storage_.views_.size());
    return &storage_.views_[view_id];
}

LocalEdgeView *LocalViewSnapshot::get_view(Edge query_edge) {
    auto it = storage_.edge_view_mapping_.find(query_edge);
    if (it != storage_.edge_view_mapping_.end()) {
        return &storage_.views_[it->second];
    }
    return nullptr;
}

uint32_t LocalViewSnapshot::get_view_id(Edge query_edge) {
    auto it = storage_.edge_view_mapping_.find(query_edge);
    if (it != storage_.edge_view_mapping_.end()) {
        return it->second;
    }
    return storage_.edge_view_mapping_.size();
}

uint32_t *LocalViewSnapshot::get_candidate_set(uint32_t u) {
    return storage_.candidates_store_[u].data();
}

uint32_t LocalViewSnapshot::get_candidate_set_size(uint32_t u) {
    return storage_.candidates_store_[u].size();
}

bool
LocalViewManager::create_view(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm, Edge data_edge) {
    build_visited_neighbor_count_  = 0;
    generate_visited_neighbor_count_ = 0;
    first_vertex_neighbor_ = 0;
    if (storage_.views_.empty()) {
        storage_.initialize(query_graph);
    } else {
        storage_.clear_view_contents();
    }
    bool is_valid = optimized_generate_local_candidates_v2(query_graph, orders, gvm, data_edge);
    if (!is_valid)
        return false;
     optimized_build_local_view(query_graph, orders, gvm);
    return true;
}

bool LocalViewManager::create_candidates_only(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm,
                                              Edge data_edge) {
    build_visited_neighbor_count_  = 0;
    generate_visited_neighbor_count_ = 0;
    first_vertex_neighbor_ = 0;
    if (storage_.candidates_store_.empty()) {
        storage_.initialize_candidates_only(query_graph);
    } else {
        storage_.clear_candidates_only();
    }
    return optimized_generate_local_candidates_v2(query_graph, orders, gvm, data_edge);
}

bool LocalViewManager::create_candidates_only_for_subgraph(const Graph *query_graph, const Graph *full_query_graph,
                                                           OrdersPerEdge &orders,
                                                           GlobalViewManager &gvm, Edge data_edge,
                                                           const std::vector<int32_t>* new2old_sub_to_full,
                                                           const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>*
                                                               sub_vertex_nlf,
                                                           const std::vector<std::vector<uint32_t>>* mapped_new_orbit_peers) {
    build_visited_neighbor_count_  = 0;
    generate_visited_neighbor_count_ = 0;
    first_vertex_neighbor_ = 0;
    if (storage_.candidates_store_.empty()) {
        storage_.initialize_candidates_only(query_graph);
    } else {
        storage_.clear_candidates_only();
    }
    return optimized_generate_local_candidates_v2_for_subgraph(query_graph, full_query_graph, orders, gvm, data_edge,
                                                                new2old_sub_to_full, sub_vertex_nlf,
                                                                mapped_new_orbit_peers);
}


bool LocalViewManager::create_candidates_only_for_subgraph_iso(const Graph *query_graph, const Graph *full_query_graph,
                                                                OrdersPerEdge &orders,
                                                                GlobalViewManager &gvm, Edge data_edge,
                                                                const std::vector<int32_t>* old2new_full_to_sub,
                                                                uint32_t sub_vertex_count) {
    std::vector<int32_t> new2old(sub_vertex_count, -1);
    if (old2new_full_to_sub) {
        for (uint32_t old_u = 0; old_u < old2new_full_to_sub->size(); ++old_u) {
            int32_t sub_u = (*old2new_full_to_sub)[old_u];
            if (sub_u >= 0 && static_cast<uint32_t>(sub_u) < sub_vertex_count) {
                new2old[sub_u] = static_cast<int32_t>(old_u);
            }
        }
    }
    return create_candidates_only_for_subgraph(query_graph, full_query_graph, orders, gvm, data_edge,
                                                &new2old, nullptr, nullptr);
}


bool LocalViewManager::create_candidates_only_with_sub_constraint(
        const Graph *query_graph,
        OrdersPerEdge &orders,
        GlobalViewManager &gvm,
        Edge data_edge,
        const LocalCandidatesSnapshot &sub_snapshot,
        const std::vector<std::vector<uint32_t>> &mapped_new_set_by_old,
        const Graph *sub_query,
        const std::vector<int32_t> *new2old_sub_to_full) {
    build_visited_neighbor_count_  = 0;
    generate_visited_neighbor_count_ = 0;
    first_vertex_neighbor_ = 0;
    if (storage_.candidates_store_.empty()) {
        storage_.initialize_candidates_only(query_graph);
    } else {
        storage_.clear_candidates_only();
    }
    return optimized_generate_local_candidates_v2_with_subgraph(
            query_graph,
            orders,
            gvm,
            data_edge,
            sub_snapshot,
            mapped_new_set_by_old,
            sub_query,
            new2old_sub_to_full);
}

LocalCandidatesSnapshot LocalViewManager::extract_candidates_snapshot(const Graph *query_graph) {
    LocalCandidatesSnapshot snapshot;
    snapshot.generate_visited_neighbor_count_ = generate_visited_neighbor_count_;
    snapshot.first_vertex_neighbor_ = first_vertex_neighbor_;
    const uint32_t n = query_graph ? query_graph->getVerticesCount() : 0;
    snapshot.candidates_store_.resize(n);
    for (uint32_t u = 0; u < n; ++u) {
        snapshot.candidates_store_[u] = std::move(storage_.candidates_store_[u]);
    }
    storage_.candidates_store_.assign(n, {});
    return snapshot;
}

void LocalViewManager::restore_candidates_from(LocalCandidatesSnapshot&& snap) {
    storage_.candidates_store_ = std::move(snap.candidates_store_);
    generate_visited_neighbor_count_ = snap.generate_visited_neighbor_count_;
    first_vertex_neighbor_ = snap.first_vertex_neighbor_;
}

void LocalViewManager::build_view_from_candidates(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm) {
    if (storage_.views_.empty()) {
        storage_.views_.resize(query_graph->getEdgesCount() * 2);
        storage_.buffer_pool_.reserve(1024 * 1024);
        storage_.edge_view_mapping_.reserve(256);
    } else {
        storage_.clear_views_only();
    }

    optimized_build_local_view(query_graph, orders, gvm);
}

void LocalViewManager::build_view_from_candidates_for_subgraph(const Graph *query_graph, OrdersPerEdge &orders,
                                                               GlobalViewManager &gvm,
                                                               const std::vector<int32_t>* new2old_sub_to_full) {
    if (storage_.views_.empty()) {
        storage_.views_.resize(query_graph->getEdgesCount() * 2);
        storage_.buffer_pool_.reserve(1024 * 1024);
        storage_.edge_view_mapping_.reserve(256);
    } else {
        storage_.clear_views_only();
    }
    optimized_build_local_view_for_subgraph(query_graph, orders, gvm, new2old_sub_to_full);
}

bool
LocalViewManager::create_view_for_subgraph(const Graph *query_graph, const Graph *full_query_graph,
                                           OrdersPerEdge &orders, GlobalViewManager &gvm,
                                           Edge data_edge, const std::vector<int32_t>* new2old_sub_to_full,
                                           const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>* sub_vertex_nlf,
                                           const std::vector<std::vector<uint32_t>>* mapped_new_orbit_peers) {
    build_visited_neighbor_count_  = 0;
    generate_visited_neighbor_count_ = 0;
    first_vertex_neighbor_ = 0;
    if (storage_.views_.empty()) {
        storage_.initialize(query_graph);
    } else {
        storage_.clear_view_contents();
    }
    bool is_valid = optimized_generate_local_candidates_v2_for_subgraph(query_graph, full_query_graph, orders, gvm, data_edge,
                                                                        new2old_sub_to_full, sub_vertex_nlf,
                                                                        mapped_new_orbit_peers);
    if (!is_valid)
        return false;
    optimized_build_local_view_for_subgraph(query_graph, orders, gvm, new2old_sub_to_full);
    return true;
}


bool LocalViewManager::optimized_generate_local_candidates(const Graph *query_graph, OrdersPerEdge &orders,
                                                           GlobalViewManager &gvm, Edge exclude_data_edge) {
    auto& indexing_order = orders.indexing_order_;
    auto& indexing_order_bn = orders.indexing_order_bn_;
    auto& indexing_order_bn_offset = orders.indexing_order_bn_offset_;
    uint32_t u0 = indexing_order[0];
    uint32_t u1 = indexing_order[1];
    uint32_t v0 = exclude_data_edge.first;
    uint32_t v1 = exclude_data_edge.second;
    if (!gvm.nlf_check(u0, v0)
        || !gvm.nlf_check(u1, v1)) {
        return false;
    }
    if (query_graph->getVertexDegree(u0) > 1) {
        uint32_t global_encoded_v0 = gvm.get_encoded_id(u0, v0);
        storage_.candidates_store_[u0].push_back(global_encoded_v0);
    }
    else {
        storage_.candidates_store_[u0].push_back(v0);
    }

    if (query_graph->getVertexDegree(u1) > 1) {
        uint32_t global_encoded_v1 = gvm.get_encoded_id(u1, v1);
        storage_.candidates_store_[u1].push_back(global_encoded_v1);
    }
    else {
        storage_.candidates_store_[u1].push_back(v1);
    }

    for (uint32_t i = 2; i < indexing_order.size(); ++i) {
        uint32_t u = indexing_order[i];

        // Skip leaf node.
        if (query_graph->getVertexDegree(u) == 1)
            continue;

        uint32_t begin = indexing_order_bn_offset[i];
        uint32_t end = indexing_order_bn_offset[i + 1];
        updated_count_ = 0;

        uint32_t flag_value = 0;
        for (uint32_t j = begin; j < end; ++j) {
            uint32_t uu = indexing_order_bn[j];
            auto gv = gvm.get_nlf_view({uu, u});
            auto reverse_gv = gvm.get_nlf_view({u, uu});

            bool flag = false;
            if (end - begin >= 2) {
                flag = (indexing_order_bn[begin] == u0 && indexing_order_bn[begin + 1] == u1);
            }
            for (auto vv : storage_.candidates_store_[uu]) {
                uint32_t vv_nbrs_count;
                auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);
                generate_visited_neighbor_count_ += vv_nbrs_count;
                if ((uu == indexing_order[0] || uu == indexing_order[1]) && flag)
                    first_vertex_neighbor_ += vv_nbrs_count;

                // If it is the first bn or the cost of binary search is lower than that of scan, then use the first approach.
                if ((j == begin) || vv_nbrs_count < 1024 || vv_nbrs_count < updated_count_ * 32) {
                    for (uint32_t k = 0; k < vv_nbrs_count; ++k) {
                        uint32_t v = vv_nbrs[k];

                        if (workspace_->flag_array_[v] == flag_value) {
                            workspace_->flag_array_[v] += 1;
                            if (flag_value == 0) {
                                workspace_->updated_flag_[updated_count_++] = v;
                            }
                        }
                    }
                }
                else {
                    for (uint32_t k = 0; k < updated_count_; ++k) {
                        uint32_t v = workspace_->updated_flag_[k];
                        if (workspace_->flag_array_[v] == flag_value) {
                            uint32_t v_nbrs_count;
                            auto v_nbrs = reverse_gv->get_neighbor(v, v_nbrs_count);
                            if (vv_nbrs_count < v_nbrs_count) {
                                auto it = std::lower_bound(vv_nbrs, vv_nbrs + vv_nbrs_count, v);
                                if (it != vv_nbrs + vv_nbrs_count && *it == v) {
                                    workspace_->flag_array_[v] += 1;
                                }
                            }
                            else {
                                auto it = std::lower_bound(v_nbrs, v_nbrs + v_nbrs_count, vv);
                                if (it != v_nbrs + v_nbrs_count && *it == vv) {
                                    workspace_->flag_array_[v] += 1;
                                }
                            }
                        }
                    }
                }
            }

            flag_value += 1;
        }
        uint32_t local_encoded_v0 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(u, v0)) {
            local_encoded_v0 = gvm.get_encoded_id(u, v0);
        }

        uint32_t local_encoded_v1 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(u, v1)) {
            local_encoded_v1 = gvm.get_encoded_id(u, v1);
        }

        for (uint32_t j = 0; j < updated_count_; ++j) {
            uint32_t v = workspace_->updated_flag_[j];
            if (workspace_->flag_array_[v] == flag_value && v != local_encoded_v0 && v != local_encoded_v1) {
                storage_.candidates_store_[u].push_back(v);
            }

            workspace_->flag_array_[v] = 0;
        }

        if (storage_.candidates_store_[u].empty())
            return false;
    }

    for (auto& candidate_set : storage_.candidates_store_) {
        std::sort(candidate_set.begin(), candidate_set.end());
    }

    return true;
}

void
LocalViewManager::optimized_build_local_view(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm) {
    auto& matching_order = orders.matching_order_;
    auto& matching_order_bn = orders.matching_order_bn_;
    auto& matching_order_bn_offset = orders.matching_order_bn_offset_;
    auto& matching_order_view_mappings = orders.matching_order_view_mappings_;
    auto& matching_order_relation_edge_type = orders.matching_order_edge_type_;
    uint32_t local_view_id = 0;
    for (uint32_t i = 0; i < matching_order.size(); ++i) {
        uint32_t u = matching_order[i];
        updated_count_ = 0;
        for (uint32_t j = matching_order_bn_offset[i]; j < matching_order_bn_offset[i + 1]; ++j) {
            uint32_t uu = matching_order_bn[j];
            RelationEdgeType type = matching_order_relation_edge_type[j];
            Edge query_edge = {uu, u};
            Edge reverse_query_edge = {u, uu};
            if (query_graph->getVertexDegree(u) == 1) {
                matching_order_view_mappings[j]= gvm.get_view_id(query_edge);
                continue;
            }
            matching_order_view_mappings[j] = local_view_id;
            storage_.edge_view_mapping_[query_edge] = local_view_id;
            auto &lv = storage_.views_[local_view_id++];
            if (type == RelationEdgeType::REGULAR) {
                if (updated_count_ == 0) {
                    for (uint32_t k = 0; k < storage_.candidates_store_[u].size(); ++k) {
                        workspace_->flag_array_[storage_.candidates_store_[u][k]] = k + 1;
                        workspace_->updated_flag_[updated_count_++] = storage_.candidates_store_[u][k];
                    }
                }
                auto gv = gvm.get_nlf_view(query_edge);
                auto reverse_gv = gvm.get_nlf_view(reverse_query_edge);
                for (auto vv : storage_.candidates_store_[uu]) {
                    uint32_t begin_pos = storage_.buffer_pool_.size();
                    uint32_t vv_nbrs_count;
                    auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);
                    build_visited_neighbor_count_ += vv_nbrs_count;
                    if (vv_nbrs_count < 1024 || vv_nbrs_count < updated_count_ * 32) {
                        for (uint32_t k = 0; k < vv_nbrs_count; ++k) {
                            uint32_t v = vv_nbrs[k];
                            if (workspace_->flag_array_[v] > 0) {
                                storage_.buffer_pool_.push_back(workspace_->flag_array_[v] - 1);
                            }
                        }
                    }
               
                    else {
                        for (auto v : storage_.candidates_store_[u]) {
                            uint32_t v_nbrs_count;
                            auto v_nbrs = reverse_gv->get_neighbor(v, v_nbrs_count);
                            if (vv_nbrs_count < v_nbrs_count) {
                                auto it = std::lower_bound(vv_nbrs, vv_nbrs + vv_nbrs_count, v);
                                if (it != vv_nbrs + vv_nbrs_count && *it == v) {
                                    storage_.buffer_pool_.push_back(workspace_->flag_array_[v] - 1);
                                }
                            }
                            else {
                                auto it = std::lower_bound(v_nbrs, v_nbrs + v_nbrs_count, vv);
                                if (it != v_nbrs + v_nbrs_count && *it == vv) {
                                    storage_.buffer_pool_.push_back(workspace_->flag_array_[v] - 1);
                                }
                            }
                        }
                    }
                    lv.trie_.emplace_back(begin_pos, storage_.buffer_pool_.size());
                    lv.cardinality_ += storage_.buffer_pool_.size() - begin_pos;
                }
            }
        }
        for (uint32_t j = 0; j < updated_count_; ++j) {
            uint32_t v = workspace_->updated_flag_[j];
            workspace_->flag_array_[v] = 0;
        }
    }
    for (uint32_t i = 0; i < local_view_id; ++i) {
        storage_.views_[i].data_ = storage_.buffer_pool_.data();
    }
    for (uint32_t u = 0; u < query_graph->getVerticesCount(); ++u) {
        if (query_graph->getVertexDegree(u) > 1) {
            for (uint32_t i = 0; i < storage_.candidates_store_[u].size(); ++i) {
                storage_.candidates_store_[u][i] = gvm.get_id(u, storage_.candidates_store_[u][i]);
            }
        }
    }
}
void
LocalViewManager::optimized_build_local_view_for_subgraph(const Graph *query_graph, OrdersPerEdge &orders,
                                                          GlobalViewManager &gvm,
                                                          const std::vector<int32_t>* new2old_sub_to_full) {
    auto& matching_order = orders.matching_order_;
    auto& matching_order_bn = orders.matching_order_bn_;
    auto& matching_order_bn_offset = orders.matching_order_bn_offset_;
    auto& matching_order_view_mappings = orders.matching_order_view_mappings_;
    auto& matching_order_relation_edge_type = orders.matching_order_edge_type_;
    uint32_t local_view_id = 0;
    for (uint32_t i = 0; i < matching_order.size(); ++i) {
        uint32_t u = matching_order[i]; 
        updated_count_ = 0;
        for (uint32_t j = matching_order_bn_offset[i]; j < matching_order_bn_offset[i + 1]; ++j) {
            uint32_t uu = matching_order_bn[j]; 
            RelationEdgeType type = matching_order_relation_edge_type[j];
            Edge query_edge = {uu, u}; 
            const uint32_t full_uu = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, uu);
            const uint32_t full_u = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, u);
            const Edge full_query_edge = {full_uu, full_u};
            const Edge full_reverse_edge = {full_u, full_uu};

            if (query_graph->getVertexDegree(u) == 1) {
                matching_order_view_mappings[j]= gvm.get_view_id(full_query_edge);
                continue;
            }

            matching_order_view_mappings[j] = local_view_id;
            storage_.edge_view_mapping_[query_edge] = local_view_id;

            auto &lv = storage_.views_[local_view_id++];

            if (type == RelationEdgeType::REGULAR) {
                if (updated_count_ == 0) {
                    for (uint32_t k = 0; k < storage_.candidates_store_[u].size(); ++k) {
                        workspace_->flag_array_[storage_.candidates_store_[u][k]] = k + 1;
                        workspace_->updated_flag_[updated_count_++] = storage_.candidates_store_[u][k];
                    }
                }
                auto gv = gvm.get_nlf_view(full_query_edge);
                auto reverse_gv = gvm.get_nlf_view(full_reverse_edge);
                for (auto vv : storage_.candidates_store_[uu]) {
                    uint32_t begin_pos = storage_.buffer_pool_.size();
                    uint32_t vv_nbrs_count;
                    auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);
                    build_visited_neighbor_count_ += vv_nbrs_count;
                    if (vv_nbrs_count < 1024 || vv_nbrs_count < updated_count_ * 32) {
                        for (uint32_t k = 0; k < vv_nbrs_count; ++k) {
                            uint32_t v = vv_nbrs[k];
                            if (workspace_->flag_array_[v] > 0) {
                                storage_.buffer_pool_.push_back(workspace_->flag_array_[v] - 1);
                            }
                        }
                    }
                    else {
                        for (auto v : storage_.candidates_store_[u]) {
                            uint32_t v_nbrs_count;
                            auto v_nbrs = reverse_gv->get_neighbor(v, v_nbrs_count);
                            if (vv_nbrs_count < v_nbrs_count) {
                                auto it = std::lower_bound(vv_nbrs, vv_nbrs + vv_nbrs_count, v);
                                if (it != vv_nbrs + vv_nbrs_count && *it == v) {
                                    storage_.buffer_pool_.push_back(workspace_->flag_array_[v] - 1);
                                }
                            }
                            else {
                                auto it = std::lower_bound(v_nbrs, v_nbrs + v_nbrs_count, vv);
                                if (it != v_nbrs + v_nbrs_count && *it == vv) {
                                    storage_.buffer_pool_.push_back(workspace_->flag_array_[v] - 1);
                                }
                            }
                        }
                    }

                    lv.trie_.emplace_back(begin_pos, storage_.buffer_pool_.size());
                    lv.cardinality_ += storage_.buffer_pool_.size() - begin_pos;
                }
            }
        }

        for (uint32_t j = 0; j < updated_count_; ++j) {
            uint32_t v = workspace_->updated_flag_[j];
            workspace_->flag_array_[v] = 0;
        }
    }

    for (uint32_t i = 0; i < local_view_id; ++i) {
        storage_.views_[i].data_ = storage_.buffer_pool_.data();
    }

    for (uint32_t u = 0; u < query_graph->getVerticesCount(); ++u) {
        if (query_graph->getVertexDegree(u) > 1) {
            const uint32_t full_u = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, u);
            for (uint32_t i = 0; i < storage_.candidates_store_[u].size(); ++i) {
                storage_.candidates_store_[u][i] = gvm.get_id(full_u, storage_.candidates_store_[u][i]);
            }
        }
    }
}

bool LocalViewManager::optimized_generate_local_candidates_v2(const Graph *query_graph, OrdersPerEdge &orders,
                                                              GlobalViewManager &gvm, Edge exclude_data_edge) {
    auto& indexing_order = orders.indexing_order_;
    auto& indexing_order_bn = orders.indexing_order_bn_;
    auto& indexing_order_bn_offset = orders.indexing_order_bn_offset_;
    uint32_t u0 = indexing_order[0];
    uint32_t u1 = indexing_order[1];
    uint32_t v0 = exclude_data_edge.first;
    uint32_t v1 = exclude_data_edge.second;
    const uint32_t full_u0 = u0;
    const uint32_t full_u1 = u1;

    if (!gvm.nlf_check(full_u0, v0)
        || !gvm.nlf_check(full_u1, v1)) {
        return false;
    }
    if (query_graph->getVertexDegree(u0) > 1) {
        uint32_t global_encoded_v0 = gvm.get_encoded_id(full_u0, v0);
        storage_.candidates_store_[u0].push_back(global_encoded_v0);
    }
    else {
        storage_.candidates_store_[u0].push_back(v0);
    }
    if (query_graph->getVertexDegree(u1) > 1) {
        uint32_t global_encoded_v1 = gvm.get_encoded_id(full_u1, v1);
        storage_.candidates_store_[u1].push_back(global_encoded_v1);
    }
    else {
        storage_.candidates_store_[u1].push_back(v1);
    }
    workspace_->candidate_set_size_[u0] = 1;
    workspace_->candidate_set_pointer_[u0] = storage_.candidates_store_[u0].data();
    workspace_->candidate_set_size_[u1] = 1;
    workspace_->candidate_set_pointer_[u1] = storage_.candidates_store_[u1].data();
    for (uint32_t i = 2; i < orders.triangle_end_index_; ++i) {
        uint32_t u = indexing_order[i];
        const uint32_t full_u = u;
        auto gv0 = gvm.get_nlf_view({full_u0, full_u});
        uint32_t v0_nbr_count;
        auto v0_nbr = gv0->get_neighbor(workspace_->candidate_set_pointer_[u0][0], v0_nbr_count);
        auto gv1 = gvm.get_nlf_view({full_u1, full_u});
        uint32_t v1_nbr_count;
        auto v1_nbr = gv1->get_neighbor(workspace_->candidate_set_pointer_[u1][0], v1_nbr_count);
        uint32_t lc = 0;
        ComputeSetIntersection::ComputeCandidates(v0_nbr, v0_nbr_count, v1_nbr, v1_nbr_count, workspace_->updated_flag_.data(), lc);
        if (lc == 0)
            return false;
        storage_.candidates_store_[u].insert(storage_.candidates_store_[u].end(), workspace_->updated_flag_.begin(), workspace_->updated_flag_.begin() + lc);
        workspace_->candidate_set_pointer_[u] = storage_.candidates_store_[u].data();
        workspace_->candidate_set_size_[u] = storage_.candidates_store_[u].size();
    }

    // ----- 阶段二：与「更新边」相邻的顶点段 [triangle_end_index_, adjacent_end_index_) -----
    // 从已匹配的 backward uu 拉 NLF 邻居作为 u 的候选；小集合时拷入 candidates_store_ 并去掉更新边的两端在 u 上的编码（避免重复占用）。
    for (uint32_t i = orders.triangle_end_index_; i < orders.adjacent_end_index_; ++i) {
        uint32_t u = indexing_order[i];

        if (query_graph->getVertexDegree(u) == 1)
            continue;

        uint32_t uu = indexing_order_bn[indexing_order_bn_offset[i]];
        const uint32_t full_u = u;
        const uint32_t full_uu = uu;
        auto gv = gvm.get_nlf_view({full_uu, full_u});
        workspace_->candidate_set_pointer_[u] = gv->get_neighbor(workspace_->candidate_set_pointer_[uu][0], workspace_->candidate_set_size_[u]);

        if (workspace_->candidate_set_size_[u] < 1024) {
            uint32_t local_encoded_v0 = std::numeric_limits<uint32_t>::max();
            if (gvm.candidate_check(full_u, v0)) {
                local_encoded_v0 = gvm.get_encoded_id(full_u, v0);
            }

            uint32_t local_encoded_v1 = std::numeric_limits<uint32_t>::max();
            if (gvm.candidate_check(full_u, v1)) {
                local_encoded_v1 = gvm.get_encoded_id(full_u, v1);
            }

            for (uint32_t j = 0; j < workspace_->candidate_set_size_[u]; ++j) {
                uint32_t v = workspace_->candidate_set_pointer_[u][j];
                if (v != local_encoded_v0 && v != local_encoded_v1) {
                    storage_.candidates_store_[u].push_back(v);
                }
            }

            workspace_->candidate_set_pointer_[u] = storage_.candidates_store_[u].data();
            workspace_->candidate_set_size_[u] = storage_.candidates_store_[u].size();
        }

        if (workspace_->candidate_set_size_[u] == 0)
            return false;
    }

    // ----- 阶段三：仅为「阶段二涉及」的顶点重排处理顺序 -----
    // 将 indexing_order[2..adjacent_end) 按当前 candidate_set_size 升序排序，使后续剪枝先处理更紧的顶点；
    // 并重建与之一致的 backward 列表 optimized_indexing_order_bn / offset（只覆盖 i < adjacent_end_index_ 的大循环）。
        std::vector<uint32_t> optimized_indexing_order(indexing_order.begin() + 2,
                                                       indexing_order.begin() + orders.adjacent_end_index_);
        std::sort(optimized_indexing_order.begin(), optimized_indexing_order.end(), [this](uint32_t a, uint32_t b) {
            return workspace_->candidate_set_size_[a] < workspace_->candidate_set_size_[b];
        });

        std::vector<bool> visited(query_graph->getVerticesCount(), false);
        std::vector<uint32_t> optimized_indexing_order_bn;
        std::vector<uint32_t> optimized_indexing_order_bn_offset;

        optimized_indexing_order_bn_offset.push_back(0);
        for (auto u : optimized_indexing_order) {
            uint32_t u_nbrs_count;
            auto u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
            for (uint32_t i = 0; i < u_nbrs_count; ++i) {
                uint32_t uu = u_nbrs[i];
                if (visited[uu]) {
                    optimized_indexing_order_bn.push_back(uu);
                }
            }
            optimized_indexing_order_bn_offset.push_back(optimized_indexing_order_bn.size());
            visited[u] = true;
        }


    // ----- 阶段四：对 indexing_order[2..] 逐顶点用 backward 邻居做 NLF 传播 + 剪枝 -----
    for (uint32_t i = 2; i < indexing_order.size(); ++i) {
        uint32_t u;
        std::vector<uint32_t> bn;
        if (i < orders.adjacent_end_index_) {
            // 与更新边「相邻」段：u 与 backward 列表来自阶段三的重排结果。
            uint32_t index = i - 2;
            u = optimized_indexing_order[index];
            uint32_t begin = optimized_indexing_order_bn_offset[index];
            uint32_t end = optimized_indexing_order_bn_offset[index + 1];
            bn.insert(bn.end(), optimized_indexing_order_bn.begin() + begin, optimized_indexing_order_bn.begin() + end);
        }
        else {
            // 其余深度：仍用原始 orders 的 indexing_order 与 indexing_order_bn。
            u = indexing_order[i];
            uint32_t begin = indexing_order_bn_offset[i];
            uint32_t end = indexing_order_bn_offset[i + 1];
            bn.insert(bn.end(), indexing_order_bn.begin() + begin, indexing_order_bn.begin() + end);
        }
        if (query_graph->getVertexDegree(u) == 1 || bn.empty())
            continue;

        const uint32_t full_u = u;

        // 更新边两端在「当前查询顶点 u」的编码；用于在 flag 传播中排除已锚定的数据顶点。
        uint32_t local_encoded_v0 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(full_u, v0)) {
            local_encoded_v0 = gvm.get_encoded_id(full_u, v0);
        }

        uint32_t local_encoded_v1 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(full_u, v1)) {
            local_encoded_v1 = gvm.get_encoded_id(full_u, v1);
        }

        if (i < orders.adjacent_end_index_) {
            // 相邻段：在 workspace_ 上打标，表示「各 backward 约束下 u 仍可能的数据顶点」。
            set_adjacent_update_candidates_flag(gvm, u, bn, local_encoded_v0, local_encoded_v1);
        }
        else {
            // 更深顶点：另一套 flag 传播（与是否紧贴本次更新边有关）。
            set_update_candidates_flag(gvm, u, bn, local_encoded_v0, local_encoded_v1);
        }

        // 剪枝：根据 NLF 与 backward 一致性收紧 C(u)；若不可能则整条局部视图失败。
        if (!prune_local_candidates(gvm, u, bn))
            return false;
    }

    // ----- 阶段五：把仍只在指针里的候选落盘到 candidates_store_，并排序 -----
    // 高度数顶点：若尚未 push 进 storage（大候选时前面可能只维护了 pointer），这里补齐；最后排序便于 build 阶段与二分。
    for (auto u : indexing_order) {
        if (query_graph->getVertexDegree(u) > 1) {
            if (storage_.candidates_store_[u].empty()) {
                for (uint32_t i = 0; i < workspace_->candidate_set_size_[u]; ++i) {
                    storage_.candidates_store_[u].push_back(workspace_->candidate_set_pointer_[u][i]);
                }
            }
            std::sort(storage_.candidates_store_[u].begin(), storage_.candidates_store_[u].end());
            workspace_->candidate_set_pointer_[u] = storage_.candidates_store_[u].data();
            workspace_->candidate_set_size_[u] = storage_.candidates_store_[u].size();
        }
    }
    return true;
}

// optimized_generate_local_candidates_v2_for_subgraph：同 v2；凡 gvm 的查询顶点/边键经 new2old 映到全图（nullptr 恒等）。
bool LocalViewManager::optimized_generate_local_candidates_v2_for_subgraph(const Graph *query_graph, const Graph *full_query_graph,
                                                                          OrdersPerEdge &orders,
                                                                          GlobalViewManager &gvm, Edge exclude_data_edge,
                                                                          const std::vector<int32_t>* new2old_sub_to_full,
                                                                          const std::vector<std::vector<std::pair<uint32_t, uint32_t>>>*
                                                                              sub_vertex_nlf,
                                                                          const std::vector<std::vector<uint32_t>>*
                                                                              mapped_new_orbit_peers) {
    auto& indexing_order = orders.indexing_order_;
    auto& indexing_order_bn = orders.indexing_order_bn_;
    auto& indexing_order_bn_offset = orders.indexing_order_bn_offset_;
    uint32_t u0 = indexing_order[0];
    uint32_t u1 = indexing_order[1];
    uint32_t v0 = exclude_data_edge.first;
    uint32_t v1 = exclude_data_edge.second;
    const uint32_t fu0 = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, u0);
    const uint32_t fu1 = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, u1);
    if (full_query_graph->getVertexDegree(fu0) > 1) {
        uint32_t global_encoded_v0 = gvm.get_encoded_id(fu0, v0);
        storage_.candidates_store_[u0].push_back(global_encoded_v0);
    }
    else {
        storage_.candidates_store_[u0].push_back(v0);
    }

    if (full_query_graph->getVertexDegree(fu1) > 1) {
        uint32_t global_encoded_v1 = gvm.get_encoded_id(fu1, v1);
        storage_.candidates_store_[u1].push_back(global_encoded_v1);
    }
    else {
        storage_.candidates_store_[u1].push_back(v1);
    }
    workspace_->candidate_set_size_[u0] = 1;
    workspace_->candidate_set_pointer_[u0] = storage_.candidates_store_[u0].data();
    workspace_->candidate_set_size_[u1] = 1;
    workspace_->candidate_set_pointer_[u1] = storage_.candidates_store_[u1].data();
    for (uint32_t i = 2; i < orders.triangle_end_index_; ++i) {
        uint32_t u = indexing_order[i];
        const uint32_t full_u = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, u);
        const uint32_t full_u0 = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, u0);
        const uint32_t full_u1 = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, u1);

        auto gv0 = gvm.get_nlf_view({full_u0, full_u});
        if (!gv0) {
            return false;
        }
        uint32_t v0_nbr_count;
        auto v0_nbr = gv0->get_neighbor(workspace_->candidate_set_pointer_[u0][0], v0_nbr_count);
        auto gv1 = gvm.get_nlf_view({full_u1, full_u});
        if (!gv1) {
            return false;
        }
        uint32_t v1_nbr_count;
        auto v1_nbr = gv1->get_neighbor(workspace_->candidate_set_pointer_[u1][0], v1_nbr_count);
        uint32_t lc = 0;
        ComputeSetIntersection::ComputeCandidates(v0_nbr, v0_nbr_count, v1_nbr, v1_nbr_count, workspace_->updated_flag_.data(), lc);
        if (lc == 0) {

            return false;
        }
        storage_.candidates_store_[u].insert(storage_.candidates_store_[u].end(), workspace_->updated_flag_.begin(), workspace_->updated_flag_.begin() + lc);
        workspace_->candidate_set_pointer_[u] = storage_.candidates_store_[u].data();
        workspace_->candidate_set_size_[u] = storage_.candidates_store_[u].size();
    }

    for (uint32_t i = orders.triangle_end_index_; i < orders.adjacent_end_index_; ++i) {
        uint32_t u = indexing_order[i];
        const uint32_t full_u = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, u);
        if (full_query_graph->getVertexDegree(full_u) == 1)
            continue;
        uint32_t uu = indexing_order_bn[indexing_order_bn_offset[i]];
        const uint32_t full_uu = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, uu);
        auto gv = gvm.get_nlf_view({full_uu, full_u});
        if (!gv) {
            return false;
        }
        workspace_->candidate_set_pointer_[u] = gv->get_neighbor(workspace_->candidate_set_pointer_[uu][0], workspace_->candidate_set_size_[u]);
        if (workspace_->candidate_set_size_[u] < 1024) {
            uint32_t local_encoded_v0 = std::numeric_limits<uint32_t>::max();
            if (gvm.candidate_check(full_u, v0)) {
                local_encoded_v0 = gvm.get_encoded_id(full_u, v0);
            }
            uint32_t local_encoded_v1 = std::numeric_limits<uint32_t>::max();
            if (gvm.candidate_check(full_u, v1)) {
                local_encoded_v1 = gvm.get_encoded_id(full_u, v1);
            }
            for (uint32_t j = 0; j < workspace_->candidate_set_size_[u]; ++j) {
                uint32_t v = workspace_->candidate_set_pointer_[u][j];
                if (v != local_encoded_v0 && v != local_encoded_v1) {
                    storage_.candidates_store_[u].push_back(v);
                }
            }
            workspace_->candidate_set_pointer_[u] = storage_.candidates_store_[u].data();
            workspace_->candidate_set_size_[u] = storage_.candidates_store_[u].size();
        }
        if (workspace_->candidate_set_size_[u] == 0) {
                return false;
        }
    }

        std::vector<uint32_t> optimized_indexing_order(indexing_order.begin() + 2,
                                                       indexing_order.begin() + orders.adjacent_end_index_);
        std::sort(optimized_indexing_order.begin(), optimized_indexing_order.end(), [this](uint32_t a, uint32_t b) {
            return workspace_->candidate_set_size_[a] < workspace_->candidate_set_size_[b];
        });

        std::vector<bool> visited(query_graph->getVerticesCount(), false);
        std::vector<uint32_t> optimized_indexing_order_bn;
        std::vector<uint32_t> optimized_indexing_order_bn_offset;
        optimized_indexing_order_bn_offset.push_back(0);
        for (auto u : optimized_indexing_order) {
            uint32_t u_nbrs_count;
            auto u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
            for (uint32_t i = 0; i < u_nbrs_count; ++i) {
                uint32_t uu = u_nbrs[i];
                if (visited[uu]) {
                    optimized_indexing_order_bn.push_back(uu);
                }
            }
            optimized_indexing_order_bn_offset.push_back(optimized_indexing_order_bn.size());
            visited[u] = true;
        }

    for (uint32_t i = 2; i < indexing_order.size(); ++i) {
        uint32_t u;
        std::vector<uint32_t> bn;
        if (i < orders.adjacent_end_index_) {
            uint32_t index = i - 2;
            u = optimized_indexing_order[index];
            uint32_t begin = optimized_indexing_order_bn_offset[index];
            uint32_t end = optimized_indexing_order_bn_offset[index + 1];
            bn.insert(bn.end(), optimized_indexing_order_bn.begin() + begin, optimized_indexing_order_bn.begin() + end);
        }
        else {
            u = indexing_order[i];
            uint32_t begin = indexing_order_bn_offset[i];
            uint32_t end = indexing_order_bn_offset[i + 1];
            bn.insert(bn.end(), indexing_order_bn.begin() + begin, indexing_order_bn.begin() + end);
        }
        const uint32_t full_u = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, u);
        if (full_query_graph->getVertexDegree(full_u) == 1 || bn.empty())
            continue;

        uint32_t local_encoded_v0 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(full_u, v0)) {
            local_encoded_v0 = gvm.get_encoded_id(full_u, v0);
        }
        uint32_t local_encoded_v1 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(full_u, v1)) {
            local_encoded_v1 = gvm.get_encoded_id(full_u, v1);
        }

        if (i < orders.adjacent_end_index_) {
            set_adjacent_update_candidates_flag_for_subgraph(gvm, u, bn, local_encoded_v0, local_encoded_v1,
                                                             new2old_sub_to_full);
        }
        else {
            set_update_candidates_flag_for_subgraph(gvm, u, bn, local_encoded_v0, local_encoded_v1, new2old_sub_to_full);
        }

        if (!prune_local_candidates_for_subgraph(gvm, u, bn, new2old_sub_to_full)) {
            return false;
        }
    }
    for (auto u : indexing_order) {
        const uint32_t full_u = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, u);
        if (full_query_graph->getVertexDegree(full_u) > 1
    ) {
            if (storage_.candidates_store_[u].empty()) {
                for (uint32_t i = 0; i < workspace_->candidate_set_size_[u]; ++i) {
                    storage_.candidates_store_[u].push_back(workspace_->candidate_set_pointer_[u][i]);
                }
            }
            std::sort(storage_.candidates_store_[u].begin(), storage_.candidates_store_[u].end());
            workspace_->candidate_set_pointer_[u] = storage_.candidates_store_[u].data();
            workspace_->candidate_set_size_[u] = storage_.candidates_store_[u].size();
        }
    }
    return true;
}

bool LocalViewManager::optimized_generate_local_candidates_v2_with_subgraph(
        const Graph *query_graph,
        OrdersPerEdge &orders,
        GlobalViewManager &gvm,
        Edge exclude_data_edge,
        const LocalCandidatesSnapshot &sub_snapshot,
        const std::vector<std::vector<uint32_t>> &mapped_new_set_by_old,
        const Graph *sub_query,
        const std::vector<int32_t> *new2old_sub_to_full) {
    auto& indexing_order = orders.indexing_order_;
    auto& indexing_order_bn = orders.indexing_order_bn_;
    auto& indexing_order_bn_offset = orders.indexing_order_bn_offset_;

    if (!query_graph || indexing_order.size() < 2) {
        return false;
    }

    const uint32_t u0 = indexing_order[0];
    const uint32_t u1 = indexing_order[1];
    const uint32_t v0 = exclude_data_edge.first;
    const uint32_t v1 = exclude_data_edge.second;
    const uint32_t full_u0 = u0;
    const uint32_t full_u1 = u1;
    const uint32_t qn = query_graph->getVerticesCount();
    const uint32_t old_limit = static_cast<uint32_t>(mapped_new_set_by_old.size());
    const uint32_t sub_k = static_cast<uint32_t>(sub_snapshot.candidates_store_.size());
    std::vector<uint8_t> seeded(qn, 0);
    for (uint32_t old_u = 0; old_u < old_limit; ++old_u) {
        const auto& mapped_new_list = mapped_new_set_by_old[old_u];
   
        auto& dst = storage_.candidates_store_[old_u];
        dst.clear();
        const bool full_is_non_leaf = (query_graph->getVertexDegree(old_u) > 1);
        const uint32_t encoded_cap = gvm.get_encoded_candidate_set_size(old_u);
        for (uint32_t mapped_new_u : mapped_new_list) {
            if (mapped_new_u >= sub_k) {
                continue;
            }
            const auto& src = sub_snapshot.candidates_store_[mapped_new_u];
            for (uint32_t cand : src) {
                uint32_t raw_v = 0;
                bool have_raw = false;

                if (mapped_new_u >= new2old_sub_to_full->size()) {
                    continue;
                }
                const int32_t full_src_i = (*new2old_sub_to_full)[mapped_new_u];
                if (full_src_i < 0) {
                    continue;
                }
                const uint32_t full_src = static_cast<uint32_t>(full_src_i);
                if (query_graph->getVertexDegree(full_src) > 1) {
                    const uint32_t src_enc_cap = gvm.get_encoded_candidate_set_size(full_src);
                    if (cand >= src_enc_cap) {
                        continue;
                    }
                    raw_v = gvm.get_id(full_src, cand);
                } else {
                    raw_v = cand;
                }
                have_raw = true;

                if (!have_raw) {
                    continue;
                }
                if (!gvm.candidate_check(old_u, raw_v)) {
                    continue;
                }
                if (full_is_non_leaf) {
                    const uint32_t enc = gvm.get_encoded_id(old_u, raw_v);
                    if (enc >= encoded_cap) {
                        continue;
                    }
                    dst.push_back(enc);
                } else {
                    dst.push_back(raw_v);
                }
            }
        }

        std::sort(dst.begin(), dst.end());
        dst.erase(std::unique(dst.begin(), dst.end()), dst.end());
        workspace_->candidate_set_pointer_[old_u] = dst.data();
        workspace_->candidate_set_size_[old_u] = static_cast<uint32_t>(dst.size());
        seeded[old_u] = !dst.empty() ? 1 : 0;
    }

    
    auto bind_anchor_singleton = [&](uint32_t anchor_u, uint32_t data_v) -> bool {
        uint32_t anchor_candidate = data_v;
        if (query_graph->getVertexDegree(anchor_u) > 1) {
            if (!gvm.candidate_check(anchor_u, data_v)) {
                return false;
            }
            anchor_candidate = gvm.get_encoded_id(anchor_u, data_v);
        }

        auto& dst = storage_.candidates_store_[anchor_u];
        if (!dst.empty() && !std::binary_search(dst.begin(), dst.end(), anchor_candidate)) {
            return false;
        }
        dst.clear();
        dst.push_back(anchor_candidate);
        workspace_->candidate_set_pointer_[anchor_u] = dst.data();
        workspace_->candidate_set_size_[anchor_u] = 1;
        seeded[anchor_u] = 1;
        return true;
    };

    if (!bind_anchor_singleton(u0, v0) || !bind_anchor_singleton(u1, v1)) {
        return false;
    }

    for (uint32_t old_u = 0; old_u < old_limit; ++old_u) {
        if (old_u == u0 || old_u == u1 || !seeded[old_u])
            continue;
        auto& dd = storage_.candidates_store_[old_u];
        if (dd.empty()) continue;
        const bool old_is_non_leaf = (query_graph->getVertexDegree(old_u) > 1);
        uint32_t rm0 = v0, rm1 = v1;
        if (old_is_non_leaf) {
            if (gvm.candidate_check(old_u, v0))
                rm0 = gvm.get_encoded_id(old_u, v0);
            else
                rm0 = std::numeric_limits<uint32_t>::max();
            if (gvm.candidate_check(old_u, v1))
                rm1 = gvm.get_encoded_id(old_u, v1);
            else
                rm1 = std::numeric_limits<uint32_t>::max();
        }
        auto it_end = std::remove_if(dd.begin(), dd.end(), [rm0, rm1](uint32_t x) {
            return x == rm0 || x == rm1;
        });
        dd.erase(it_end, dd.end());
        workspace_->candidate_set_pointer_[old_u] = dd.data();
        workspace_->candidate_set_size_[old_u] = static_cast<uint32_t>(dd.size());
    }

    for (uint32_t i = 2; i < orders.triangle_end_index_; ++i) {
        uint32_t u = indexing_order[i];
        if (u >= seeded.size() || seeded[u]) {
            continue;
        }
        const uint32_t full_u = u;

        auto gv0 = gvm.get_nlf_view({full_u0, full_u});
        uint32_t v0_nbr_count;
        auto v0_nbr = gv0->get_neighbor(workspace_->candidate_set_pointer_[u0][0], v0_nbr_count);

        auto gv1 = gvm.get_nlf_view({full_u1, full_u});
        uint32_t v1_nbr_count;
        auto v1_nbr = gv1->get_neighbor(workspace_->candidate_set_pointer_[u1][0], v1_nbr_count);

        uint32_t lc = 0;
        ComputeSetIntersection::ComputeCandidates(v0_nbr, v0_nbr_count, v1_nbr, v1_nbr_count,
                                                  workspace_->updated_flag_.data(), lc);

        if (lc == 0) {
            return false;
        }

        auto& dst = storage_.candidates_store_[u];
        dst.clear();
        dst.insert(dst.end(), workspace_->updated_flag_.begin(), workspace_->updated_flag_.begin() + lc);
        workspace_->candidate_set_pointer_[u] = dst.data();
        workspace_->candidate_set_size_[u] = static_cast<uint32_t>(dst.size());
    }

    for (uint32_t i = orders.triangle_end_index_; i < orders.adjacent_end_index_; ++i) {
        uint32_t u = indexing_order[i];
        if (u >= seeded.size() || seeded[u]) {
            continue;
        }

        if (query_graph->getVertexDegree(u) == 1) {
            continue;
        }

        uint32_t uu = indexing_order_bn[indexing_order_bn_offset[i]];
        const uint32_t full_u = u;
        const uint32_t full_uu = uu;
        auto gv = gvm.get_nlf_view({full_uu, full_u});
        workspace_->candidate_set_pointer_[u] =
                gv->get_neighbor(workspace_->candidate_set_pointer_[uu][0], workspace_->candidate_set_size_[u]);

        if (workspace_->candidate_set_size_[u] < 1024) {
            uint32_t local_encoded_v0 = std::numeric_limits<uint32_t>::max();
            if (gvm.candidate_check(full_u, v0)) {
                local_encoded_v0 = gvm.get_encoded_id(full_u, v0);
            }

            uint32_t local_encoded_v1 = std::numeric_limits<uint32_t>::max();
            if (gvm.candidate_check(full_u, v1)) {
                local_encoded_v1 = gvm.get_encoded_id(full_u, v1);
            }

            auto& dst = storage_.candidates_store_[u];
            dst.clear();
            for (uint32_t j = 0; j < workspace_->candidate_set_size_[u]; ++j) {
                uint32_t v = workspace_->candidate_set_pointer_[u][j];
                if (v != local_encoded_v0 && v != local_encoded_v1) {
                    dst.push_back(v);
                }
            }

            workspace_->candidate_set_pointer_[u] = dst.data();
            workspace_->candidate_set_size_[u] = static_cast<uint32_t>(dst.size());
        }

        if (workspace_->candidate_set_size_[u] == 0) {
            return false;
        }
    }

    std::vector<uint32_t> optimized_indexing_order(indexing_order.begin() + 2,
                                                   indexing_order.begin() + orders.adjacent_end_index_);
    std::sort(optimized_indexing_order.begin(), optimized_indexing_order.end(), [this](uint32_t a, uint32_t b) {
        return workspace_->candidate_set_size_[a] < workspace_->candidate_set_size_[b];
    });

    std::vector<bool> visited(query_graph->getVerticesCount(), false);
    std::vector<uint32_t> optimized_indexing_order_bn;
    std::vector<uint32_t> optimized_indexing_order_bn_offset;
    optimized_indexing_order_bn_offset.push_back(0);
    for (auto u : optimized_indexing_order) {
        uint32_t u_nbrs_count;
        auto u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        for (uint32_t i = 0; i < u_nbrs_count; ++i) {
            uint32_t uu = u_nbrs[i];
            if (visited[uu]) {
                optimized_indexing_order_bn.push_back(uu);
            }
        }
        optimized_indexing_order_bn_offset.push_back(optimized_indexing_order_bn.size());
        visited[u] = true;
    }

    for (uint32_t i = 2; i < indexing_order.size(); ++i) {
        uint32_t u;
        std::vector<uint32_t> bn;
        if (i < orders.adjacent_end_index_) {
            uint32_t index = i - 2;
            u = optimized_indexing_order[index];
            uint32_t begin = optimized_indexing_order_bn_offset[index];
            uint32_t end = optimized_indexing_order_bn_offset[index + 1];
            bn.insert(bn.end(), optimized_indexing_order_bn.begin() + begin, optimized_indexing_order_bn.begin() + end);
        } else {
            u = indexing_order[i];
            uint32_t begin = indexing_order_bn_offset[i];
            uint32_t end = indexing_order_bn_offset[i + 1];
            bn.insert(bn.end(), indexing_order_bn.begin() + begin, indexing_order_bn.begin() + end);
        }

        if (u < seeded.size() && seeded[u]) {
            continue;
        }
        if (query_graph->getVertexDegree(u) == 1 || bn.empty()) {
            continue;
        }
        const uint32_t full_u = u;
        uint32_t local_encoded_v0 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(full_u, v0)) {
            local_encoded_v0 = gvm.get_encoded_id(full_u, v0);
        }
        uint32_t local_encoded_v1 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(full_u, v1)) {
            local_encoded_v1 = gvm.get_encoded_id(full_u, v1);
        }
        if (i < orders.adjacent_end_index_) {
            set_adjacent_update_candidates_flag(gvm, u, bn, local_encoded_v0, local_encoded_v1);
        } else {
            set_update_candidates_flag(gvm, u, bn, local_encoded_v0, local_encoded_v1);
        }
        if (!prune_local_candidates(gvm, u, bn)) {
            return false;
        }
    }
    for (auto u : indexing_order) {
        if (query_graph->getVertexDegree(u) > 1) {
            if (storage_.candidates_store_[u].empty()) {
                for (uint32_t i = 0; i < workspace_->candidate_set_size_[u]; ++i) {
                    storage_.candidates_store_[u].push_back(workspace_->candidate_set_pointer_[u][i]);
                }
            }
            std::sort(storage_.candidates_store_[u].begin(), storage_.candidates_store_[u].end());
            workspace_->candidate_set_pointer_[u] = storage_.candidates_store_[u].data();
            workspace_->candidate_set_size_[u] = static_cast<uint32_t>(storage_.candidates_store_[u].size());
        }
    }
    return true;
}

bool LocalViewManager::prune_local_candidates(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn) {
    uint32_t flag_value = 1;
    uint32_t current_valid_count = updated_count_;

    for (auto uu : bn) {
        auto gv = gvm.get_nlf_view({uu, u});
        auto reverse_gv = gvm.get_nlf_view({u, uu});

        uint32_t target_valid_count = current_valid_count;
        current_valid_count = 0;

        for (uint32_t j = 0; j < workspace_->candidate_set_size_[uu]; ++j) {
            uint32_t vv = workspace_->candidate_set_pointer_[uu][j];
            uint32_t vv_nbrs_count;
            auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);

            // If the cost of binary search is lower than that of scan, then use the first approach.
            if (vv_nbrs_count < 1024 || vv_nbrs_count < updated_count_ * 32) {
                for (uint32_t k = 0; k < vv_nbrs_count; ++k) {
                    uint32_t v = vv_nbrs[k];

                    if (workspace_->flag_array_[v] == flag_value) {
                        workspace_->flag_array_[v] += 1;
                        current_valid_count += 1;
                    }
                }
            }
            else {
                for (uint32_t k = 0; k < updated_count_; ++k) {
                    uint32_t v = workspace_->updated_flag_[k];
                    if (workspace_->flag_array_[v] == flag_value) {
                        uint32_t v_nbrs_count;
                        auto v_nbrs = reverse_gv->get_neighbor(v, v_nbrs_count);
                        if (vv_nbrs_count < v_nbrs_count) {
                            auto it = std::lower_bound(vv_nbrs, vv_nbrs + vv_nbrs_count, v);
                            if (it != vv_nbrs + vv_nbrs_count && *it == v) {
                                workspace_->flag_array_[v] += 1;
                                current_valid_count += 1;
                            }
                        }
                        else {
                            auto it = std::lower_bound(v_nbrs, v_nbrs + v_nbrs_count, vv);
                            if (it != v_nbrs + v_nbrs_count && *it == vv) {
                                workspace_->flag_array_[v] += 1;
                                current_valid_count += 1;
                            }
                        }
                    }
                }
            }

            if (current_valid_count == target_valid_count)
                break;
        }

        flag_value += 1;
    }

    bool push = storage_.candidates_store_[u].empty();
    uint32_t local_pos = 0;
    for (uint32_t j = 0; j < updated_count_; ++j) {
        uint32_t v = workspace_->updated_flag_[j];
        if (workspace_->flag_array_[v] == flag_value) {
            if (push)
                storage_.candidates_store_[u].push_back(v);
            else
                storage_.candidates_store_[u][local_pos++] = v;
        }

        workspace_->flag_array_[v] = 0;
    }

    if (!push)
        storage_.candidates_store_[u].resize(local_pos);

    workspace_->candidate_set_pointer_[u] = storage_.candidates_store_[u].data();
    workspace_->candidate_set_size_[u] = storage_.candidates_store_[u].size();

    if (workspace_->candidate_set_size_[u] == 0)
        return false;

    return true;
}

void
LocalViewManager::set_adjacent_update_candidates_flag(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                                      uint32_t encoded_v0, uint32_t encoded_v1) {
    bool is_flag_based = false;
    if (workspace_->candidate_set_size_[u] < 1024) {
        is_flag_based = true;
    }
    else {
        for (auto uu : bn) {
            if (workspace_->candidate_set_size_[u] < 32 * workspace_->candidate_set_size_[uu]) {
                is_flag_based = true;
                break;
            }
        }
    }

    if (is_flag_based) {
        updated_count_ = 0;
        for (uint32_t i = 0; i < workspace_->candidate_set_size_[u]; ++i) {
            uint32_t v = workspace_->candidate_set_pointer_[u][i];
            if (v != encoded_v0 && v != encoded_v1) {
                workspace_->flag_array_[v] = 1;
                workspace_->updated_flag_[updated_count_++] = v;
            }
        }
    }
    else {
        uint32_t uu = select_bn_with_minimum_degree_sum(gvm, u, bn);
        auto gv = gvm.get_nlf_view({uu, u});

        updated_count_ = 0;
        for (uint32_t i = 0; i < workspace_->candidate_set_size_[uu]; ++i) {
            uint32_t vv = workspace_->candidate_set_pointer_[uu][i];

            uint32_t vv_nbrs_count;
            auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);

            uint32_t lc = 0;
            ComputeSetIntersection::ComputeCandidates(workspace_->candidate_set_pointer_[u], workspace_->candidate_set_size_[u],
                                                      vv_nbrs, vv_nbrs_count, workspace_->si_buffer_.data(), lc);


            for (uint32_t j = 0; j < lc; ++j) {
                uint32_t v = workspace_->si_buffer_[j];

                if (workspace_->flag_array_[v] == 0 && v != encoded_v0 && v != encoded_v1) {
                    workspace_->flag_array_[v] = 1;
                    workspace_->updated_flag_[updated_count_++] = v;
                }
            }
            if (workspace_->candidate_set_size_[u] == updated_count_)
                break;
        }
    }
}

void LocalViewManager::set_update_candidates_flag(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                                  uint32_t encoded_v0, uint32_t encoded_v1) {

    uint32_t selected_bn = select_bn_with_minimum_degree_sum(gvm, u, bn);
    auto gv = gvm.get_nlf_view({selected_bn, u});
    updated_count_ = 0;
    for (uint32_t i = 0; i < workspace_->candidate_set_size_[selected_bn]; ++i) {
        uint32_t vv = workspace_->candidate_set_pointer_[selected_bn][i];
        uint32_t vv_nbrs_count;
        auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);

        for (uint32_t j = 0; j < vv_nbrs_count; ++j) {
            uint32_t v = vv_nbrs[j];
            if (workspace_->flag_array_[v] == 0 && v != encoded_v0 && v != encoded_v1) {
                workspace_->flag_array_[v] = 1;
                workspace_->updated_flag_[updated_count_++] = v;
            }
        }
    }
}

uint32_t
LocalViewManager::select_bn_with_minimum_degree_sum(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn) {
    std::sort(bn.begin(), bn.end(), [this](const uint32_t a, const uint32_t b) -> bool {
        return workspace_->candidate_set_size_[a] < workspace_->candidate_set_size_[b];
    });

    uint64_t min_degree_sum = std::numeric_limits<uint64_t>::max();
    uint32_t selected_bn = 0;
    uint32_t selected_idx = 0;
    for (uint32_t i = 0; i < bn.size(); ++i) {
        uint32_t uu = bn[i];
        auto gv = gvm.get_nlf_view({uu, u});
        uint64_t local_degree_sum = 0;
        if (workspace_->candidate_set_size_[uu] < min_degree_sum) {
            for (uint32_t j = 0; j < workspace_->candidate_set_size_[uu]; ++j) {
                uint32_t vv = workspace_->candidate_set_pointer_[uu][j];
                local_degree_sum += gv->get_neighbor_num(vv);
            }

            if (local_degree_sum < min_degree_sum) {
                min_degree_sum = local_degree_sum;
                selected_bn = uu;
                selected_idx = i;
            }
        }
    }

    bn.erase(bn.begin() + selected_idx);
    return selected_bn;
}

uint32_t
LocalViewManager::select_bn_with_minimum_degree_sum_for_subgraph(GlobalViewManager &gvm, uint32_t u,
                                                                 std::vector<uint32_t> &bn,
                                                                 const std::vector<int32_t>* new2old_sub_to_full) {
    const uint32_t full_u = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, u);
    std::sort(bn.begin(), bn.end(), [this](const uint32_t a, const uint32_t b) -> bool {
        return workspace_->candidate_set_size_[a] < workspace_->candidate_set_size_[b];
    });
    uint64_t min_degree_sum = std::numeric_limits<uint64_t>::max();
    uint32_t selected_bn = 0;
    uint32_t selected_idx = 0;
    for (uint32_t i = 0; i < bn.size(); ++i) {
        uint32_t uu = bn[i];
        const uint32_t full_uu = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, uu);
        auto gv = gvm.get_nlf_view({full_uu, full_u});
        uint64_t local_degree_sum = 0;
        if (workspace_->candidate_set_size_[uu] < min_degree_sum) {
            for (uint32_t j = 0; j < workspace_->candidate_set_size_[uu]; ++j) {
                uint32_t vv = workspace_->candidate_set_pointer_[uu][j];
                local_degree_sum += gv->get_neighbor_num(vv);
            }
            if (local_degree_sum < min_degree_sum) {
                min_degree_sum = local_degree_sum;
                selected_bn = uu;
                selected_idx = i;
            }
        }
    }
    bn.erase(bn.begin() + selected_idx);
    return selected_bn;
}

void
LocalViewManager::set_adjacent_update_candidates_flag_for_subgraph(GlobalViewManager &gvm, uint32_t u,
                                                                   std::vector<uint32_t> &bn,
                                                                   uint32_t encoded_v0, uint32_t encoded_v1,
                                                                   const std::vector<int32_t>* new2old_sub_to_full) {
    bool is_flag_based = false;
    if (workspace_->candidate_set_size_[u] < 1024) {
        is_flag_based = true;
    }
    else {
        for (auto uu : bn) {
            if (workspace_->candidate_set_size_[u] < 32 * workspace_->candidate_set_size_[uu]) {
                is_flag_based = true;
                break;
            }
        }
    }
    if (is_flag_based) {
        updated_count_ = 0;
        for (uint32_t i = 0; i < workspace_->candidate_set_size_[u]; ++i) {
            uint32_t v = workspace_->candidate_set_pointer_[u][i];
            if (v != encoded_v0 && v != encoded_v1) {
                workspace_->flag_array_[v] = 1;
                workspace_->updated_flag_[updated_count_++] = v;
            }
        }
    }
    else {
        uint32_t uu = select_bn_with_minimum_degree_sum_for_subgraph(gvm, u, bn, new2old_sub_to_full);
        const uint32_t full_u = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, u);
        const uint32_t full_uu = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, uu);
        auto gv = gvm.get_nlf_view({full_uu, full_u});

        updated_count_ = 0;
        for (uint32_t i = 0; i < workspace_->candidate_set_size_[uu]; ++i) {
            uint32_t vv = workspace_->candidate_set_pointer_[uu][i];

            uint32_t vv_nbrs_count;
            auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);

            uint32_t lc = 0;
            ComputeSetIntersection::ComputeCandidates(workspace_->candidate_set_pointer_[u], workspace_->candidate_set_size_[u],
                                                      vv_nbrs, vv_nbrs_count, workspace_->si_buffer_.data(), lc);

            for (uint32_t j = 0; j < lc; ++j) {
                uint32_t v = workspace_->si_buffer_[j];

                if (workspace_->flag_array_[v] == 0 && v != encoded_v0 && v != encoded_v1) {
                    workspace_->flag_array_[v] = 1;
                    workspace_->updated_flag_[updated_count_++] = v;
                }
            }
            if (workspace_->candidate_set_size_[u] == updated_count_)
                break;
        }
    }
}

void LocalViewManager::set_update_candidates_flag_for_subgraph(GlobalViewManager &gvm, uint32_t u,
                                                               std::vector<uint32_t> &bn,
                                                               uint32_t encoded_v0, uint32_t encoded_v1,
                                                               const std::vector<int32_t>* new2old_sub_to_full) {
    uint32_t selected_bn = select_bn_with_minimum_degree_sum_for_subgraph(gvm, u, bn, new2old_sub_to_full);
    const uint32_t full_u = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, u);
    const uint32_t full_selected_bn = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, selected_bn);
    auto gv = gvm.get_nlf_view({full_selected_bn, full_u});
    updated_count_ = 0;
    for (uint32_t i = 0; i < workspace_->candidate_set_size_[selected_bn]; ++i) {
        uint32_t vv = workspace_->candidate_set_pointer_[selected_bn][i];
        uint32_t vv_nbrs_count;
        auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);

        for (uint32_t j = 0; j < vv_nbrs_count; ++j) {
            uint32_t v = vv_nbrs[j];
            if (workspace_->flag_array_[v] == 0 && v != encoded_v0 && v != encoded_v1) {
                workspace_->flag_array_[v] = 1;
                workspace_->updated_flag_[updated_count_++] = v;
            }
        }
    }
}

bool LocalViewManager::prune_local_candidates_for_subgraph(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                                           const std::vector<int32_t>* new2old_sub_to_full) {
    const uint32_t full_u = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, u);
    uint32_t flag_value = 1;
    uint32_t current_valid_count = updated_count_;
    for (auto uu : bn) {
        const uint32_t full_uu = map_sub_query_vertex_for_full_gvm(new2old_sub_to_full, uu);
        auto gv = gvm.get_nlf_view({full_uu, full_u});
        auto reverse_gv = gvm.get_nlf_view({full_u, full_uu});
        if (!gv || !reverse_gv) {
            return false;
        }

        uint32_t target_valid_count = current_valid_count;
        current_valid_count = 0;

        for (uint32_t j = 0; j < workspace_->candidate_set_size_[uu]; ++j) {
            uint32_t vv = workspace_->candidate_set_pointer_[uu][j];
            uint32_t vv_nbrs_count;
            auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);

            if (vv_nbrs_count < 1024 || vv_nbrs_count < updated_count_ * 32) {
                for (uint32_t k = 0; k < vv_nbrs_count; ++k) {
                    uint32_t v = vv_nbrs[k];

                    if (workspace_->flag_array_[v] == flag_value) {
                        workspace_->flag_array_[v] += 1;
                        current_valid_count += 1;
                    }
                }
            }
            else {
                for (uint32_t k = 0; k < updated_count_; ++k) {
                    uint32_t v = workspace_->updated_flag_[k];
                    if (workspace_->flag_array_[v] == flag_value) {
                        uint32_t v_nbrs_count;
                        auto v_nbrs = reverse_gv->get_neighbor(v, v_nbrs_count);
                        if (vv_nbrs_count < v_nbrs_count) {
                            auto it = std::lower_bound(vv_nbrs, vv_nbrs + vv_nbrs_count, v);
                            if (it != vv_nbrs + vv_nbrs_count && *it == v) {
                                workspace_->flag_array_[v] += 1;
                                current_valid_count += 1;
                            }
                        }
                        else {
                            auto it = std::lower_bound(v_nbrs, v_nbrs + v_nbrs_count, vv);
                            if (it != v_nbrs + v_nbrs_count && *it == vv) {
                                workspace_->flag_array_[v] += 1;
                                current_valid_count += 1;
                            }
                        }
                    }
                }
            }
            if (current_valid_count == target_valid_count)
                break;
        }
        flag_value += 1;
    }

    bool push = storage_.candidates_store_[u].empty();
    uint32_t local_pos = 0;
    for (uint32_t j = 0; j < updated_count_; ++j) {
        uint32_t v = workspace_->updated_flag_[j];
        if (workspace_->flag_array_[v] == flag_value) {
            if (push)
                storage_.candidates_store_[u].push_back(v);
            else
                storage_.candidates_store_[u][local_pos++] = v;
        }
        workspace_->flag_array_[v] = 0;
    }

    if (!push)
        storage_.candidates_store_[u].resize(local_pos);
    workspace_->candidate_set_pointer_[u] = storage_.candidates_store_[u].data();
    workspace_->candidate_set_size_[u] = storage_.candidates_store_[u].size();
    if (workspace_->candidate_set_size_[u] == 0) {
        return false;
    }
    return true;
}