
#include <chrono>
#include <functional>
#include <cassert>
#include "search_engine.h"
#include "computesetintersection.h"
#include "streaming_config.h"
#include <limits>

void SearchEngine::initialize(const Graph *query_graph, const Graph *data_graph) {
    uint32_t n = query_graph->getVerticesCount();
    uint32_t N = data_graph->getVerticesCount();
    embedding_.resize(n);
    encoded_embedding_.resize(n);
    local_idx_.resize(n);
    local_candidates_store_.resize(n);
    encoded_local_candidates_store_.resize(n);
    local_candidates_buffer1_.resize(n);
    for (auto& buffer : local_candidates_buffer1_) {
        buffer.resize(N);
    }
    local_candidates_buffer2_.resize(n);
    for (auto& buffer : local_candidates_buffer2_) {
        buffer.resize(N);
    }
    visited_ = new bool [data_graph->getVerticesCount()];
    std::fill(visited_, visited_ + data_graph->getVerticesCount(), false);
    reset_performance_counters();
    clear_phase3_block_cache();
}

void SearchEngine::release() {
    embedding_.clear();
    encoded_embedding_.clear();
    local_idx_.clear();
    local_candidates_store_.clear();
    encoded_local_candidates_store_.clear();
    local_candidates_buffer1_.clear();
    local_candidates_buffer2_.clear();
    delete []visited_;
    clear_phase3_block_cache();
}

uint64_t SearchEngine::search_on_reduced_query(const Graph *query_graph, OrdersPerEdge &orders, LocalViewManager &lvm,
                                               GlobalViewManager &gvm) {
    uint64_t result_count = 0;
    auto& order = orders.matching_order_;
    auto& bn_offset = orders.matching_order_bn_offset_;
    auto& bn = orders.matching_order_bn_;
    auto& view_mapping = orders.matching_order_view_mappings_;
    uint32_t target_depth = order.size();
    uint32_t start_vertex = order[0];
    Edge data_edge = {lvm.get_candidate_set(orders.indexing_order_[0])[0], lvm.get_candidate_set(orders.indexing_order_[1])[0]};
    embedding_[orders.indexing_order_[0]] = data_edge.first;
    embedding_[orders.indexing_order_[1]] = data_edge.second;
    visited_[data_edge.first] = true;
    visited_[data_edge.second] = true;
    uint32_t* seeds;
    local_idx_[0].first = 0;
    if (query_graph->getVertexDegree(start_vertex) == 1) {
        local_idx_[0].second = compute_local_candidates_for_reduced_query(query_graph, 0, order,
                                                                                      bn_offset, bn,
                                                                                      view_mapping, lvm, gvm);
        seeds = local_candidates_store_[0];
    }
    else {
        seeds = lvm.get_candidate_set(start_vertex);
        local_idx_[0].second = lvm.get_candidate_set_size(start_vertex);
    }
    if (target_depth == 1) {
#if EXECUTION_MODE == 1
        for (local_idx_[0].first = 0; local_idx_[0].first < local_idx_[0].second; ++local_idx_[0].first) {
            uint32_t v = (*seeds)[local_idx_[0].first];
            embedding_[start_vertex] = v;
        }
#endif
        visited_[data_edge.first] = false;
        visited_[data_edge.second] = false;
        return local_idx_[0].second;
    }
#ifdef ENABLE_PERFORMANCE_COUNTERS
    std::vector<uint64_t> result_count_vec(target_depth, 0);
#endif
    for (local_idx_[0].first = 0; local_idx_[0].first < local_idx_[0].second; local_idx_[0].first++) {
#ifdef ENABLE_PERFORMANCE_COUNTERS
        partial_result_count_ += 1;
        result_count_vec[0] = result_count;
#endif
        uint32_t seed = seeds[local_idx_[0].first];
        if (visited_[seed])
            continue;
        embedding_[start_vertex] = seed;
        encoded_embedding_[start_vertex] = local_idx_[0].first;
        visited_[seed] = true;
        uint32_t current_depth = 1;
        local_idx_[current_depth].first = 0;
        local_idx_[current_depth].second = compute_local_candidates_for_reduced_query(query_graph, current_depth, order,
                                                                                      bn_offset, bn,
                                                                                      view_mapping, lvm, gvm);                   
        if (target_depth == 2) {
#if EXECUTION_MODE == 1
            for (;local_idx_[current_depth].first < local_idx_[current_depth].second; local_idx_[current_depth].first++) {
                uint32_t vv = local_candidates_store_[current_depth][local_idx_[current_depth].first];
                embedding_[last_vertex] = vv;
            }
#endif
            result_count += local_idx_[current_depth].second;
            visited_[seed] = false;

            if (result_count >= target_number)
                return result_count;
        } else {
            while (true) {
                while (local_idx_[current_depth].first < local_idx_[current_depth].second) {
#ifdef ENABLE_PERFORMANCE_COUNTERS
                    partial_result_count_ += 1;
                    result_count_vec[current_depth] = result_count;
#endif
                    uint32_t u = order[current_depth];
                    uint32_t encoded_v = encoded_local_candidates_store_[current_depth][local_idx_[current_depth].first];
                    uint32_t v = local_candidates_store_[current_depth][local_idx_[current_depth].first++];
                    encoded_embedding_[u] = encoded_v;
                    embedding_[u] = v;
                    visited_[v] = true;

                    uint32_t next_depth = current_depth + 1;
                    local_idx_[next_depth].first = 0;
                    local_idx_[next_depth].second = compute_local_candidates_for_reduced_query(query_graph, next_depth,
                                                                                               order,
                                                                                               bn_offset, bn,
                                                                                               view_mapping, lvm,
                                                                                               gvm);
                    if (local_idx_[next_depth].second == 0) {
#ifdef ENABLE_PERFORMANCE_COUNTERS
                        lc_empty_count_ += 1;
                        invalid_partial_result_count_ += 1;
#endif
                        visited_[v] = false;
                    } else if (next_depth == target_depth - 1) {
                        result_count += local_idx_[next_depth].second;
                        if (result_count >= target_number) {
                            for (uint32_t x = 0; x < target_depth; ++x) {
                                visited_[embedding_[order[x]]] = false;
                            }
                            visited_[data_edge.first] = false;
                            visited_[data_edge.second] = false;
                            return result_count;
                        }
                        visited_[v] = false;

                        if(g_exit) {
                            return result_count;
                        }
                    } else {
                        current_depth += 1;
                    }
                }
                if (g_exit) {
                    return result_count;
                }

                current_depth -= 1;
                visited_[embedding_[order[current_depth]]] = false;
#ifdef ENABLE_PERFORMANCE_COUNTERS
                if (result_count == result_count_vec[current_depth]) {
                    invalid_partial_result_count_ += 1;
                }
#endif
                if (current_depth < 1) {
                    break;
                }
            }
        }
    }

    visited_[data_edge.first] = false;
    visited_[data_edge.second] = false;
    return result_count;
}

uint64_t SearchEngine::block_reuse_search_on_reduced_query(const Graph *query_graph,
                                                           OrdersPerEdge &orders,
                                                           LocalViewManager &lvm,
                                                           GlobalViewManager &gvm,
                                                           const QueryBCTreeDecomposition& query_bctree) {
    phase3_block_reuse_entry_count_ += 1;
    auto& cache = phase3_block_cache_;
    const uint32_t INVALID = std::numeric_limits<uint32_t>::max();
    uint64_t result_count = 0;
    auto& order = orders.matching_order_;
    auto& bn_offset = orders.matching_order_bn_offset_;
    auto& bn = orders.matching_order_bn_;
    auto& view_mapping = orders.matching_order_view_mappings_;
    const auto& blocks = query_bctree.blocks;
    uint32_t target_depth = order.size();

    std::vector<char> excluded_suffix_block(query_bctree.block_count, 0);
    const uint32_t anchor_qv0 = orders.indexing_order_[0];
    const uint32_t anchor_qv1 = orders.indexing_order_[1];
    if (anchor_qv0 < blocks.size() && !blocks[anchor_qv0].empty()) {
        excluded_suffix_block[blocks[anchor_qv0][0]] = 1;
    }
    if (anchor_qv1 < blocks.size() && !blocks[anchor_qv1].empty()) {
        excluded_suffix_block[blocks[anchor_qv1][0]] = 1;
    }
    std::vector<uint32_t> block_inner_count_by_id(query_bctree.block_count, 0);
    uint32_t suffix_inner_total = 0;
    for (uint32_t block_id = 0; block_id < query_bctree.block_count; ++block_id) {
        if (excluded_suffix_block[block_id]) {
            continue;
        }
        uint32_t inner_cnt = query_bctree.blocks[block_id].size() - 1;
        block_inner_count_by_id[block_id] = inner_cnt;
        suffix_inner_total += inner_cnt;
    }
    uint32_t target_base_depth = target_depth - suffix_inner_total;
    std::vector<uint32_t> block_start_depth_by_id(query_bctree.block_count, target_depth);
    uint32_t suffix_cursor = target_base_depth;
    for (uint32_t block_id = 0; block_id < query_bctree.block_count; ++block_id) {
        if (excluded_suffix_block[block_id]) {
            continue;
        }
        block_start_depth_by_id[block_id] = suffix_cursor;
        suffix_cursor += block_inner_count_by_id[block_id];
    }

    auto get_block_count_on_suffix = [&](uint32_t block_id) -> uint64_t {
        if (block_id >= block_inner_count_by_id.size()) {
            return 1;
        }
        uint32_t inner_count = block_inner_count_by_id[block_id];
        if (inner_count == 0) {
            return 1;
        }
        uint32_t begin_depth = block_start_depth_by_id[block_id];
        uint32_t end_depth = begin_depth + inner_count;
        uint64_t block_result_count = 0;
        uint32_t current_depth_2 = begin_depth;
        local_idx_[current_depth_2].first = 0;
        local_idx_[current_depth_2].second = compute_local_candidates_for_reduced_query(query_graph, current_depth_2, order,
                                                                                         bn_offset, bn,
                                                                                         view_mapping, lvm, gvm);
        while (true) {
            while (local_idx_[current_depth_2].first < local_idx_[current_depth_2].second) {
                uint32_t u_blk = order[current_depth_2];
                uint32_t encoded_v = encoded_local_candidates_store_[current_depth_2][local_idx_[current_depth_2].first];
                uint32_t v_blk = local_candidates_store_[current_depth_2][local_idx_[current_depth_2].first++];
                if (visited_[v_blk]) {
                    continue;
                }
                encoded_embedding_[u_blk] = encoded_v;
                embedding_[u_blk] = v_blk;
                visited_[v_blk] = true;
                uint32_t next_depth_2 = current_depth_2 + 1;
                if (next_depth_2 >= end_depth) {
                    block_result_count += 1;
                    visited_[v_blk] = false;
                    embedding_[u_blk] = INVALID;
                    encoded_embedding_[u_blk] = INVALID;
                    if (g_exit) {
                        return block_result_count;
                    }
                } else {
                    local_idx_[next_depth_2].first = 0;
                    local_idx_[next_depth_2].second = compute_local_candidates_for_reduced_query(query_graph, next_depth_2, order,
                                                                                                  bn_offset, bn,
                                                                                                  view_mapping, lvm, gvm);
                    if (local_idx_[next_depth_2].second == 0) {
                        visited_[v_blk] = false;
                        embedding_[u_blk] = INVALID;
                        encoded_embedding_[u_blk] = INVALID;
                    } else if (next_depth_2 == end_depth - 1) {
                        block_result_count += local_idx_[next_depth_2].second;
                        visited_[v_blk] = false;
                        embedding_[u_blk] = INVALID;
                        encoded_embedding_[u_blk] = INVALID;
                        if (g_exit) {
                            return block_result_count;
                        }
                    } else {
                        current_depth_2 += 1;
                    }
                }
            }
            if (g_exit) {
                return block_result_count;
            }
            if (current_depth_2 <= begin_depth) {
                break;
            }
            current_depth_2 -= 1;
            uint32_t qv_back = order[current_depth_2];
            if (qv_back < embedding_.size() && embedding_[qv_back] != INVALID) {
                visited_[embedding_[qv_back]] = false;
                embedding_[qv_back] = INVALID;
                encoded_embedding_[qv_back] = INVALID;
            }
        }
        return block_result_count;
    };

    uint32_t start_vertex = order[0];
    if (target_base_depth < 2) {
        phase3_block_reuse_fallback_plain_count_ += 1;
        return search_on_reduced_query(query_graph, orders, lvm, gvm);
    }
    phase3_block_reuse_skeleton_run_count_ += 1;
    Edge data_edge = {lvm.get_candidate_set(orders.indexing_order_[0])[0], lvm.get_candidate_set(orders.indexing_order_[1])[0]};
    embedding_[orders.indexing_order_[0]] = data_edge.first;
    embedding_[orders.indexing_order_[1]] = data_edge.second;
    visited_[data_edge.first] = true;
    visited_[data_edge.second] = true;
    uint32_t* seeds;
    local_idx_[0].first = 0;
    if (query_graph->getVertexDegree(start_vertex) == 1) {
        local_idx_[0].second = compute_local_candidates_for_reduced_query(query_graph, 0, order,
                                                                          bn_offset, bn,
                                                                          view_mapping, lvm, gvm);
        seeds = local_candidates_store_[0];
    } else {
        seeds = lvm.get_candidate_set(start_vertex);
        local_idx_[0].second = lvm.get_candidate_set_size(start_vertex);
    }
    if (target_depth == 1) {
        visited_[data_edge.first] = false;
        visited_[data_edge.second] = false;
        return local_idx_[0].second;
    }
    std::vector<uint64_t> base_path_multiplier(target_base_depth, 1);
    for (local_idx_[0].first = 0; local_idx_[0].first < local_idx_[0].second; local_idx_[0].first++) {


        uint32_t seed = seeds[local_idx_[0].first];
        if (visited_[seed])
            continue;

        embedding_[start_vertex] = seed;
        encoded_embedding_[start_vertex] = local_idx_[0].first;

        visited_[seed] = true;

        uint64_t seed_block_mult = 1;
        if (query_bctree.is_cut_vertex[start_vertex]) {
            const uint32_t block_id = blocks[start_vertex][0];
            uint32_t cached_mult = 0;
            bool hit = false;
            auto outer_it = cache.find(block_id);
            if (outer_it != cache.end()) {
                auto inner_it = outer_it->second.find(seed);
                if (inner_it != outer_it->second.end()) {
                    hit = true;
                    cached_mult = inner_it->second;
                }
            }
            if (hit) {
                seed_block_mult = cached_mult;
            } else {
                seed_block_mult = get_block_count_on_suffix(block_id);
                cache[block_id][seed] = seed_block_mult;
            }
            if (seed_block_mult == 0) {
                visited_[seed] = false;
                continue;
            }
        }

        uint32_t current_depth = 1;
        uint32_t current_depth_base_order = current_depth;
        base_path_multiplier[0] = seed_block_mult;
        if (target_base_depth > 1) {
            base_path_multiplier[1] = seed_block_mult;
        }
        local_idx_[current_depth].first = 0;
        local_idx_[current_depth].second = compute_local_candidates_for_reduced_query(query_graph, current_depth, order,
                                                                                      bn_offset, bn,
                                                                                      view_mapping, lvm, gvm);
                                      
        if (target_depth == 2) {
            result_count += local_idx_[current_depth].second;
            visited_[seed] = false;

            if (result_count >= target_number)
                return result_count;
        } else {
            while (true) {
        
                current_depth = current_depth_base_order;
                if (current_depth >= local_idx_.size()) {
                    phase3_skeleton_invalid_depth_break_count_ += 1;
                    break;
                }
                while (local_idx_[current_depth].first < local_idx_[current_depth].second) {

                    uint32_t u = order[current_depth];
                    uint32_t encoded_v = encoded_local_candidates_store_[current_depth][local_idx_[current_depth].first];
                    uint32_t v = local_candidates_store_[current_depth][local_idx_[current_depth].first++];
                    encoded_embedding_[u] = encoded_v;
                    embedding_[u] = v;
                    visited_[v] = true;
                    uint64_t branch_multiplier = base_path_multiplier[current_depth_base_order];
                    if(query_bctree.is_cut_vertex[u]){
                        phase3_cut_vertex_branch_count_ += 1;
                        if (u >= blocks.size() || blocks[u].empty()) {
                            phase3_cut_skip_no_phase3_block_count_ += 1;
                            visited_[v] = false;
                            continue;
                        }
                        const uint32_t block_id = blocks[u][0];
                        const uint32_t entry_data_v = v;
                        auto outer_it = cache.find(block_id);
                        uint64_t cached_block_count = 0;
                        bool phase3_block_cache_hit = false;
                        if (outer_it != cache.end()) {
                            auto inner_it = outer_it->second.find(entry_data_v);
                            if (inner_it != outer_it->second.end()) {
                                phase3_block_cache_hit = true;
                                cached_block_count = inner_it->second;
                            }
                        }
                        if (phase3_block_cache_hit) {
                            phase3_block_cache_hit_count_ += 1;
                            branch_multiplier *= cached_block_count;
                        } else {
                            phase3_block_cache_miss_count_ += 1;
                            phase3_block_inner_dfs_run_count_ += 1;
                            uint64_t block_result_count = get_block_count_on_suffix(block_id);
                            cache[block_id][entry_data_v] = block_result_count;
                            phase3_block_cache_store_count_ += 1;
                            branch_multiplier *= block_result_count;
                        }
                        if (branch_multiplier == 0) {

                            visited_[v] = false;
                            continue;
                        }
                    }

                    uint32_t next_depth_base_order = current_depth_base_order + 1;
                    const uint32_t next_depth = next_depth_base_order;
                    if (next_depth >= local_idx_.size()) {
                        visited_[v] = false;
                        continue;
                    }

                    local_idx_[next_depth].first = 0;
                    local_idx_[next_depth].second = compute_local_candidates_for_reduced_query(query_graph, next_depth,
                                                                                               order,
                                                                                               bn_offset, bn,
                                                                                               view_mapping, lvm,
                                                                                               gvm);

                    if (local_idx_[next_depth].second == 0) {

                        visited_[v] = false;
                    } else if (next_depth_base_order == target_base_depth - 1) {
                        result_count += branch_multiplier * static_cast<uint64_t>(local_idx_[next_depth].second);

                        if (result_count >= target_number) {
                            for (uint32_t x = 0; x < target_base_depth; ++x) {
                                visited_[embedding_[order[x]]] = false;
                            }
                            visited_[data_edge.first] = false;
                            visited_[data_edge.second] = false;
                            return result_count;
                        }
                        visited_[v] = false;

                        if(g_exit) {
                            return result_count;
                        }
                    } else {
                        current_depth_base_order += 1;
                        if (current_depth_base_order < base_path_multiplier.size()) {
                            base_path_multiplier[current_depth_base_order] = branch_multiplier;
                        }
                        current_depth = current_depth_base_order;
                    }
                }

                if (g_exit) {
                    return result_count;
                }

                current_depth_base_order -= 1;
                visited_[embedding_[order[current_depth_base_order]]] = false;

                if (current_depth_base_order < 1) {
                    break;
                }
            }
        }
    }

    visited_[data_edge.first] = false;
    visited_[data_edge.second] = false;
    return result_count;
}

std::vector<std::vector<uint32_t>> SearchEngine::special_search_on_reduced_query(const Graph *query_graph,
                                                                                 OrdersPerEdge &orders,
                                                                                 LocalViewManager &lvm,
                                                                                 GlobalViewManager &gvm,
                                                                                 uint64_t max_results,
                                                                                 const std::vector<int32_t>* new2old) {
    std::vector<std::vector<uint32_t>> results;
    if (max_results == 0) return results;

    auto& order = orders.matching_order_;
    auto& bn_offset = orders.matching_order_bn_offset_;
    auto& bn = orders.matching_order_bn_;
    auto& view_mapping = orders.matching_order_view_mappings_;

    const uint32_t target_depth = order.size();
    const uint32_t start_vertex = order[0];

    const uint32_t INVALID = std::numeric_limits<uint32_t>::max();
    const uint32_t n = query_graph->getVerticesCount();
    std::fill(embedding_.begin(), embedding_.end(), INVALID);

    auto emit_embedding = [&]() {
        std::vector<uint32_t> emb;
        emb.reserve(n);
        for (uint32_t u = 0; u < n; ++u) emb.push_back(embedding_[u]);
        results.push_back(std::move(emb));
    };

    const uint32_t anchor_u0 = orders.indexing_order_[0];
    const uint32_t anchor_u1 = orders.indexing_order_[1];
    if (lvm.get_candidate_set_size(anchor_u0) == 0 || lvm.get_candidate_set_size(anchor_u1) == 0) {
        return results;
    }
    Edge data_edge = {lvm.get_candidate_set(anchor_u0)[0], lvm.get_candidate_set(anchor_u1)[0]};
    embedding_[anchor_u0] = data_edge.first;
    embedding_[anchor_u1] = data_edge.second;

    visited_[data_edge.first] = true;
    visited_[data_edge.second] = true;

    auto cleanup_and_return = [&]() {
        for (uint32_t x = 0; x < target_depth; ++x) {
            uint32_t mapped = embedding_[order[x]];
            if (mapped != INVALID) visited_[mapped] = false;
        }
        visited_[data_edge.first] = false;
        visited_[data_edge.second] = false;
        return results;
    };

    uint32_t* seeds;
    local_idx_[0].first = 0;
    if (query_graph->getVertexDegree(start_vertex) == 1) {
        local_idx_[0].second = compute_local_candidates_for_reduced_query_with_new2old(query_graph, 0, order, bn_offset, bn,
                                                                      view_mapping, lvm, gvm, *new2old);
        seeds = local_candidates_store_[0];
    } else {
        seeds = lvm.get_candidate_set(start_vertex);
        local_idx_[0].second = lvm.get_candidate_set_size(start_vertex);
    }

    if (target_depth == 1) {
        for (local_idx_[0].first = 0; local_idx_[0].first < local_idx_[0].second; ++local_idx_[0].first) {
            uint32_t v = seeds[local_idx_[0].first];
            if (visited_[v]) continue;
            embedding_[start_vertex] = v;
            emit_embedding();
            if (results.size() >= max_results) break;
        }
        visited_[data_edge.first] = false;
        visited_[data_edge.second] = false;
        return results;
    }

    for (local_idx_[0].first = 0; local_idx_[0].first < local_idx_[0].second; local_idx_[0].first++) {
        uint32_t seed = seeds[local_idx_[0].first];
        if (visited_[seed]) continue;

        embedding_[start_vertex] = seed;
        encoded_embedding_[start_vertex] = local_idx_[0].first;
        visited_[seed] = true;

        uint32_t current_depth = 1;
        local_idx_[current_depth].first = 0;
        local_idx_[current_depth].second = new2old
            ? compute_local_candidates_for_reduced_query_with_new2old(query_graph, current_depth, order, bn_offset, bn,
                                                                      view_mapping, lvm, gvm, *new2old)
            : compute_local_candidates_for_reduced_query(query_graph, current_depth, order, bn_offset, bn,
                                                         view_mapping, lvm, gvm);
        if (target_depth == 2) {
            const uint32_t last_u = order[1];
            for (; local_idx_[current_depth].first < local_idx_[current_depth].second; ++local_idx_[current_depth].first) {
                uint32_t vv = local_candidates_store_[current_depth][local_idx_[current_depth].first];
                embedding_[last_u] = vv;
                emit_embedding();
                if (results.size() >= max_results) break;
            }
            visited_[seed] = false;
            if (results.size() >= max_results) {
                visited_[data_edge.first] = false;
                visited_[data_edge.second] = false;
                return results;
            }
            continue;
        }

        while (true) {
            while (local_idx_[current_depth].first < local_idx_[current_depth].second) {
                uint32_t u = order[current_depth];
                uint32_t encoded_v = encoded_local_candidates_store_[current_depth][local_idx_[current_depth].first];
                uint32_t v = local_candidates_store_[current_depth][local_idx_[current_depth].first++];
                encoded_embedding_[u] = encoded_v;
                embedding_[u] = v;
                visited_[v] = true;

                uint32_t next_depth = current_depth + 1;
                local_idx_[next_depth].first = 0;
                local_idx_[next_depth].second = new2old
                    ? compute_local_candidates_for_reduced_query_with_new2old(query_graph, next_depth, order, bn_offset,
                                                                              bn, view_mapping, lvm, gvm, *new2old)
                    : compute_local_candidates_for_reduced_query(query_graph, next_depth, order, bn_offset, bn,
                                                                 view_mapping, lvm, gvm);

                if (local_idx_[next_depth].second == 0) {
                    visited_[v] = false;
                } else if (next_depth == target_depth - 1) {
                    const uint32_t last_u = order[next_depth];
                    for (uint32_t j = 0; j < local_idx_[next_depth].second; ++j) {
                        uint32_t vv = local_candidates_store_[next_depth][j];
                        embedding_[last_u] = vv;
                        emit_embedding();
                        if (results.size() >= max_results) {
                            visited_[v] = false;
                            return cleanup_and_return();
                        }
                    }
                    visited_[v] = false;
                    if (g_exit) return cleanup_and_return();
                } else {
                    current_depth += 1;
                }
            }

            if (g_exit) return cleanup_and_return();

            current_depth -= 1;
            uint32_t mapped = embedding_[order[current_depth]];
            if (mapped != INVALID) visited_[mapped] = false;
            if (current_depth < 1) break;
        }

        visited_[seed] = false;
        if (results.size() >= max_results) break;
    }
    visited_[data_edge.first] = false;
    visited_[data_edge.second] = false;
    return results;
}



uint64_t SearchEngine::constrained_search_on_reduced_query(const Graph *query_graph,
                                                           OrdersPerEdge &orders,
                                                           LocalViewManager &lvm,
                                                           GlobalViewManager &gvm,
                                                           const std::vector<std::pair<uint32_t, uint32_t>> &forced_assignments) {
    uint64_t result_count = 0;

    auto& order = orders.matching_order_;
    auto& bn_offset = orders.matching_order_bn_offset_;
    auto& bn = orders.matching_order_bn_;
    auto& view_mapping = orders.matching_order_view_mappings_;

    if (order.empty()) return 0;

    const uint32_t INVALID = std::numeric_limits<uint32_t>::max();
    const uint32_t n = query_graph->getVerticesCount();
    std::fill(embedding_.begin(), embedding_.end(), INVALID);

    std::vector<uint32_t> forced_data_by_query(n, INVALID);
    for (const auto &kv : forced_assignments) {
        if (kv.first >= n) return 0;
        if (forced_data_by_query[kv.first] != INVALID && forced_data_by_query[kv.first] != kv.second) {
            return 0;
        }
        forced_data_by_query[kv.first] = kv.second;
    }

    uint32_t target_depth = order.size();
    uint32_t start_vertex = order[0];

    Edge data_edge = {lvm.get_candidate_set(orders.indexing_order_[0])[0], lvm.get_candidate_set(orders.indexing_order_[1])[0]};
    embedding_[orders.indexing_order_[0]] = data_edge.first;
    embedding_[orders.indexing_order_[1]] = data_edge.second;
    uint32_t forced_idx0 = forced_data_by_query[orders.indexing_order_[0]];
    uint32_t forced_idx1 = forced_data_by_query[orders.indexing_order_[1]];
    if ((forced_idx0 != INVALID && forced_idx0 != data_edge.first) ||
        (forced_idx1 != INVALID && forced_idx1 != data_edge.second)) {
        return 0;
    }

    visited_[data_edge.first] = true;
    visited_[data_edge.second] = true;
    uint32_t* seeds;
    local_idx_[0].first = 0;
    local_idx_[0].second = compute_local_candidates_for_reduced_query(query_graph, 0, order,
                                                                       bn_offset, bn,
                                                                       view_mapping, lvm, gvm);
    seeds = local_candidates_store_[0];

    if (target_depth == 1) {
        for (local_idx_[0].first = 0; local_idx_[0].first < local_idx_[0].second; ++local_idx_[0].first) {
            uint32_t v = seeds[local_idx_[0].first];
            if (visited_[v]) continue;
            if (forced_data_by_query[start_vertex] != INVALID && forced_data_by_query[start_vertex] != v) continue;
            embedding_[start_vertex] = v;
            result_count += 1;
        }
        visited_[data_edge.first] = false;
        visited_[data_edge.second] = false;
        return result_count;
    }

    for (local_idx_[0].first = 0; local_idx_[0].first < local_idx_[0].second; local_idx_[0].first++) {
        uint32_t seed = seeds[local_idx_[0].first];
        if (visited_[seed]) continue;
        if (forced_data_by_query[start_vertex] != INVALID && forced_data_by_query[start_vertex] != seed) continue;

        embedding_[start_vertex] = seed;
        encoded_embedding_[start_vertex] = local_idx_[0].first;
        visited_[seed] = true;

        uint32_t current_depth = 1;
        local_idx_[current_depth].first = 0;
        local_idx_[current_depth].second = compute_local_candidates_for_reduced_query(query_graph, current_depth, order,
                                                                                       bn_offset, bn,
                                                                                       view_mapping, lvm, gvm);
        if (target_depth == 2) {
            uint32_t last_u = order[current_depth];
            for (; local_idx_[current_depth].first < local_idx_[current_depth].second; local_idx_[current_depth].first++) {
                uint32_t vv = local_candidates_store_[current_depth][local_idx_[current_depth].first];
                if (forced_data_by_query[last_u] != INVALID && forced_data_by_query[last_u] != vv) continue;
                embedding_[last_u] = vv;
                result_count += 1;
            }
            visited_[seed] = false;
            if (result_count >= target_number) {
                visited_[data_edge.first] = false;
                visited_[data_edge.second] = false;
                return result_count;
            }
        } else {
            while (true) {
                while (local_idx_[current_depth].first < local_idx_[current_depth].second) {
#ifdef ENABLE_PERFORMANCE_COUNTERS
                    partial_result_count_ += 1;
#endif
                    uint32_t u = order[current_depth];
                    uint32_t encoded_v = encoded_local_candidates_store_[current_depth][local_idx_[current_depth].first];
                    uint32_t v = local_candidates_store_[current_depth][local_idx_[current_depth].first++];

                    if (forced_data_by_query[u] != INVALID && forced_data_by_query[u] != v) {
                        continue;
                    }

                    encoded_embedding_[u] = encoded_v;
                    embedding_[u] = v;
                    visited_[v] = true;

                    uint32_t next_depth = current_depth + 1;
                    local_idx_[next_depth].first = 0;
                    local_idx_[next_depth].second = compute_local_candidates_for_reduced_query(query_graph, next_depth,
                                                                                                order,
                                                                                                bn_offset, bn,
                                                                                                view_mapping, lvm,
                                                                                                gvm);
                    if (local_idx_[next_depth].second == 0) {
#ifdef ENABLE_PERFORMANCE_COUNTERS
                        lc_empty_count_ += 1;
                        invalid_partial_result_count_ += 1;
#endif
                        visited_[v] = false;
                    } else if (next_depth == target_depth - 1) {
                        uint32_t last_u = order[next_depth];
                        for (uint32_t j = 0; j < local_idx_[next_depth].second; ++j) {
                            uint32_t vv = local_candidates_store_[next_depth][j];
                            if (forced_data_by_query[last_u] != INVALID && forced_data_by_query[last_u] != vv) continue;
                            embedding_[last_u] = vv;
                            result_count += 1;
                            if (result_count >= target_number) {
                                for (uint32_t x = 0; x < target_depth; ++x) {
                                    uint32_t mapped = embedding_[order[x]];
                                    if (mapped != INVALID) visited_[mapped] = false;
                                }
                                visited_[data_edge.first] = false;
                                visited_[data_edge.second] = false;
                                return result_count;
                            }
                        }
                        visited_[v] = false;
                        if (g_exit) {
                            visited_[data_edge.first] = false;
                            visited_[data_edge.second] = false;
                            return result_count;
                        }
                    } else {
                        current_depth += 1;
                    }
                }

                if (g_exit) {
                    visited_[data_edge.first] = false;
                    visited_[data_edge.second] = false;
                    return result_count;
                }

                current_depth -= 1;
                uint32_t mapped = embedding_[order[current_depth]];
                if (mapped != INVALID) visited_[mapped] = false;
#ifdef ENABLE_PERFORMANCE_COUNTERS
                invalid_partial_result_count_ += 1;
#endif
                if (current_depth < 1) {
                    break;
                }
            }
        }
    }

    visited_[data_edge.first] = false;
    visited_[data_edge.second] = false;
    return result_count;
}

uint64_t SearchEngine::constrained_search_with_subemb(const Graph *query_graph,
                                                      OrdersPerEdge &orders,
                                                      LocalViewManager &lvm,
                                                      GlobalViewManager &gvm,
                                                      const std::vector<std::pair<uint32_t, uint32_t>> &forced_pairs_template,
                                                      const std::vector<uint32_t> &emb) {
    uint64_t result_count = 0;

    auto& order = orders.matching_order_;
    auto& bn_offset = orders.matching_order_bn_offset_;
    auto& bn = orders.matching_order_bn_;
    auto& view_mapping = orders.matching_order_view_mappings_;

    if (order.empty()) return 0;
    for (const auto &p : forced_pairs_template) {
        const uint32_t data_v = emb[p.second];
        embedding_[p.first] = data_v;
        visited_[data_v] = true;
        const uint32_t u = p.first;
        const uint32_t cand_sz = lvm.get_candidate_set_size(u);
        uint32_t *cand = lvm.get_candidate_set(u);
        const bool non_leaf_u = (query_graph->getVertexDegree(u) > 1);
        const uint32_t enc_cap = gvm.get_encoded_candidate_set_size(u);
        bool found = false;
        for (uint32_t i = 0; i < cand_sz; ++i) {
            if (cand[i] == data_v) {
                encoded_embedding_[u] = i;
                found = true;
                break;
            }
            if (non_leaf_u && cand[i] < enc_cap && gvm.get_id(u, cand[i]) == data_v) {
                encoded_embedding_[u] = i;
                found = true;
                break;
            }
        }
        if (!found) {
            for (const auto &pp : forced_pairs_template) {
                const uint32_t vv = emb[pp.second];
                visited_[vv] = false;
            }
            return 0;
        }
    }

    uint32_t target_depth = order.size();
    uint32_t start_depth = static_cast<uint32_t>(forced_pairs_template.size() - 2);
    if (start_depth > target_depth) start_depth = target_depth;
    uint32_t start_vertex = order[start_depth];
    const uint32_t effective_depth = target_depth - start_depth;
    uint32_t* seeds;
    local_idx_[start_depth].first = 0;
    local_idx_[start_depth].second = compute_local_candidates_for_reduced_query(query_graph, start_depth, order,
                                                                                bn_offset, bn,
                                                                                view_mapping, lvm, gvm);
    seeds = local_candidates_store_[start_depth];
    if (effective_depth == 1) {
        for (local_idx_[start_depth].first = 0; local_idx_[start_depth].first < local_idx_[start_depth].second; ++local_idx_[start_depth].first) {
            uint32_t v = seeds[local_idx_[start_depth].first];
            if (visited_[v]) continue;
            embedding_[start_vertex] = v;
            result_count += 1;
        }
        for (const auto &p : forced_pairs_template) {
            const uint32_t data_v = emb[p.second];
            visited_[data_v] = false;
        }
        return result_count;
    }
    for (local_idx_[start_depth].first = 0; local_idx_[start_depth].first < local_idx_[start_depth].second; local_idx_[start_depth].first++) {
        uint32_t seed = seeds[local_idx_[start_depth].first];
        if (visited_[seed]) continue;

        embedding_[start_vertex] = seed;
        encoded_embedding_[start_vertex] = local_idx_[start_depth].first;
        visited_[seed] = true;

        uint32_t current_depth = start_depth + 1;
        local_idx_[current_depth].first = 0;
        local_idx_[current_depth].second = compute_local_candidates_for_reduced_query(query_graph, current_depth, order,
                                                                                       bn_offset, bn,
                                                                                       view_mapping, lvm, gvm);
        if (effective_depth == 2) {
            uint32_t last_u = order[current_depth];
            for (; local_idx_[current_depth].first < local_idx_[current_depth].second; local_idx_[current_depth].first++) {
                uint32_t vv = local_candidates_store_[current_depth][local_idx_[current_depth].first];
                embedding_[last_u] = vv;
                result_count += 1;
            }
            visited_[seed] = false;
            if (result_count >= target_number) {
                for (const auto &p : forced_pairs_template) {
                    const uint32_t data_v = emb[p.second];
                    visited_[data_v] = false;
                }
                return result_count;
            }
        } else {
            while (true) {
                while (local_idx_[current_depth].first < local_idx_[current_depth].second) {
#ifdef ENABLE_PERFORMANCE_COUNTERS
                    partial_result_count_ += 1;
#endif
                    uint32_t u = order[current_depth];
                    uint32_t encoded_v = encoded_local_candidates_store_[current_depth][local_idx_[current_depth].first];
                    uint32_t v = local_candidates_store_[current_depth][local_idx_[current_depth].first++];


                    encoded_embedding_[u] = encoded_v;
                    embedding_[u] = v;
                    visited_[v] = true;

                    uint32_t next_depth = current_depth + 1;
                    local_idx_[next_depth].first = 0;
                    local_idx_[next_depth].second = compute_local_candidates_for_reduced_query(query_graph, next_depth,
                                                                                                order,
                                                                                                bn_offset, bn,
                                                                                                view_mapping, lvm,
                                                                                                gvm);
                    if (local_idx_[next_depth].second == 0) {
#ifdef ENABLE_PERFORMANCE_COUNTERS
                        lc_empty_count_ += 1;
                        invalid_partial_result_count_ += 1;
#endif
                        visited_[v] = false;
                    } else if (next_depth == target_depth - 1) {
                        uint32_t last_u = order[next_depth];
                        for (uint32_t j = 0; j < local_idx_[next_depth].second; ++j) {
                            uint32_t vv = local_candidates_store_[next_depth][j];
                            embedding_[last_u] = vv;
                            result_count += 1;
                            if (result_count >= target_number) {
                                for (uint32_t x = start_depth; x <= current_depth; ++x) {
                                    visited_[embedding_[order[x]]] = false;
                                }
                                for (const auto &p : forced_pairs_template) {
                                    const uint32_t data_v = emb[p.second];
                                    visited_[data_v] = false;
                                }
                                return result_count;
                            }
                        }
                        visited_[v] = false;
                        if (g_exit) {
                            for (uint32_t x = start_depth; x < current_depth; ++x) {
                                visited_[embedding_[order[x]]] = false;
                            }
                            for (const auto &p : forced_pairs_template) {
                                const uint32_t data_v = emb[p.second];
                                visited_[data_v] = false;
                            }
                            return result_count;
                        }
                    } else {
                        current_depth += 1;
                    }
                }

                if (g_exit) {
                    for (uint32_t x = start_depth; x <= current_depth; ++x) {
                        visited_[embedding_[order[x]]] = false;
                    }
                    for (const auto &p : forced_pairs_template) {
                        const uint32_t data_v = emb[p.second];
                        visited_[data_v] = false;
                    }
                    return result_count;
                }

                current_depth -= 1;
                visited_[embedding_[order[current_depth]]] = false;
#ifdef ENABLE_PERFORMANCE_COUNTERS
                invalid_partial_result_count_ += 1;
#endif
                if (current_depth < start_depth + 1) {
                    break;
                }
            }
        }
    }

    for (const auto &p : forced_pairs_template) {
        const uint32_t data_v = emb[p.second];
        visited_[data_v] = false;
    }
    return result_count;
}

uint32_t SearchEngine::compute_local_candidates_for_reduced_query(const Graph *query_graph, uint32_t depth,
                                                                  std::vector<uint32_t> &order,
                                                                  std::vector<uint32_t> &bn_offset,
                                                                  std::vector<uint32_t> &bn,
                                                                  std::vector<uint32_t> &view_mapping,
                                                                  LocalViewManager &lvm, GlobalViewManager &gvm) {
    uint32_t bn_begin = bn_offset[depth];
    uint32_t bn_end = bn_offset[depth + 1];
    bool one_bn = bn_end - bn_begin == 1;
    uint32_t *lc1 = nullptr;
    uint32_t lc_count1 = 0;
    uint32_t *lc2 = nullptr;
    uint32_t lc_count2 = 0;

    if (bn_end - bn_begin == 0) {
        lc1 = lvm.get_candidate_set(order[depth]);
        lc_count1 = lvm.get_candidate_set_size(order[depth]);
        lc2 = local_candidates_buffer1_[depth].data();

        lc_count2 = 0;

        for (uint32_t i = 0; i < lc_count1; ++i) {
            uint32_t v = lc1[i];

            if (!visited_[v]) {
                local_candidates_buffer2_[depth][lc_count2] = i;
                lc2[lc_count2++] = v;
            }
#ifdef ENABLE_PERFORMANCE_COUNTERS
            else {
                iso_conflict_count_ += 1;
            }
#endif
        }
        encoded_local_candidates_store_[depth] = local_candidates_buffer2_[depth].data();
        local_candidates_store_[depth] = lc2;
        return lc_count2;
    }
    else {
        uint32_t u = bn[bn_begin];
        uint32_t v = embedding_[u];
        uint32_t encoded_v = encoded_embedding_[u];

        uint32_t view_id = view_mapping[bn_begin];

        lc1 = query_graph->getVertexDegree(order[depth]) == 1 ? gvm.get_view(view_id)->get_neighbor(v, lc_count1)
                : lvm.get_view(view_id)->get_neighbors(encoded_v, lc_count1);

        lc2 = local_candidates_buffer1_[depth].data();

        if (lc_count1 == 0) {
            goto EXIT;
        }
        if (bn_begin + 1 < bn_end) {
            bn_begin += 1;
            u = bn[bn_begin];
            encoded_v = encoded_embedding_[u];
            view_id = view_mapping[bn_begin];

            lc2 = lvm.get_view(view_id)->get_neighbors(encoded_v,lc_count2);

            uint32_t temp_count;
            uint32_t *temp_buffer = local_candidates_buffer1_[depth].data();

            ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, lc2, lc_count2, temp_buffer, temp_count);

            if (temp_count == 0) {
                lc_count1 = 0;
                goto EXIT;
            }

            lc1 = temp_buffer;
            lc_count1 = temp_count;
            temp_buffer = local_candidates_buffer2_[depth].data();

            for (bn_begin += 1; bn_begin < bn_end; ++bn_begin) {
                u = bn[bn_begin];
                encoded_v = encoded_embedding_[u];
                view_id = view_mapping[bn_begin];

                lc2 =  lvm.get_view(view_id)->get_neighbors(encoded_v,lc_count2);
                ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, lc2, lc_count2, temp_buffer, temp_count);

                if (temp_count == 0) return 0;
                std::swap(temp_buffer, lc1);
                std::swap(temp_count, lc_count1);
            }

            lc2 = temp_buffer;
        }
    }

EXIT:
    if (lc_count1 == 0) {
#ifdef ENABLE_PERFORMANCE_COUNTERS
       si_empty_count_ += 1;
#endif
        return 0;
    }

    if (query_graph->getVertexDegree(order[depth]) == 1) {
        lc_count2 = 0;

        for (uint32_t i = 0; i < lc_count1; ++i) {
            uint32_t v = lc1[i];

            if (!visited_[v]) {
                lc2[lc_count2++] = v;
            }
#ifdef ENABLE_PERFORMANCE_COUNTERS
            else {
                iso_conflict_count_ += 1;
            }
#endif
        }
        encoded_local_candidates_store_[depth] = local_candidates_buffer2_[depth].data();
        local_candidates_store_[depth] = lc2;
        return lc_count2;
    }
    else {
        lc_count2 = 0;
        uint32_t* candidate_set = lvm.get_candidate_set(order[depth]);
        uint32_t* encoded_buffer = one_bn ? local_candidates_buffer2_[depth].data() : lc1;
        for (uint32_t i = 0; i < lc_count1; ++i) {
            uint32_t encoded_v = lc1[i];
            uint32_t v = candidate_set[encoded_v];

            if (!visited_[v]) {
                encoded_buffer[lc_count2] = encoded_v;
                lc2[lc_count2++] = v;
            }
#ifdef ENABLE_PERFORMANCE_COUNTERS
            else {
                iso_conflict_count_ += 1;
            }
#endif
        }

        encoded_local_candidates_store_[depth] = encoded_buffer;
        local_candidates_store_[depth] = lc2;
        return lc_count2;
    }
}

uint32_t SearchEngine::compute_local_candidates_for_reduced_query_with_new2old(const Graph *query_graph, uint32_t depth,
                                                                              std::vector<uint32_t> &order,
                                                                              std::vector<uint32_t> &bn_offset,
                                                                              std::vector<uint32_t> &bn,
                                                                              std::vector<uint32_t> &view_mapping,
                                                                              LocalViewManager &lvm, GlobalViewManager &gvm,
                                                                              const std::vector<int32_t> &new2old) {
    uint32_t bn_begin = bn_offset[depth];
    uint32_t bn_end = bn_offset[depth + 1];
    bool one_bn = bn_end - bn_begin == 1;
    uint32_t *lc1 = nullptr;
    uint32_t lc_count1 = 0;
    uint32_t *lc2 = nullptr;
    uint32_t lc_count2 = 0;

    if (bn_end - bn_begin == 0) {
        lc1 = lvm.get_candidate_set(order[depth]);
        lc_count1 = lvm.get_candidate_set_size(order[depth]);
        lc2 = local_candidates_buffer1_[depth].data();

        lc_count2 = 0;

        for (uint32_t i = 0; i < lc_count1; ++i) {
            uint32_t v = lc1[i];

            if (!visited_[v]) {
                local_candidates_buffer2_[depth][lc_count2] = i;
                lc2[lc_count2++] = v;
            }
#ifdef ENABLE_PERFORMANCE_COUNTERS
            else {
                iso_conflict_count_ += 1;
            }
#endif
        }
        encoded_local_candidates_store_[depth] = local_candidates_buffer2_[depth].data();
        local_candidates_store_[depth] = lc2;
        return lc_count2;
    }
    else {
        uint32_t u = bn[bn_begin];
        uint32_t v = embedding_[u];
        uint32_t encoded_v = encoded_embedding_[u];

        uint32_t view_id = view_mapping[bn_begin];
        if (query_graph->getVertexDegree(order[depth]) == 1) {
            const uint32_t sub_back = bn[bn_begin];
            const uint32_t sub_cur = order[depth];
            assert(sub_back < new2old.size() && sub_cur < new2old.size());
            const int32_t full_back = new2old[sub_back];
            const int32_t full_cur = new2old[sub_cur];
            assert(full_back >= 0 && full_cur >= 0);
            view_id = gvm.get_view_id({static_cast<uint32_t>(full_back), static_cast<uint32_t>(full_cur)});
        }

        lc1 = query_graph->getVertexDegree(order[depth]) == 1 ? gvm.get_view(view_id)->get_neighbor(v, lc_count1)
                : lvm.get_view(view_id)->get_neighbors(encoded_v, lc_count1);

        lc2 = local_candidates_buffer1_[depth].data();

        if (lc_count1 == 0) {
            goto EXIT_MAPPED;
        }

        if (bn_begin + 1 < bn_end) {
            bn_begin += 1;
            u = bn[bn_begin];
            encoded_v = encoded_embedding_[u];
            view_id = view_mapping[bn_begin];

            lc2 = lvm.get_view(view_id)->get_neighbors(encoded_v,lc_count2);

            uint32_t temp_count;
            uint32_t *temp_buffer = local_candidates_buffer1_[depth].data();

            ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, lc2, lc_count2, temp_buffer, temp_count);

            if (temp_count == 0) {
                lc_count1 = 0;
                goto EXIT_MAPPED;
            }

            lc1 = temp_buffer;
            lc_count1 = temp_count;
            temp_buffer = local_candidates_buffer2_[depth].data();

            for (bn_begin += 1; bn_begin < bn_end; ++bn_begin) {
                u = bn[bn_begin];
                encoded_v = encoded_embedding_[u];
                view_id = view_mapping[bn_begin];

                lc2 =  lvm.get_view(view_id)->get_neighbors(encoded_v,lc_count2);
                ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, lc2, lc_count2, temp_buffer, temp_count);

                if (temp_count == 0) return 0;
                std::swap(temp_buffer, lc1);
                std::swap(temp_count, lc_count1);
            }

            lc2 = temp_buffer;
        }
    }

EXIT_MAPPED:
    if (lc_count1 == 0) {
#ifdef ENABLE_PERFORMANCE_COUNTERS
       si_empty_count_ += 1;
#endif
        return 0;
    }

    if (query_graph->getVertexDegree(order[depth]) == 1) {
        lc_count2 = 0;

        for (uint32_t i = 0; i < lc_count1; ++i) {
            uint32_t v = lc1[i];

            if (!visited_[v]) {
                lc2[lc_count2++] = v;
            }
#ifdef ENABLE_PERFORMANCE_COUNTERS
            else {
                iso_conflict_count_ += 1;
            }
#endif
        }
        encoded_local_candidates_store_[depth] = local_candidates_buffer2_[depth].data();
        local_candidates_store_[depth] = lc2;
        return lc_count2;
    }
    else {
        lc_count2 = 0;
        uint32_t* candidate_set = lvm.get_candidate_set(order[depth]);
        uint32_t* encoded_buffer = one_bn ? local_candidates_buffer2_[depth].data() : lc1;
        for (uint32_t i = 0; i < lc_count1; ++i) {
            uint32_t encoded_v = lc1[i];
            uint32_t v = candidate_set[encoded_v];

            if (!visited_[v]) {
                encoded_buffer[lc_count2] = encoded_v;
                lc2[lc_count2++] = v;
            }
#ifdef ENABLE_PERFORMANCE_COUNTERS
            else {
                iso_conflict_count_ += 1;
            }
#endif
        }

        encoded_local_candidates_store_[depth] = encoded_buffer;
        local_candidates_store_[depth] = lc2;
        return lc_count2;
    }
}

   
