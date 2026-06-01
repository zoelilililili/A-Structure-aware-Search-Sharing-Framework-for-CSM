#include <queue>
#include <computesetintersection.h>
#include "order_manager.h"
#include "graphoperations.h"
void OrderManager::initialize(const Graph *query_graph) {}
void OrderManager::release() {
    automorphisms_.clear();
    orders_.clear();
    special_orders_.clear();
    block_reuse_orders_.clear();
    edge_orders_mapping_.clear();
    label_edge_mapping_.clear();
    label_automorphism_mapping_.clear();
    automorphism_edges_.clear();
}

OrdersPerEdge *OrderManager::get_orders(Edge edge) {
    auto it = edge_orders_mapping_.find(edge);
    if (it != edge_orders_mapping_.end()) {
        return &orders_[it->second];
    }

    return nullptr;
}

std::vector<uint32_t> *OrderManager::get_mapped_automorphism(LabelTriple label_triple) {
    auto it = label_automorphism_mapping_.find(label_triple);
    if (it != label_automorphism_mapping_.end()) {
        return &it->second;
    }

    return nullptr;
}

std::vector<Edge> *OrderManager::get_mapped_edges(LabelTriple label_triple) {
    auto it = label_edge_mapping_.find(label_triple);
    if (it != label_edge_mapping_.end()) {
        return &it->second;
    }
    return nullptr;
}

std::tuple<Edge, OrdersPerEdge*, uint32_t> OrderManager::get_automorphism_meta(uint32_t id) {
    return std::make_tuple(automorphism_edges_[id][0], &orders_[id], (uint32_t)(automorphism_edges_[id].size()));
}

std::tuple<Edge, OrdersPerEdge*, uint32_t> OrderManager::get_block_reuse_automorphism_meta(uint32_t id) {
    OrdersPerEdge* order = &block_reuse_orders_[id];
    return std::make_tuple(automorphism_edges_[id][0], order,
                           static_cast<uint32_t>(automorphism_edges_[id].size()));
}

std::tuple<Edge, OrdersPerEdge*, uint32_t> OrderManager::get_automorphism_meta_special(uint32_t automorphism_id, uint32_t sub_id) {
    const uint64_t key = (static_cast<uint64_t>(sub_id) << 32) | static_cast<uint64_t>(automorphism_id);
    auto it = special_orders_.find(key);
    if (it != special_orders_.end()) {
        // Edge 代表边仍沿用 automorphism 桶的第一个代表边；automorphism_size 沿用该桶大小。
        return std::make_tuple(automorphism_edges_[automorphism_id][0], &it->second,
                               static_cast<uint32_t>(automorphism_edges_[automorphism_id].size()));
    }
    return get_automorphism_meta(automorphism_id);
}

void OrderManager::create_block_reuse_matching_order_for_reduced_graph(const Graph *query_graph,
                                                                        Edge edge,
                                                                        OrdersPerEdge &order,
                                                                        const QueryBCTreeDecomposition& query_bctree) {
    create_matching_order_for_reduced_graph(query_graph, edge, order);
    const auto baseline_order = order.matching_order_;
    std::vector<uint32_t> outside_order;
    std::vector<std::vector<uint32_t>> block_inner_orders(query_bctree.block_count);
    outside_order.reserve(baseline_order.size());
    auto get_owned_block_id = [&](uint32_t qv) -> int32_t {
        if (qv < query_bctree.blocks.size() && !query_bctree.blocks[qv].empty()) {
            return static_cast<int32_t>(query_bctree.blocks[qv][0]);
        }
        return -1;
    };
    const int32_t edge_u_block = get_owned_block_id(edge.first);
    const int32_t edge_v_block = get_owned_block_id(edge.second);
    std::vector<char> force_outside_block(query_bctree.block_count, 0);
    if (edge_u_block >= 0 && static_cast<uint32_t>(edge_u_block) < force_outside_block.size()) {
        force_outside_block[edge_u_block] = 1;
    }
    if (edge_v_block >= 0 && static_cast<uint32_t>(edge_v_block) < force_outside_block.size()) {
        force_outside_block[edge_v_block] = 1;
    }

    for (uint32_t qv : baseline_order) {
        bool in_phase3_block = (qv < query_bctree.blocks.size() && !query_bctree.blocks[qv].empty());
        bool is_cut = (qv < query_bctree.is_cut_vertex.size() && query_bctree.is_cut_vertex[qv]);
        if (!in_phase3_block || is_cut) {
            outside_order.push_back(qv);
            continue;
        }
        uint32_t block_id = query_bctree.blocks[qv][0];
        if (block_id < force_outside_block.size() && force_outside_block[block_id]) {
            outside_order.push_back(qv);
            continue;
        }
        if (block_id < block_inner_orders.size()) {
            block_inner_orders[block_id].push_back(qv);
        } else {
            outside_order.push_back(qv);
        }
    }
    std::vector<uint32_t> new_matching_order;
    new_matching_order.reserve(baseline_order.size());
    new_matching_order.insert(new_matching_order.end(), outside_order.begin(), outside_order.end());
    for (uint32_t block_id = 0; block_id < block_inner_orders.size(); ++block_id) {
        new_matching_order.insert(new_matching_order.end(),
                                  block_inner_orders[block_id].begin(),
                                  block_inner_orders[block_id].end());
    }
    auto& updated_matching_order = order.matching_order_;
    auto& updated_matching_order_bn = order.matching_order_bn_;
    auto& updated_matching_order_bn_offset = order.matching_order_bn_offset_;
    auto& updated_matching_order_view_mapping = order.matching_order_view_mappings_;
    auto& updated_matching_order_edge_type = order.matching_order_edge_type_;
    updated_matching_order = std::move(new_matching_order);
    updated_matching_order_bn.clear();
    updated_matching_order_bn_offset.clear();
    updated_matching_order_view_mapping.clear();
    updated_matching_order_edge_type.clear();
    updated_matching_order_bn_offset.push_back(0);
    for (uint32_t depth = 0; depth < updated_matching_order.size(); ++depth) {
        uint32_t u = updated_matching_order[depth];
        for (uint32_t prev = 0; prev < depth; ++prev) {
            uint32_t p = updated_matching_order[prev];
            if (query_graph->checkEdgeExistence(u, p)) {
                updated_matching_order_bn.push_back(p);
                updated_matching_order_view_mapping.push_back(0);
                updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
            }
        }
        if (depth > 0 && query_graph->getVertexDegree(u) == 1 && query_graph->checkEdgeExistence(u, edge.first)) {
            updated_matching_order_bn.push_back(edge.first);
            updated_matching_order_view_mapping.push_back(0);
            updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
        }
        if (depth > 0 && query_graph->getVertexDegree(u) == 1 && query_graph->checkEdgeExistence(u, edge.second)) {
            updated_matching_order_bn.push_back(edge.second);
            updated_matching_order_view_mapping.push_back(0);
            updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
        }
        updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
    }
}

void OrderManager::create_block_reuse_orders(const Graph *query_graph, const QueryBCTreeDecomposition& query_bctree) {
    if (automorphism_edges_.empty()) {
        detect_automorphism_edges(query_graph);
    }

    block_reuse_orders_.clear();
    for (auto& automorphism : automorphism_edges_) {
        Edge edge = automorphism.front();
        OrdersPerEdge order;
        create_indexing_order_for_reduced_graph(query_graph, edge, order);
        create_block_reuse_matching_order_for_reduced_graph(query_graph, edge, order, query_bctree);
        block_reuse_orders_.push_back(std::move(order));
    }
}

void OrderManager::detect_automorphism_edges(const Graph *query_graph) {
    GraphOperations::compute_automorphism(query_graph, automorphisms_);

    uint32_t n = query_graph->getVerticesCount();
    spp::sparse_hash_set<Edge> selected;

    for (uint32_t u = 0; u < n; ++u) {
        uint32_t u_nbr_count;
        auto u_nbr = query_graph->getVertexNeighbors(u, u_nbr_count);

        for (uint32_t i = 0; i < u_nbr_count; ++i) {
            uint32_t uu = u_nbr[i];
            Edge e = {u, uu};
            if (!selected.contains(e)) {
                selected.insert(e);
                automorphism_edges_.push_back({e});
                for (auto &embedding: automorphisms_) {
                    Edge mapped_e = {embedding[u], embedding[uu]};
                    if (!selected.contains(mapped_e)) {
                        selected.insert(mapped_e);
                        automorphism_edges_.back().push_back(mapped_e);
                    }
                }
            }
        }
    }
}

void OrderManager::create_indexing_order_for_reduced_graph(const Graph *query_graph, Edge edge, OrdersPerEdge &order) {
    uint32_t n = query_graph->getVerticesCount();
    auto& indexing_order = order.indexing_order_;
    auto& indexing_order_bn = order.indexing_order_bn_;
    auto& indexing_order_bn_offset = order.indexing_order_bn_offset_;
    indexing_order = {edge.first, edge.second};
    indexing_order_bn = {edge.first};
    indexing_order_bn_offset = {0, 0, 1};
    std::vector<bool> visited(n);
    visited[edge.first] = true;
    visited[edge.second] = true;
    std::vector<uint32_t> adjacent_to_both;
    std::vector<uint32_t> adjacent_to_one;
    std::vector<uint32_t> adjacent_to_none;
    for (uint32_t u = 0; u < n; ++u) {
        if (u == edge.first || u == edge.second)
            continue;
        bool f1 = query_graph->checkEdgeExistence(u, edge.first);
        bool f2 = query_graph->checkEdgeExistence(u, edge.second);
        if (f1 && f2) {
            adjacent_to_both.push_back(u);
        }
        else if (!f1 && !f2) {
            adjacent_to_none.push_back(u);
        }
        else {
            adjacent_to_one.push_back(u);
        }
    }

    auto generate_function = [query_graph, &visited, &indexing_order, &indexing_order_bn, &indexing_order_bn_offset]
            (std::vector<uint32_t>& target_vertex) {
        for (uint32_t i = 0; i < target_vertex.size(); ++i) {
            std::vector<uint32_t> current_vertex_bn;
            std::vector<uint32_t> selected_vertex_bn;
            uint32_t selected_vertex;

            for (auto u : target_vertex) {
                if (!visited[u]) {
                    current_vertex_bn.clear();
                    for (auto uu : indexing_order) {
                        if (query_graph->checkEdgeExistence(u, uu)) {
                            current_vertex_bn.push_back(uu);
                        }
                    }
                    if (current_vertex_bn.size() > selected_vertex_bn.size()) {
                        current_vertex_bn.swap(selected_vertex_bn);
                        selected_vertex = u;
                    }
                }
            }
            indexing_order.push_back(selected_vertex);
            indexing_order_bn.insert(indexing_order_bn.end(), selected_vertex_bn.begin(), selected_vertex_bn.end());
            indexing_order_bn_offset.push_back(indexing_order_bn.size());
            visited[selected_vertex] = true;
        }
    };
    order.triangle_end_index_ = adjacent_to_both.size() + 2;
    order.adjacent_end_index_ = adjacent_to_both.size() + adjacent_to_one.size() + 2;

    generate_function(adjacent_to_both);
    generate_function(adjacent_to_one);
    generate_function(adjacent_to_none);
}

void OrderManager::create_matching_order_for_reduced_graph(const Graph *query_graph, Edge edge, OrdersPerEdge &order) {
    /**
     * 1. Get all connected components.
     * 2. Generate the matching order for each connected component.
     * 3. Merge the connected component.
     */
     uint32_t n = query_graph->getVerticesCount();
     std::vector<bool> visited(n, false);
     visited[edge.first] = true;
     visited[edge.second] = true;

     std::vector<std::vector<uint32_t>> connected_components;
     std::queue<uint32_t> q;

     // 1. Generate connected components by conducting a BFS from each vertex.
     for (uint32_t u = 0; u < n; ++u) {
         if (!visited[u]) {
             std::vector<uint32_t> component;
             component.push_back(u);
             q.push(u);
             visited[u] = true;

             while (!q.empty()) {
                 uint32_t uu = q.front();
                 q.pop();
                 uint32_t uu_nbrs_count;
                 auto uu_nbrs = query_graph->getVertexNeighbors(uu, uu_nbrs_count);

                 for (uint32_t i = 0; i < uu_nbrs_count; ++i) {
                     uint32_t uuu = uu_nbrs[i];
                     if (!visited[uuu]) {
                         q.push(uuu);
                         component.push_back(uuu);
                         visited[uuu] = true;
                     }
                 }
             }

             connected_components.emplace_back(component);
         }
     }

     // 2. Generate the matching order for each connected component.
     std::vector<std::vector<uint32_t>> matching_orders;

     // map the old vertex id to new vertex id
     std::vector<uint32_t> reverse_mapping(n, 0);

     std::vector<Graph*> graphs;

     for (auto& component : connected_components) {
         if (component.size() == 1 || component.size() == 2) {
             matching_orders.push_back(component);
             graphs.push_back(nullptr);
         }
         else {
             for (int i = 0; i < component.size(); ++i) {
                reverse_mapping[component[i]] = i;
             }

             // Create vertex list (vertex_id, vertex_label_id)  and edge list (begin_vertex_id, end_vertex_id, edge_label_id)
             std::vector<std::pair<uint32_t, uint32_t>> vertex_list;
             std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> edge_list;

             for (auto u: component) {
                 uint32_t u_id = reverse_mapping[u];
                 uint32_t u_label = query_graph->getVertexLabel(u);
                 vertex_list.emplace_back(u_id, u_label);

                 uint32_t u_nbrs_count;
                 auto u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);

                 for (uint32_t i = 0; i < u_nbrs_count; ++i) {
                     uint32_t uu = u_nbrs[i];
                     if (uu != edge.first && uu != edge.second) {
                         uint32_t uu_id = reverse_mapping[uu];
     // [无向边核对] 组件内边只存 u_id<uu_id 一次，避免无向重复；注意 getEdgeLabelByVertex 参数此处为组件重编号，需与 Graph::loadGraphFromMemory 约定一致。
                         if (u_id < uu_id) {
                             edge_list.emplace_back(u_id, uu_id, query_graph->getEdgeLabelByVertex(u_id, uu_id));
                         }
                     }
                 }
             }

             Graph *graph = new Graph(false);
             graph->is_edge_labeled = true;
             graph->loadGraphFromMemory(vertex_list, edge_list);

             std::vector<uint32_t> matching_order;
             generate_matching_order_with_RI(graph, matching_order);
             matching_orders.push_back(matching_order);
             graphs.push_back(graph);
         }
     }

     // 3. Merge the matching orders into one based on 1) the 2-core size; and 2) the component size.
     auto& updated_matching_order = order.matching_order_;
     auto& updated_matching_order_bn = order.matching_order_bn_;
     auto& updated_matching_order_bn_offset = order.matching_order_bn_offset_;
     auto& updated_matching_order_view_mapping = order.matching_order_view_mappings_;
     auto& updated_matching_order_edge_type = order.matching_order_edge_type_;

     updated_matching_order_bn_offset = {0};

     std::vector<uint32_t> core_size;
     std::vector<uint32_t> graph_size;

     // Initialize the core size and graph size for each subgraph.
     for (uint32_t i = 0; i < graphs.size(); ++i) {
         Graph* graph = graphs[i];
         if (graph != nullptr) {
             uint32_t size = 0;
             std::vector<int> core(graph->getVerticesCount(), 0);
             GraphOperations::getKCore(graph, core.data());

             for (auto core_value: core) {
                 if (core_value >= 2) {
                     size += 1;
                 }
             }
             core_size.push_back(size);
             graph_size.push_back(graph->getVerticesCount());
         }
         else {
             core_size.push_back(0);
             graph_size.push_back(connected_components[i].size());
         }
     }

     std::fill(visited.begin(), visited.end(), false);

     for (uint32_t i = 0; i < graphs.size(); ++i) {
         // Pick the graph.
         uint32_t selected_core_value = 0;
         uint32_t selected_graph_size = 0;
         uint32_t selected_graph = 0;

         for (uint32_t j = 0; j < graphs.size(); ++j) {
             if (!visited[j]) {
                 if ((core_size[j] > selected_core_value)
                     || (core_size[j] == selected_core_value && graph_size[j] > selected_graph_size)) {
                     selected_graph = j;
                     selected_core_value = core_size[j];
                     selected_graph_size = graph_size[j];
                 }
             }
         }

         visited[selected_graph] = true;

         // Update matching order and backward neighbors.
         if (graph_size[selected_graph] == 1) {
             uint32_t u = connected_components[selected_graph][0];
             updated_matching_order.push_back(u);
             if (query_graph->getVertexDegree(u) == 1) {
                 if (query_graph->checkEdgeExistence(u, edge.first)) {
                     updated_matching_order_bn.push_back(edge.first);
                 } else {
                     updated_matching_order_bn.push_back(edge.second);
                 }
             }
             updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
             updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
         } else if (graph_size[selected_graph] == 2) {
             if (query_graph->getVertexDegree(connected_components[selected_graph][0]) < query_graph->getVertexDegree(connected_components[selected_graph][1])) {
                 std::swap(connected_components[selected_graph][0], connected_components[selected_graph][1]);
             }

             updated_matching_order.insert(updated_matching_order.end(), connected_components[selected_graph].begin(), connected_components[selected_graph].end());
             // Insert the first vertex.
             updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());

             // Insert the second vertex.
             updated_matching_order_bn.push_back(connected_components[selected_graph][0]);
             updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
             updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
         } else {
             auto graph = graphs[selected_graph];
             std::vector<bool> vertex_visited(matching_orders[selected_graph].size(), false);
             for (auto u: matching_orders[selected_graph]) {
                 uint32_t u_nbrs_count;
                 auto u_nbrs = graph->getVertexNeighbors(u, u_nbrs_count);

                 for (uint32_t j = 0; j < u_nbrs_count; ++j) {
                     uint32_t uu = u_nbrs[j];
                     if (vertex_visited[uu]) {
                         updated_matching_order_bn.push_back(connected_components[selected_graph][uu]);
                         updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
                     }
                 }

                 updated_matching_order.push_back(connected_components[selected_graph][u]);
                 updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
                 vertex_visited[u] = true;
             }
         }
     }

    updated_matching_order_view_mapping.resize(updated_matching_order_bn.size());

     // Release graph.
     for (auto graph : graphs) {
         delete graph;
     }
}

void OrderManager::generate_matching_order_with_RI(const Graph *graph, std::vector<uint32_t> &matching_order) {
    uint32_t n = graph->getVerticesCount();
    std::vector<bool> visited(n, false);
    uint32_t selected_vertex = 0;
    uint32_t selected_vertex_selectivity = graph->getVertexDegree(selected_vertex);
    for (ui u = 1; u < n; ++u) {
        uint32_t u_selectivity = graph->getVertexDegree(u);
        if (u_selectivity > selected_vertex_selectivity) {
            selected_vertex = u;
            selected_vertex_selectivity = u_selectivity;
        }
    }
    matching_order.push_back(selected_vertex);
    visited[selected_vertex] = true;
    std::vector<uint32_t> tie_vertices;
    std::vector<uint32_t> temp;
    for (uint32_t i = 1; i < n; ++i) {
        selected_vertex_selectivity = 0;
        for (uint32_t u = 0; u < n; ++u) {
            if (!visited[u]) {
                uint32_t u_selectivity = 0;
                for (auto uu : matching_order) {
                    if (graph->checkEdgeExistence(u, uu)) {
                        u_selectivity += 1;
                    }
                }
                if (u_selectivity > selected_vertex_selectivity) {
                    selected_vertex_selectivity = u_selectivity;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                } else if (u_selectivity == selected_vertex_selectivity) {
                    tie_vertices.push_back(u);
                }
            }
        }
        if (tie_vertices.size() != 1) {
            temp.swap(tie_vertices);
            tie_vertices.clear();
            uint32_t count = 0;
            std::vector<uint32_t> u_fn;
            for (auto u : temp) {
                uint32_t un_count;
                auto un = graph->getVertexNeighbors(u, un_count);
                for (uint32_t j = 0; j < un_count; ++j) {
                    if (!visited[un[j]]) {
                        u_fn.push_back(un[j]);
                    }
                }
                uint32_t cur_count = 0;
                for (auto uu : matching_order) {
                    uint32_t uun_count;
                    auto uun = graph->getVertexNeighbors(uu, uun_count);
                    uint32_t common_neighbor_count = 0;
                    ComputeSetIntersection::ComputeCandidates(uun, uun_count, u_fn.data(), (uint32_t)u_fn.size(), common_neighbor_count);
                    if (common_neighbor_count > 0) {
                        cur_count += 1;
                    }
                }
                u_fn.clear();
                if (cur_count > count) {
                    count = cur_count;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                }
                else if (cur_count == count){
                    tie_vertices.push_back(u);
                }
            }
        }
        if (tie_vertices.size() != 1) {
            temp.swap(tie_vertices);
            tie_vertices.clear();
            uint32_t count = 0;
            std::vector<uint32_t> u_fn;
            for (auto u : temp) {
                uint32_t un_count;
                auto un = graph->getVertexNeighbors(u, un_count);
                for (uint32_t j = 0; j < un_count; ++j) {
                    if (!visited[un[j]]) {
                        u_fn.push_back(un[j]);
                    }
                }
                uint32_t cur_count = 0;
                for (auto uu : u_fn) {
                    bool valid = true;
                    for (auto uuu : matching_order) {
                        if (graph->checkEdgeExistence(uu, uuu)) {
                            valid = false;
                            break;
                        }
                    }
                    if (valid) {
                        cur_count += 1;
                    }
                }

                u_fn.clear();
                if (cur_count > count) {
                    count = cur_count;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                }
                else if (cur_count == count){
                    tie_vertices.push_back(u);
                }
            }
        }
        matching_order.push_back(tie_vertices[0]);
        visited[tie_vertices[0]] = true;
        tie_vertices.clear();
        temp.clear();
    }
}

void OrderManager::create_orders(const Graph *query_graph) {
    detect_automorphism_edges(query_graph);
    uint32_t count = 0;
    for (auto& automorphism : automorphism_edges_) {
        Edge edge = automorphism.front();
        orders_.emplace_back(OrdersPerEdge());

        create_indexing_order_for_reduced_graph(query_graph, edge, orders_.back());
        create_matching_order_for_reduced_graph(query_graph, edge, orders_.back());

        LabelTriple label_triple = { query_graph->getVertexLabel(edge.first),
                                     query_graph->getEdgeLabelByVertex(edge.first, edge.second), query_graph->getVertexLabel(edge.second)};
        {
            auto it = label_edge_mapping_.find(label_triple);
            if (it == label_edge_mapping_.end()) {
                auto temp_it = label_edge_mapping_.emplace(label_triple, std::vector<Edge>());
                it = temp_it.first;
            }

            for (auto e: automorphism) {
                edge_orders_mapping_.insert({e, orders_.size() - 1});
                it->second.push_back(e);
            }
        }
        {
            auto it = label_automorphism_mapping_.find(label_triple);
            if (it == label_automorphism_mapping_.end()) {
                auto temp_it = label_automorphism_mapping_.emplace(label_triple, std::vector<uint32_t>());
                it = temp_it.first;
            }
            it->second.push_back(orders_.size() - 1);
            count += 1;
        }
    }
}

void OrderManager::create_orders_no_automorphism(const Graph *query_graph) {
    uint32_t n = query_graph->getVerticesCount();
    spp::sparse_hash_set<Edge> selected;

    for (uint32_t u = 0; u < n; ++u) {
        uint32_t u_nbr_count;
        auto u_nbr = query_graph->getVertexNeighbors(u, u_nbr_count);
        for (uint32_t i = 0; i < u_nbr_count; ++i) {
            uint32_t uu = u_nbr[i];
            Edge e = {u, uu};
            if (!selected.contains(e)) {
                selected.insert(e);
                automorphism_edges_.push_back({e});
            }
        }
    }

    for (auto& automorphism : automorphism_edges_) {
        Edge edge = automorphism.front();
        orders_.emplace_back(OrdersPerEdge());

        create_indexing_order_for_reduced_graph(query_graph, edge, orders_.back());
        create_matching_order_for_reduced_graph(query_graph, edge, orders_.back());

        LabelTriple label_triple = { query_graph->getVertexLabel(edge.first),
                                     query_graph->getEdgeLabelByVertex(edge.first, edge.second),
                                     query_graph->getVertexLabel(edge.second)};
        {
            auto it = label_edge_mapping_.find(label_triple);
            if (it == label_edge_mapping_.end()) {
                auto temp_it = label_edge_mapping_.emplace(label_triple, std::vector<Edge>());
                it = temp_it.first;
            }

            for (auto e: automorphism) {
                edge_orders_mapping_.insert({e, orders_.size() - 1});
                it->second.push_back(e);
            }
        }
        {
            auto it = label_automorphism_mapping_.find(label_triple);
            if (it == label_automorphism_mapping_.end()) {
                auto temp_it = label_automorphism_mapping_.emplace(label_triple, std::vector<uint32_t>());
                it = temp_it.first;
            }
            it->second.push_back(orders_.size() - 1);
        }
    }
}

void OrderManager::create_special_orders(const Graph *query_graph, const std::vector<std::vector<int32_t>> &sub_new2old) {
    if (automorphism_edges_.empty()) {
        detect_automorphism_edges(query_graph);
    }
    special_orders_.clear();
    if (!sub_new2old.empty()) {
        special_orders_.reserve(static_cast<size_t>(automorphism_edges_.size()) * sub_new2old.size());
    }
    const uint32_t sub_count = static_cast<uint32_t>(sub_new2old.size());
    for (uint32_t automorphism_id = 0; automorphism_id < automorphism_edges_.size(); ++automorphism_id) {
        const Edge edge = automorphism_edges_[automorphism_id].front();
        for (uint32_t sub_id = 0; sub_id < sub_count; ++sub_id) {
            OrdersPerEdge special;
            create_indexing_order_for_reduced_graph(query_graph, edge, special);
            create_special_matching_order_for_reduced_graph(query_graph, edge, special, sub_new2old, sub_id);

            const uint64_t key = (static_cast<uint64_t>(sub_id) << 32) | static_cast<uint64_t>(automorphism_id);
            special_orders_.emplace(key, std::move(special));
        }
    }
}

void OrderManager::create_special_matching_order_for_reduced_graph(const Graph *query_graph,
                                                                   Edge edge,
                                                                   OrdersPerEdge &order,
                                                                   const std::vector<std::vector<int32_t>> &sub_new2old,
                                                                   uint32_t sub_id) {
     uint32_t n = query_graph->getVerticesCount();
     std::vector<bool> visited(n, false);
     std::vector<char> is_forced(n, 0);
     std::vector<uint32_t> forced_prefix;
     if (sub_id < sub_new2old.size()) {
         forced_prefix.reserve(sub_new2old[sub_id].size());
         for (int32_t old_u : sub_new2old[sub_id]) {
             if (old_u >= 0 && static_cast<uint32_t>(old_u) < n) {
                 const uint32_t u = static_cast<uint32_t>(old_u);
                 visited[u] = true;
                 is_forced[u] = 1;
                 // 与原 reduced-query 语义保持一致：锚边两端点由 data_edge 锚定，不应出现在 matching_order_（也不计入 forced 前缀长度）。
                 if (u == edge.first || u == edge.second) continue;
                 forced_prefix.push_back(u);
             }
         }
     }
     std::vector<std::vector<uint32_t>> connected_components;
     std::queue<uint32_t> q;
     for (uint32_t u = 0; u < n; ++u) {
         if (!visited[u]) {
             std::vector<uint32_t> component;
             component.push_back(u);
             q.push(u);
             visited[u] = true;
             while (!q.empty()) {
                 uint32_t uu = q.front();
                 q.pop();
                 uint32_t uu_nbrs_count;
                 auto uu_nbrs = query_graph->getVertexNeighbors(uu, uu_nbrs_count);
                 for (uint32_t i = 0; i < uu_nbrs_count; ++i) {
                     uint32_t uuu = uu_nbrs[i];
                     if (!visited[uuu]) {
                         q.push(uuu);
                         component.push_back(uuu);
                         visited[uuu] = true;
                     }
                 }
             }
             connected_components.emplace_back(component);
         }
     }
     std::vector<std::vector<uint32_t>> matching_orders;
     std::vector<uint32_t> reverse_mapping(n, 0);
     std::vector<Graph*> graphs;
     for (auto& component : connected_components) {
         if (component.size() == 1 || component.size() == 2) {
             matching_orders.push_back(component);
             graphs.push_back(nullptr);
         }
         else {
             for (int i = 0; i < component.size(); ++i) {
                reverse_mapping[component[i]] = i;
             }
             std::vector<std::pair<uint32_t, uint32_t>> vertex_list;
             std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> edge_list;
             for (auto u: component) {
                 uint32_t u_id = reverse_mapping[u];
                 uint32_t u_label = query_graph->getVertexLabel(u);
                 vertex_list.emplace_back(u_id, u_label);
                 uint32_t u_nbrs_count;
                 auto u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
                 for (uint32_t i = 0; i < u_nbrs_count; ++i) {
                     uint32_t uu = u_nbrs[i];
                     if (uu != edge.first && uu != edge.second) {
                         uint32_t uu_id = reverse_mapping[uu];
                         if (u_id < uu_id) {
                             edge_list.emplace_back(u_id, uu_id, query_graph->getEdgeLabelByVertex(u_id, uu_id));
                         }
                     }
                 }
             }
             Graph *graph = new Graph(false);
             graph->is_edge_labeled = true;
             graph->loadGraphFromMemory(vertex_list, edge_list);
             std::vector<uint32_t> matching_order;
             generate_matching_order_with_RI(graph, matching_order);
             matching_orders.push_back(matching_order);
             graphs.push_back(graph);
         }
     }
     auto& updated_matching_order = order.matching_order_;
     auto& updated_matching_order_bn = order.matching_order_bn_;
     auto& updated_matching_order_bn_offset = order.matching_order_bn_offset_;
     auto& updated_matching_order_view_mapping = order.matching_order_view_mappings_;
     auto& updated_matching_order_edge_type = order.matching_order_edge_type_;
     updated_matching_order.clear();
     updated_matching_order_bn.clear();
     updated_matching_order_bn_offset.clear();
     updated_matching_order_view_mapping.clear();
     updated_matching_order_edge_type.clear();
     updated_matching_order_bn_offset = {0};
     for (uint32_t fu : forced_prefix) {
         updated_matching_order.push_back(fu);
         updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
     }
     auto append_adjacent_forced_prefix_bn = [&](uint32_t old_u) {
         uint32_t nbr_count = 0;
         auto nbrs = query_graph->getVertexNeighbors(old_u, nbr_count);
         for (uint32_t k = 0; k < nbr_count; ++k) {
             const uint32_t uu = nbrs[k];
             if (uu == edge.first || uu == edge.second) continue;
             if (uu < is_forced.size() && is_forced[uu]) {
                 updated_matching_order_bn.push_back(uu);
                 updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
             }
         }
     };
     std::vector<uint32_t> core_size;
     std::vector<uint32_t> graph_size;
     for (uint32_t i = 0; i < graphs.size(); ++i) {
         Graph* graph = graphs[i];
         if (graph != nullptr) {
             uint32_t size = 0;
             std::vector<int> core(graph->getVerticesCount(), 0);
             GraphOperations::getKCore(graph, core.data());

             for (auto core_value: core) {
                 if (core_value >= 2) {
                     size += 1;
                 }
             }
             core_size.push_back(size);
             graph_size.push_back(graph->getVerticesCount());
         }
         else {
             core_size.push_back(0);
             graph_size.push_back(connected_components[i].size());
         }
     }

     std::fill(visited.begin(), visited.end(), false);

     for (uint32_t i = 0; i < graphs.size(); ++i) {
         // Pick the graph.
         uint32_t selected_core_value = 0;
         uint32_t selected_graph_size = 0;
         uint32_t selected_graph = 0;

         for (uint32_t j = 0; j < graphs.size(); ++j) {
             if (!visited[j]) {
                 if ((core_size[j] > selected_core_value)
                     || (core_size[j] == selected_core_value && graph_size[j] > selected_graph_size)) {
                     selected_graph = j;
                     selected_core_value = core_size[j];
                     selected_graph_size = graph_size[j];
                 }
             }
         }

         visited[selected_graph] = true;

         if (graph_size[selected_graph] == 1) {
             uint32_t u = connected_components[selected_graph][0];
             updated_matching_order.push_back(u);
             append_adjacent_forced_prefix_bn(u);
             if (query_graph->getVertexDegree(u) == 1) {
                 if (query_graph->checkEdgeExistence(u, edge.first)) {
                     if (!is_forced[edge.first]) {
                         updated_matching_order_bn.push_back(edge.first);
                         updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
                     }
                 } else {
                     if (!is_forced[edge.second]) {
                         updated_matching_order_bn.push_back(edge.second);
                         updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
                     }
                 }
             }
             updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
         } else if (graph_size[selected_graph] == 2) {
             if (query_graph->getVertexDegree(connected_components[selected_graph][0]) < query_graph->getVertexDegree(connected_components[selected_graph][1])) {
                 std::swap(connected_components[selected_graph][0], connected_components[selected_graph][1]);
             }

             updated_matching_order.push_back(connected_components[selected_graph][0]);
             append_adjacent_forced_prefix_bn(connected_components[selected_graph][0]);
             updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());

             updated_matching_order.push_back(connected_components[selected_graph][1]);
             append_adjacent_forced_prefix_bn(connected_components[selected_graph][1]);
             updated_matching_order_bn.push_back(connected_components[selected_graph][0]);
             updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
             updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
         } else {
             auto graph = graphs[selected_graph];
             std::vector<bool> vertex_visited(matching_orders[selected_graph].size(), false);
             for (auto u: matching_orders[selected_graph]) {
                 const uint32_t old_u = connected_components[selected_graph][u];
                 append_adjacent_forced_prefix_bn(old_u);
                 uint32_t u_nbrs_count;
                 auto u_nbrs = graph->getVertexNeighbors(u, u_nbrs_count);
                 for (uint32_t j = 0; j < u_nbrs_count; ++j) {
                     uint32_t uu = u_nbrs[j];
                     if (vertex_visited[uu]) {
                         updated_matching_order_bn.push_back(connected_components[selected_graph][uu]);
                         updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
                     }
                 }
                 updated_matching_order.push_back(old_u);
                 updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
                 vertex_visited[u] = true;
             }
         }
     }

    updated_matching_order_view_mapping.resize(updated_matching_order_bn.size());

     for (auto graph : graphs) {
         delete graph;
     }
}


void OrderManager::create_special_orders2(const Graph *query_graph, const std::vector<std::vector<std::vector<int32_t>>> &iso_instance_vertices) {
    if (automorphism_edges_.empty()) {
        detect_automorphism_edges(query_graph);
    }
    special_iso_orders_.clear();
    if (!iso_instance_vertices.empty()) {
        special_iso_orders_.reserve(static_cast<size_t>(automorphism_edges_.size()) * iso_instance_vertices.size());
    }
    const uint32_t iso_count = static_cast<uint32_t>(iso_instance_vertices.size());
    for (uint32_t automorphism_id = 0; automorphism_id < automorphism_edges_.size(); ++automorphism_id) {
        const Edge edge = automorphism_edges_[automorphism_id].front();
        for (uint32_t iso_group_id = 0; iso_group_id < iso_count; ++iso_group_id) {
            OrdersPerEdge special;
            create_indexing_order_for_reduced_graph(query_graph, edge, special);
            create_special_matching_order_for_reduced_graph_iso(query_graph, edge, special, iso_instance_vertices, iso_group_id);

            const uint64_t key = (static_cast<uint64_t>(iso_group_id) << 32) | static_cast<uint64_t>(automorphism_id);
            special_iso_orders_.emplace(key, std::move(special));
        }
    }
}

void OrderManager::create_special_matching_order_for_reduced_graph_iso(const Graph *query_graph,
                                                                        Edge edge,
                                                                        OrdersPerEdge &order,
                                                                        const std::vector<std::vector<std::vector<int32_t>>> &iso_instance_vertices,
                                                                        uint32_t iso_group_id) {
    uint32_t n = query_graph->getVerticesCount();
    std::vector<bool> visited(n, false);
    std::vector<char> is_forced(n, 0);
    std::vector<uint32_t> forced_prefix;

    if (iso_group_id < iso_instance_vertices.size()) {
        const auto& instances = iso_instance_vertices[iso_group_id];
        int32_t chosen_instance = -1;
        for (int32_t inst = 0; inst < static_cast<int32_t>(instances.size()); ++inst) {
            bool has_u = false, has_v = false;
            for (int32_t vid : instances[inst]) {
                if (static_cast<uint32_t>(vid) == edge.first) has_u = true;
                if (static_cast<uint32_t>(vid) == edge.second) has_v = true;
            }
            if (has_u && has_v) { chosen_instance = inst; break; }
        }
        if (chosen_instance < 0 && !instances.empty()) {
            chosen_instance = 0;
        }
        if (chosen_instance >= 0) {
            const auto& old_vertices = instances[chosen_instance];
            forced_prefix.reserve(old_vertices.size());
            for (int32_t old_u : old_vertices) {
                if (old_u >= 0 && static_cast<uint32_t>(old_u) < n) {
                    const uint32_t u = static_cast<uint32_t>(old_u);
                    visited[u] = true;
                    is_forced[u] = 1;
                    if (u == edge.first || u == edge.second) continue;
                    forced_prefix.push_back(u);
                }
            }
        }
    }
    std::vector<std::vector<uint32_t>> connected_components;
    std::queue<uint32_t> q;
    for (uint32_t u = 0; u < n; ++u) {
        if (!visited[u]) {
            std::vector<uint32_t> component;
            component.push_back(u);
            q.push(u);
            visited[u] = true;
            while (!q.empty()) {
                uint32_t uu = q.front();
                q.pop();
                uint32_t uu_nbrs_count;
                auto uu_nbrs = query_graph->getVertexNeighbors(uu, uu_nbrs_count);
                for (uint32_t i = 0; i < uu_nbrs_count; ++i) {
                    uint32_t uuu = uu_nbrs[i];
                    if (!visited[uuu]) {
                        q.push(uuu);
                        component.push_back(uuu);
                        visited[uuu] = true;
                    }
                }
            }
            connected_components.emplace_back(component);
        }
    }
    std::vector<std::vector<uint32_t>> matching_orders;
    std::vector<uint32_t> reverse_mapping(n, 0);
    std::vector<Graph*> graphs;
    for (auto& component : connected_components) {
        if (component.size() == 1 || component.size() == 2) {
            matching_orders.push_back(component);
            graphs.push_back(nullptr);
        }
        else {
            for (int i = 0; i < component.size(); ++i) {
               reverse_mapping[component[i]] = i;
            }
            std::vector<std::pair<uint32_t, uint32_t>> vertex_list;
            std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> edge_list;
            for (auto u: component) {
                uint32_t u_id = reverse_mapping[u];
                uint32_t u_label = query_graph->getVertexLabel(u);
                vertex_list.emplace_back(u_id, u_label);
                uint32_t u_nbrs_count;
                auto u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
                for (uint32_t i = 0; i < u_nbrs_count; ++i) {
                    uint32_t uu = u_nbrs[i];
                    if (uu != edge.first && uu != edge.second) {
                        uint32_t uu_id = reverse_mapping[uu];
                        if (u_id < uu_id) {
                            edge_list.emplace_back(u_id, uu_id, query_graph->getEdgeLabelByVertex(u_id, uu_id));
                        }
                    }
                }
            }
            Graph *graph = new Graph(false);
            graph->is_edge_labeled = true;
            graph->loadGraphFromMemory(vertex_list, edge_list);
            std::vector<uint32_t> matching_order;
            generate_matching_order_with_RI(graph, matching_order);
            matching_orders.push_back(matching_order);
            graphs.push_back(graph);
        }
    }
    auto& updated_matching_order = order.matching_order_;
    auto& updated_matching_order_bn = order.matching_order_bn_;
    auto& updated_matching_order_bn_offset = order.matching_order_bn_offset_;
    auto& updated_matching_order_view_mapping = order.matching_order_view_mappings_;
    auto& updated_matching_order_edge_type = order.matching_order_edge_type_;
    updated_matching_order.clear();
    updated_matching_order_bn.clear();
    updated_matching_order_bn_offset.clear();
    updated_matching_order_view_mapping.clear();
    updated_matching_order_edge_type.clear();
    updated_matching_order_bn_offset = {0};
    for (uint32_t fu : forced_prefix) {
        updated_matching_order.push_back(fu);
        updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
    }
    auto append_adjacent_forced_prefix_bn = [&](uint32_t old_u) {
        uint32_t nbr_count = 0;
        auto nbrs = query_graph->getVertexNeighbors(old_u, nbr_count);
        for (uint32_t k = 0; k < nbr_count; ++k) {
            const uint32_t uu = nbrs[k];
            if (uu == edge.first || uu == edge.second) continue;
            if (uu < is_forced.size() && is_forced[uu]) {
                updated_matching_order_bn.push_back(uu);
                updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
            }
        }
    };
    std::vector<uint32_t> core_size;
    std::vector<uint32_t> graph_size;
    for (uint32_t i = 0; i < graphs.size(); ++i) {
        Graph* graph = graphs[i];
        if (graph != nullptr) {
            uint32_t size = 0;
            std::vector<int> core(graph->getVerticesCount(), 0);
            GraphOperations::getKCore(graph, core.data());

            for (auto core_value: core) {
                if (core_value >= 2) {
                    size += 1;
                }
            }
            core_size.push_back(size);
            graph_size.push_back(graph->getVerticesCount());
        }
        else {
            core_size.push_back(0);
            graph_size.push_back(connected_components[i].size());
        }
    }

    std::fill(visited.begin(), visited.end(), false);

    for (uint32_t i = 0; i < graphs.size(); ++i) {
        uint32_t selected_core_value = 0;
        uint32_t selected_graph_size = 0;
        uint32_t selected_graph = 0;

        for (uint32_t j = 0; j < graphs.size(); ++j) {
            if (!visited[j]) {
                if ((core_size[j] > selected_core_value)
                    || (core_size[j] == selected_core_value && graph_size[j] > selected_graph_size)) {
                    selected_graph = j;
                    selected_core_value = core_size[j];
                    selected_graph_size = graph_size[j];
                }
            }
        }

        visited[selected_graph] = true;

        if (graph_size[selected_graph] == 1) {
            uint32_t u = connected_components[selected_graph][0];
            updated_matching_order.push_back(u);
            append_adjacent_forced_prefix_bn(u);
            if (query_graph->getVertexDegree(u) == 1) {
                if (query_graph->checkEdgeExistence(u, edge.first)) {
                    if (!is_forced[edge.first]) {
                        updated_matching_order_bn.push_back(edge.first);
                        updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
                    }
                } else {
                    if (!is_forced[edge.second]) {
                        updated_matching_order_bn.push_back(edge.second);
                        updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
                    }
                }
            }
            updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
        } else if (graph_size[selected_graph] == 2) {
            if (query_graph->getVertexDegree(connected_components[selected_graph][0]) < query_graph->getVertexDegree(connected_components[selected_graph][1])) {
                std::swap(connected_components[selected_graph][0], connected_components[selected_graph][1]);
            }

            updated_matching_order.push_back(connected_components[selected_graph][0]);
            append_adjacent_forced_prefix_bn(connected_components[selected_graph][0]);
            updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());

            updated_matching_order.push_back(connected_components[selected_graph][1]);
            append_adjacent_forced_prefix_bn(connected_components[selected_graph][1]);
            updated_matching_order_bn.push_back(connected_components[selected_graph][0]);
            updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
            updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
        } else {
            auto graph = graphs[selected_graph];
            std::vector<bool> vertex_visited(matching_orders[selected_graph].size(), false);
            for (auto u: matching_orders[selected_graph]) {
                const uint32_t old_u = connected_components[selected_graph][u];
                append_adjacent_forced_prefix_bn(old_u);
                uint32_t u_nbrs_count;
                auto u_nbrs = graph->getVertexNeighbors(u, u_nbrs_count);
                for (uint32_t j = 0; j < u_nbrs_count; ++j) {
                    uint32_t uu = u_nbrs[j];
                    if (vertex_visited[uu]) {
                        updated_matching_order_bn.push_back(connected_components[selected_graph][uu]);
                        updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
                    }
                }
                updated_matching_order.push_back(old_u);
                updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
                vertex_visited[u] = true;
            }
        }
    }

    updated_matching_order_view_mapping.resize(updated_matching_order_bn.size());

    for (auto graph : graphs) {
        delete graph;
    }
}

std::tuple<Edge, OrdersPerEdge*, uint32_t> OrderManager::get_automorphism_meta_special2(uint32_t automorphism_id, uint32_t iso_group_id) {
    const uint64_t key = (static_cast<uint64_t>(iso_group_id) << 32) | static_cast<uint64_t>(automorphism_id);
    auto it = special_iso_orders_.find(key);
    if (it != special_iso_orders_.end()) {
        return std::make_tuple(automorphism_edges_[automorphism_id][0], &it->second,
                               static_cast<uint32_t>(automorphism_edges_[automorphism_id].size()));
    }
    return get_automorphism_meta(automorphism_id);
}


void OrderManager::print_info() {
    printf("Order Info:\n");
    printf("The number of automorphisms: %zu\n", automorphisms_.size());

    for (auto& automorphism : automorphisms_) {
        for (auto u : automorphism) {
            printf("%u ", u);
        }
        printf("\n");
    }

    printf("-----\n");
    printf("The number of disjoint set: %zu\n", automorphism_edges_.size());
    for (uint32_t i = 0; i < automorphism_edges_.size(); ++i) {
        printf("Edge set: ");
        for (auto e : automorphism_edges_[i]) {
            printf("(%u, %u), ", e.first, e.second);
        }
        printf("\n");
        printf("Index order: ");
        for (auto u : orders_[i].indexing_order_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Index order bn offset: ");
        for (auto u : orders_[i].indexing_order_bn_offset_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Index order bn: ");
        for (auto u : orders_[i].indexing_order_bn_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Matching order: ");
        for (auto u : orders_[i].matching_order_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Matching order bn offset: ");
        for (auto u : orders_[i].matching_order_bn_offset_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Matching order bn: ");
        for (auto u : orders_[i].matching_order_bn_) {
            printf("%u ", u);
        }
        printf("\n");
        printf("-----\n");
    }
}
