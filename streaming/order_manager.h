//
// Created by sunsx on 28/05/21.
//

#ifndef RAPIDMATCH_ORDER_MANAGER_H
#define RAPIDMATCH_ORDER_MANAGER_H

#include "graph/graph.h"
#include "streaming_type.h"
#include "query_bctree.h"

class OrderManager {
private:
    std::vector<std::vector<uint32_t>> automorphisms_;
    std::vector<OrdersPerEdge> orders_;
    spp::sparse_hash_map<uint64_t, OrdersPerEdge> special_orders_;
    spp::sparse_hash_map<uint64_t, OrdersPerEdge> special_iso_orders_;
    std::vector<OrdersPerEdge> block_reuse_orders_;
    spp::sparse_hash_map<Edge, uint32_t> edge_orders_mapping_;
    spp::sparse_hash_map<LabelTriple, std::vector<Edge>> label_edge_mapping_;
    spp::sparse_hash_map<LabelTriple, std::vector<uint32_t>> label_automorphism_mapping_;
    std::vector<std::vector<Edge>> automorphism_edges_;

private:
    void detect_automorphism_edges(const Graph *query_graph);
    void create_indexing_order_for_reduced_graph(const Graph *query_graph, Edge edge, OrdersPerEdge &order);
    void create_matching_order_for_reduced_graph(const Graph *query_graph, Edge edge, OrdersPerEdge &order);
    void create_block_reuse_matching_order_for_reduced_graph(const Graph *query_graph,
                                                             Edge edge,
                                                             OrdersPerEdge &order,
                                                             const QueryBCTreeDecomposition& query_bctree);
    void create_special_matching_order_for_reduced_graph(const Graph *query_graph,
                                                        Edge edge,
                                                        OrdersPerEdge &order,
                                                        const std::vector<std::vector<int32_t>> &sub_new2old,
                                                        uint32_t sub_id);
    void create_special_matching_order_for_reduced_graph_iso(const Graph *query_graph,
                                                             Edge edge,
                                                             OrdersPerEdge &order,
                                                             const std::vector<std::vector<std::vector<int32_t>>> &iso_instance_vertices,
                                                             uint32_t iso_group_id);
    void generate_matching_order_with_RI(const Graph* graph, std::vector<uint32_t>& matching_order);
public:
    OrderManager() {}
    ~OrderManager() { release(); }

    void initialize(const Graph *query_graph);
    void release();

    OrdersPerEdge* get_orders(Edge edge);
    OrdersPerEdge* get_orders_forced_prefix(Edge edge);
    std::vector<Edge>* get_mapped_edges(LabelTriple label_triple);
    std::vector<uint32_t>* get_mapped_automorphism(LabelTriple labelTriple);
    std::tuple<Edge, OrdersPerEdge*, uint32_t> get_automorphism_meta(uint32_t id);
    std::tuple<Edge, OrdersPerEdge*, uint32_t> get_block_reuse_automorphism_meta(uint32_t id);
    std::tuple<Edge, OrdersPerEdge*, uint32_t> get_automorphism_meta_special(uint32_t automorphism_id, uint32_t sub_id);
    std::tuple<Edge, OrdersPerEdge*, uint32_t> get_automorphism_meta_special2(uint32_t automorphism_id, uint32_t iso_group_id);
    void create_orders(const Graph *query_graph);
    void create_orders_no_automorphism(const Graph *query_graph);
    void create_block_reuse_orders(const Graph *query_graph, const QueryBCTreeDecomposition& query_bctree);
    const std::vector<OrdersPerEdge>& get_block_reuse_orders() const { return block_reuse_orders_; }
    const std::vector<std::vector<Edge>>& get_automorphism_edges() const { return automorphism_edges_; }
    void create_special_orders(const Graph *query_graph, const std::vector<std::vector<int32_t>> &sub_new2old);
    void create_special_orders2(const Graph *query_graph, const std::vector<std::vector<std::vector<int32_t>>> &iso_instance_vertices);
    void print_info();
};
#endif 
