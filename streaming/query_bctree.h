#ifndef RAPIDFLOW_QUERY_BCTREE_H
#define RAPIDFLOW_QUERY_BCTREE_H

#include <cstdint>
#include <iostream>
#include <vector>
class Graph;
struct QueryBCTreeDecomposition {
    uint32_t query_vertex_count = 0;
    uint32_t block_count = 0;
    std::vector<std::vector<uint32_t>> blocks;
    std::vector<std::vector<uint32_t>> block_neighbor_block_indices;
    std::vector<bool> is_cut_vertex;
    std::vector<uint32_t> vertex_block_membership_count;
    std::vector<char> is_inner_non_cut;

    void clear();
};

QueryBCTreeDecomposition build_query_bctree_decomposition(const Graph* query_graph);

void print_query_bctree_decomposition(const QueryBCTreeDecomposition& d, std::ostream& os);

#endif
