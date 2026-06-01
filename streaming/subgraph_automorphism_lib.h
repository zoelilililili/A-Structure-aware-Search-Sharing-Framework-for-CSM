#pragma once

#include <cstdint>
#include <vector>

class Graph;

namespace subgraph_automorphism_lib {

struct SubAutomorphismResult {
    uint32_t sub_id = 0;
    std::vector<uint32_t> old_vertices_sorted;
    uint64_t moved_old_vertices_mask = 0;
    std::vector<std::vector<uint32_t>> automorphisms;
    std::vector<std::pair<uint32_t, uint32_t>> vertex_list;
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> edge_list;
};
void build_maximal_automorphic_induced_subgraphs_in_memory(const Graph& query_graph,
                                                          uint64_t max_autos,
                                                          std::vector<SubAutomorphismResult>& out);

} 

