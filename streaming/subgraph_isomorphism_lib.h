#pragma once
#include <cstdint>
#include <vector>
#include <tuple>
#include <utility>
class Graph;
namespace subgraph_isomorphism_lib {
struct SubIsomorphismGroupResult {
    uint32_t group_id = 0;                                     
    std::vector<std::pair<uint32_t, uint32_t>> vertex_list;     
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> edge_list; 

    std::vector<std::vector<uint32_t>> instances_mapping;
};

void build_maximal_isomorphic_induced_subgraphs_in_memory(const Graph& query_graph,
                                                          std::vector<SubIsomorphismGroupResult>& out,
                                                          const std::vector<uint64_t>* phase1_moved_masks = nullptr);

}