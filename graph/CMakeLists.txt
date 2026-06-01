set(GRAPH_SRC
        graph.cpp)
add_library(graph
        ${GRAPH_SRC})
target_link_libraries(graph utility_core)