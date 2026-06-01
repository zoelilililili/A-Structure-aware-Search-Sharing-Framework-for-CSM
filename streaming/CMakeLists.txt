set(STREAMING_ENGINE_FILE_LIST
        global_view_manager.cpp
        local_view_manager.cpp
        order_manager.cpp
        query_bctree.cpp
        search_engine.cpp
        subgraph_automorphism_lib.cpp
        subgraph_isomorphism_lib.cpp
        streaming_engine.cpp
        )
add_executable(RapidFlow.out ${STREAMING_ENGINE_FILE_LIST})
target_link_libraries(RapidFlow.out graph utility_core)
if (WIN32)
  target_link_libraries(RapidFlow.out dbghelp)
endif()
option(BUILD_SUBGRAPH_AUTOMORPHISM_TOOL "Build SubgraphAutomorphismTool.out" OFF)
if (BUILD_SUBGRAPH_AUTOMORPHISM_TOOL)
  add_executable(SubgraphAutomorphismTool.out subgraph_automorphism_tool.cpp)
  target_link_libraries(SubgraphAutomorphismTool.out graph utility_core)
endif ()