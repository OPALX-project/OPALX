set (CLASSIC_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
add_subdirectory(Algorithms)
add_subdirectory(MemoryManagement)
add_subdirectory(Solvers)

set(OPAL_SRCS "${OPAL_SRCS}" PARENT_SCOPE)
