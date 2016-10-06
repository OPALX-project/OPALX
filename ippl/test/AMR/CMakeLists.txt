FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test AMR found in ${_relPath}")
MESSAGE (STATUS "MPI COMPILER ${MPI_UNDERLYING_COMPILER}")

set (CMAKE_CXX_FLAGS
      "-O3 ${CMAKE_CXX_FLAGS}"
    )

include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/ippl/src
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
)

SET (IPPL_LIBS ippl)

add_executable (testAmrPartBunch testAmrPartBunch.cpp)


target_link_libraries (testAmrPartBunch ${IPPL_LIBS})


