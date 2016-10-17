FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test AMR found in ${_relPath}")
MESSAGE (STATUS "MPI COMPILER ${MPI_UNDERLYING_COMPILER}")

set (CMAKE_CXX_FLAGS
      "-O3 ${CMAKE_CXX_FLAGS} -DBL_USE_FORTRAN_MPI"
    )

include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/ippl/src
    ${CCSE_INCLUDE_DIR}
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CCSE_LIBRARY_DIR}
)

SET (IPPL_LIBS ippl)


add_executable (testAmrPartBunch testAmrPartBunch.cpp)
add_executable (testSolver testSolver.cpp)

# Boxlib has a circular dependency between -lcboxlib and -lfboxlib --> resolve by: -lcboxlib -lfboxlib -lcboxlib
target_link_libraries (testAmrPartBunch ${IPPL_LIBS} -lmpi -lmpi_mpifh -lgfortran ${CCSE_LIBRARIES} -lcboxlib)
target_link_libraries (testSolver ${IPPL_LIBS} -lmpi -lmpi_mpifh -lgfortran ${CCSE_LIBRARIES} -lcboxlib)


