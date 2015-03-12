FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test FFT found in ${_relPath}")
MESSAGE (STATUS "MPI COMPILER ${MPI_UNDERLYING_COMPILER}")

set (CMAKE_CXX_FLAGS
      "-O3 ${CMAKE_CXX_FLAGS}"
    )

include_directories (
    ${CMAKE_SOURCE_DIR}/src
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
)

SET (IPPL_LIBS ippl)

IF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_MIC)
include_directories (${DKS_INCLUDE_DIR})
link_directories (${DKS_LIBRARY_DIR})
SET (IPPL_LIBS ${IPPL_LIBS} dks)
ENDIF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_MIC)


add_executable (TestFFT TestFFT.cpp)
add_executable (TestRC TestRC.cpp)

target_link_libraries (TestFFT ${IPPL_LIBS})
target_link_libraries (TestRC ${IPPL_LIBS})


