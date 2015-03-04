FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test ippltoyfdtd found in ${_relPath}")

set (IPPLTOYFDTD_SRCS
    ToyFDTD-DataSink/DataSink.cpp
)

IF(NOT ${MPI_UNDERLYING_COMPILER} STREQUAL "icpc")
  set (CMAKE_CXX_FLAGS
      "-std=c++0x ${CMAKE_CXX_FLAGS}"
      )
ENDIF(NOT ${MPI_UNDERLYING_COMPILER} STREQUAL "icpc")

include_directories (
    ${CMAKE_SOURCE_DIR}/src
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}/ToyFDTD-DataSink
    ${CMAKE_SOURCE_DIR}/src
)

SET (IPPL_LIBS ippl)

IF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_MIC)
include_directories (${DKS_INCLUDE_DIR})
link_directories (${DKS_LIBRARY_DIR})
SET (IPPL_LIBS ${IPPL_LIBS} dks)
ENDIF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_MIC)

add_executable (ipplToyFDTD ${IPPLTOYFDTD_SRCS} ToyFDTD-DataSink/ipplToyFDTD2.cpp)
target_link_libraries (ipplToyFDTD ${IPPL_LIBS})

