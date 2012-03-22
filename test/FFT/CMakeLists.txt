message("adding test ippltoyfdtd found in ${CMAKE_CURRENT_SOURCE_DIR}")

set (OLD_CMAKE_CXX_FLAGS
    ${CMAKE_CXX_FLAGS}
)


IF(${MPI_UNDERLYING_COMPILER} STREQUAL "g++")
  set (CMAKE_CXX_FLAGS
      "-std=c++0x ${CMAKE_CXX_FLAGS}"
      )
ENDIF(${MPI_UNDERLYING_COMPILER} STREQUAL "g++")

include_directories (
    ${CMAKE_SOURCE_DIR}/src                    
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src                 
)

add_executable (TestFFT TestFFT.cpp)
target_link_libraries (TestFFT ippl)

#set (CMAKE_CXX_FLAGS
#    ${OLD_CMAKE_CXX_FLAGS}
#)
