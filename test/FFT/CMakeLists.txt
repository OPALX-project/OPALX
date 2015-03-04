FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test FFT found in ${_relPath}")

#IF(${MPI_UNDERLYING_COMPILER} STREQUAL "g++")
  set (CMAKE_CXX_FLAGS
      "-std=c++0x ${CMAKE_CXX_FLAGS}"
      )
#ENDIF(${MPI_UNDERLYING_COMPILER} STREQUAL "g++")

include_directories (
    ${CMAKE_SOURCE_DIR}/src                    
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src                 
)

add_executable (TestFFT TestFFT.cpp)
add_executable (TestRC TestRC.cpp)
target_link_libraries (TestFFT ippl)
target_link_libraries (TestRC ippl)

