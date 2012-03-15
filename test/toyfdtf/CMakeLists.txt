message("adding test ippltoyfdtd found in ${CMAKE_CURRENT_SOURCE_DIR}")

set (IPPLTOYFDTD_SRCS
    ToyFDTD-DataSink/DataSink.cpp
)

set (OLD_CMAKE_CXX_FLAGS
    ${CMAKE_CXX_FLAGS}
)

set (CMAKE_CXX_FLAGS
    "-std=c++0x ${CMAKE_CXX_FLAGS}"
)

include_directories (
    ${CMAKE_SOURCE_DIR}/src                    
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}/ToyFDTD-DataSink
    ${CMAKE_SOURCE_DIR}/src                 
)

add_executable (ipplToyFDTD ${IPPLTOYFDTD_SRCS} ToyFDTD-DataSink/ipplToyFDTD2.cpp)
target_link_libraries (ipplToyFDTD ippl)

#set (CMAKE_CXX_FLAGS
#    ${OLD_CMAKE_CXX_FLAGS}
#)
