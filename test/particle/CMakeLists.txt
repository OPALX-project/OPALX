FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test FFT found in ${_relPath}")

set (CMAKE_CXX_FLAGS
    "-std=c++0x ${CMAKE_CXX_FLAGS}"
    )


include_directories (
    ${CMAKE_SOURCE_DIR}/src                    
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src                 
)

add_executable (PIC3d PIC3d.cpp)
target_link_libraries (PIC3d ippl)
add_executable (p3m3d p3m3d.cpp)
target_link_libraries (p3m3d ippl)

