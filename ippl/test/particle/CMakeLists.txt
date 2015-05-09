FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test FFT found in ${_relPath}")

include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/ippl/src
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
)

SET (IPPL_LIBS ippl)

IF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_OPENMP)
   
   include_directories (${DKS_INCLUDE_DIR})
   link_directories (${DKS_LIBRARY_DIR})
   #SET (IPPL_LIBS ${IPPL_LIBS} dks)
   SET (IPPL_LIBS ${IPPL_LIBS} ${DKS_LIBRARY_DIR}/libdks.a)

ENDIF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_OPENMP)

set_source_files_properties(SOURCE PIC3d.cpp p3m3d.cpp
                            PROPERTY COMPILE_FLAGS
                            -fno-tree-vrp)
add_executable (PIC3d PIC3d.cpp)
target_link_libraries (PIC3d ${IPPL_LIBS})

add_executable (p3m3d p3m3d.cpp)
target_link_libraries (p3m3d ${IPPL_LIBS})

