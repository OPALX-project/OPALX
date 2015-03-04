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

SET (IPPL_LIBS ippl)

IF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_MIC)
include_directories (${DKS_INCLUDE_DIR})
link_directories (${DKS_LIBRARY_DIR})
SET (IPPL_LIBS ${IPPL_LIBS} dks)
ENDIF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_MIC)


add_executable (PIC3d PIC3d.cpp)
target_link_libraries (PIC3d ${IPPL_LIBS})

add_executable (p3m3d p3m3d.cpp)
target_link_libraries (p3m3d ${IPPL_LIBS})

