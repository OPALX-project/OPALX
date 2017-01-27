FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test FFT found in ${_relPath}")

include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/ippl/src
    ${CLASSIC_SOURCE_DIR}
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
)

set (CMAKE_CXX_FLAGS
        "${IPPL_CMAKE_CXX_FLAGS} -std=c++11 -DUSEH5FEDV2 -DPARALLEL_IO ${CMAKE_CXX_FLAGS}"
    )

SET (IPPL_LIBS ippl)
SET (COMPILE_FLAGS ${OPAL_CXX_FLAGS})

IF (ENABLE_DKS)
   
   include_directories (${DKS_INCLUDE_DIR})
   link_directories (${DKS_LIBRARY_DIR})
   #SET (IPPL_LIBS ${IPPL_LIBS} dks)
   SET (IPPL_LIBS ${IPPL_LIBS} ${DKS_LIBRARY_DIR}/libdks.a)

ENDIF (ENABLE_DKS)

set_source_files_properties(SOURCE PIC3d.cpp p3m3d.cpp
                            PROPERTY COMPILE_FLAGS
                            -fno-tree-vrp)

#add_executable (PIC3d PIC3d.cpp)
#target_link_libraries (PIC3d ${IPPL_LIBS})

add_executable (p3m3dHeating p3m3dHeating.cpp)
target_link_libraries (p3m3dHeating ${IPPL_LIBS} ${H5Hut_LIBRARY} ${HDF5_LIBRARY} z ${CMAKE_DL_LIBS})

#add_executable (p3m3d p3m3d.cpp)
#target_link_libraries (p3m3d ${IPPL_LIBS})

#add_executable (p3m3dleapfrog p3m3dleapfrog.cpp)
#target_link_libraries (p3m3dleapfrog ${IPPL_LIBS})

#add_executable (p3m3dTwoStream p3m3dTwoStream.cpp)
#target_link_libraries (p3m3dTwoStream ${IPPL_LIBS})

#add_executable (p3m3dTwoStreamParallel p3m3dTwoStreamParallel.cpp)
#target_link_libraries (p3m3dTwoStreamParallel ${IPPL_LIBS})

#add_executable (p3m3dRegressionTests p3m3dRegressionTests.cpp)
#target_link_libraries (p3m3dRegressionTests ${IPPL_LIBS})

#add_executable (p3m3dMicrobunching p3m3dMicrobunching.cpp)
#target_link_libraries (p3m3dMicrobunching ${IPPL_LIBS} ${H5Hut_LIBRARY} ${HDF5_LIBRARY} z ${CMAKE_DL_LIBS})

#add_executable (p3m3d_periodic p3m3d_periodic.cpp)
#target_link_libraries (p3m3d_periodic ${IPPL_LIBS})

#add_executable (pm pm.cpp)
#target_link_libraries (pm ${IPPL_LIBS})

#add_executable (pp pp.cpp)
#target_link_libraries (pp ${IPPL_LIBS})

#add_executable (pp_periodic pp_periodic.cpp)
#target_link_libraries (pp_periodic ${IPPL_LIBS})

#add_executable (pm_periodic pm_periodic.cpp)
#target_link_libraries (pm_periodic ${IPPL_LIBS})

#add_executable (domaintest mini_exmpl_domain.cpp)
#target_link_libraries (domaintest ${IPPL_LIBS})

#add_executable (interpolation_test interpolation_test.cpp)
#target_link_libraries (interpolation_test ${IPPL_LIBS})

#add_executable (test-scatter-1 test-scatter-1.cpp)
#target_link_libraries (test-scatter-1 ${IPPL_LIBS})



