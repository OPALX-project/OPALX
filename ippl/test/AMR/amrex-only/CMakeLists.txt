FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test AMR found in ${_relPath}")

set (CMAKE_CXX_FLAGS
    "-O3 ${CMAKE_CXX_FLAGS}"
    )

set (AMREX_FLAGS "${AMREX_Fortran_FLAGS} ${AMREX_Fortran_DEFINITIONS}")
    
set (CMAKE_CXX_FLAGS
    "${IPPL_CMAKE_CXX_FLAGS} -std=c++11 ${CMAKE_CXX_FLAGS} ${AMREX_FLAGS}"
)

add_definitions(${AMREX_DEFINES})

include_directories (
    ${CMAKE_SOURCE_DIR}/ippl/src
    ${IPPL_INCLUDE_DIR}
    ${AMREX_INCLUDE_DIR}
)
    
link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${IPPL_LIBRARY_DIR}
    ${AMREX_LIBRARY_DIR}
    ${HYPRE_LIBRARY_DIR}
)

SET (LIBS
    ${IPPL_LIBRARY}
    ippl
    ${MPI_CXX_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
    ${Fortran_LIBS}
    ${AMREX_LIBRARIES}
    ${AMREX_EXTRA_CXX_LINK_LINE}
    ${AMREX_EXTRA_Fortran_LINK_LINE}
    m
    z
)

if ( ${AMREX_DIM} EQUAL 3 )
    add_executable( testDeposition
                    testDeposition.cpp)
    
    add_executable( AssignMultiLevelDensity
                    AssignMultiLevelDensity.cpp)

#    add_executable( testGeometry
#                    testGeometry.cpp)
    
    add_executable( testDepositionFail
                    testDepositionFail.cpp
                    ../Tagging_nd.f90
                    ../AmrOpal_F.h)
    
    add_executable( testSolverFail
                    testSolverFail.cpp
                    ../Tagging_nd.f90
                    ../AmrOpal_F.h
                    ../Solver.cpp)
    
    # 4. March 2017,
    # http://stackoverflow.com/questions/13429656/cmake-how-to-copy-contents-of-a-directory-into-build-directory-after-make
    add_custom_command(TARGET testDeposition testDepositionFail POST_BUILD
                       COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/inputs $<TARGET_FILE_DIR:testDeposition>)

    target_link_libraries (testDeposition ${LIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
                        
    target_link_libraries (AssignMultiLevelDensity ${LIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

#    target_link_libraries (testGeometry ${LIBS}
#                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testDepositionFail ${LIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testSolverFail ${LIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

set (EXECS
    testDeposition
    AssignMultiLevelDensity
    testDepositionFail
    testSolverFail
)

install (TARGETS ${EXECS} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/libexec/amr/amrex-only")

endif ( ${AMREX_DIM} EQUAL 3)
