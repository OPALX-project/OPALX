FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test AMR found in ${_relPath}")
#     MESSAGE (STATUS "MPI COMPILER ${MPI_UNDERLYING_COMPILER}")
    
set (CMAKE_CXX_FLAGS
    "-O3 ${CMAKE_CXX_FLAGS} -DBL_USE_FORTRAN_MPI -DUNIQUE_PTR -DHAVE_AMR_MG_SOLVER"
    )
    
set (CMAKE_CXX_FLAGS
    "${IPPL_CMAKE_CXX_FLAGS} -std=c++11 -DUSEH5FEDV2 -DPARALLEL_IO -DIPPL_AMR ${CMAKE_CXX_FLAGS}"
)

set (CMAKE_CXX_FLAGS "-DHAVE_AMR_MG_SOLVER ${CMAKE_CXX_FLAGS}")

add_definitions(${AMREX_DEFINES})

add_compile_options(-Wno-deprecated-declarations)

include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/ippl/src
    ${CMAKE_SOURCE_DIR}/src/Classic/
    ${CMAKE_SOURCE_DIR}/src/Classic/Utilities/
    ${CMAKE_SOURCE_DIR}/src/Utilities/
    ${IPPL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${AMREX_INCLUDE_DIR}
    ${HYPRE_INCLUDE_DIR}
    ${BOOST_INCLUDE_DIR}
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
)
    
link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../ippl-amr/
    ${CMAKE_SOURCE_DIR}/src
    ${IPPL_LIBRARY_DIR}
    ${AMREX_LIBRARY_DIR}
    ${HYPRE_LIBRARY_DIR}
    ${BOOST_LIBRARY_DIR}
    ${Trilinos_LIBRARY_DIRS}
    ${Trilinos_TPL_LIBRARY_DIRS}
)
    
SET (LIBS
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARIES}
    ${IPPL_LIBRARY}
    ippl
    ${Boost_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
    ${Fortran_LIBS}
    ${AMREX_LIBRARIES}
    m
    z
)

set(F90_source_files ../Tagging_nd.f90)

set (TESTLIBS amropal multigrid)

set(EXECS "")

IF ( ${AMREX_DIM} EQUAL 2 )

    add_executable( testRestrictionMatrix testRestrictionMatrix.cpp)
    
    add_executable( testInterpolationMatrix testInterpolationMatrix.cpp)
    
    add_executable( testBoundaryMatrix testBoundaryMatrix.cpp)
    
    add_executable( testSmootherMatrix testSmootherMatrix.cpp)
    
    add_executable( testPoissonMatrix testPoissonMatrix.cpp)
    
    add_executable( testNoFine testNoFine.cpp)
    
    add_executable( testWithFine testWithFine.cpp)
    
    add_executable( testSpecialPoissonMatrix testSpecialPoissonMatrix.cpp)
    
    add_executable( testLagrange testLagrange.cpp)
    
ENDIF ( ${AMREX_DIM} EQUAL 2 )

# add_executable( testAmrMultigrid
#                 testAmrMultigrid.cpp
#                 ${F90_source_files})

# add_executable( testSolverComparison
#                 testSolverComparison.cpp
#                 ../Solver.cpp
#                 ${F90_source_files})

IF ( ${AMREX_DIM} EQUAL 3 )
    add_executable( testCopy testCopy.cpp)

    target_link_libraries (testCopy ${TESTLIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    add_executable( testBitPattern testBitPattern.cpp)
    
    target_link_libraries (testBitPattern ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    add_executable( trilinos2ippl trilinos2ippl.cpp)
    
    target_link_libraries (trilinos2ippl ${TESTLIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    set( EXECS
         ${EXECS}
         testBitPattern
         trilinos2ippl
         testCopy
    )
ENDIF ( ${AMREX_DIM} EQUAL 3 )
                        

IF ( ${AMREX_DIM} EQUAL 2 )
    target_link_libraries (testRestrictionMatrix ${TESTLIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testInterpolationMatrix ${TESTLIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testBoundaryMatrix ${TESTLIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testSmootherMatrix ${TESTLIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testPoissonMatrix ${TESTLIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testNoFine ${TESTLIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testWithFine ${TESTLIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testSpecialPoissonMatrix ${TESTLIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testLagrange ${TESTLIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    set( EXECS
         ${EXECS}
         testRestrictionMatrix
         testInterpolationMatrix
         testBoundaryMatrix
         testSmootherMatrix
         testPoissonMatrix
         testNoFine
         testWithFine
         testSpecialPoissonMatrix
         testLagrange
    )
ENDIF ( ${AMREX_DIM} EQUAL 2 )

# target_link_libraries (testAmrMultigrid ${TESTLIBS}
#                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
#                            ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

# target_link_libraries (testSolverComparison ${TESTLIBS}
#                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
#                            ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})


set( EXECS
     ${EXECS}
#      testAmrMultigrid
#      testSolverComparison
)

install (TARGETS ${EXECS} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/libexec/amr/trilinos")
