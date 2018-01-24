FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test AMR found in ${_relPath}")
#     MESSAGE (STATUS "MPI COMPILER ${MPI_UNDERLYING_COMPILER}")
    
set (CMAKE_CXX_FLAGS
    "-O3 ${CMAKE_CXX_FLAGS} -DBL_USE_FORTRAN_MPI -DUNIQUE_PTR -DHAVE_AMR_MG_SOLVER"
    )
    
set (CMAKE_CXX_FLAGS
    "${IPPL_CMAKE_CXX_FLAGS} -std=c++11 -DUSEH5FEDV2 -DPARALLEL_IO -DIPPL_AMR ${CMAKE_CXX_FLAGS}"
)

add_definitions(${AMREX_DEFINES})

# FIXME Why do we need this?
# remove digits from hostname: edison03 -> edison
STRING(REGEX MATCH "[^0-9]*" HOSTNAME_BASE "${HOSTNAME}")

IF (${HOSTNAME_BASE} MATCHES "daint")
    set(MPI_CXX_LIBRARIES "-lmpich -lmpichcxx -lfmpich -lmpichf90")
ENDIF (${HOSTNAME_BASE} MATCHES "daint")
    
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
    ${HDF5_LIBRARY}
    ${IPPL_LIBRARY}
    ippl
    ${Boost_LIBRARIES}
    ${Trilinos_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
    ${Fortran_LIBS}
    ${AMREX_LIBRARIES}
    m
    z
)

set(F90_source_files ../Tagging_nd.f90)

set(EXECS "")

# add_executable( testBoundaryVector testBoundaryVector.cpp)

# target_link_libraries (testBoundaryVector ${LIBS}
#                        ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
#                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

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

add_executable( testAmrMultigrid
                testAmrMultigrid.cpp
                ${F90_source_files}
                ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                ${CMAKE_SOURCE_DIR}/src/Classic/Utilities/ClassicException.cpp
                ${CMAKE_SOURCE_DIR}/src/Utilities/OpalException.cpp
                AmrMultiGridCore.h
                BelosBottomSolver.cpp
                AmesosBottomSolver.cpp
                AmrSmoother.cpp
                AmrMultiGrid.cpp)

add_executable( testSolverComparison
                testSolverComparison.cpp
                ../Solver.cpp
                ${F90_source_files}
                ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                ${CMAKE_SOURCE_DIR}/src/Classic/Utilities/ClassicException.cpp
                ${CMAKE_SOURCE_DIR}/src/Utilities/OpalException.cpp
                BelosBottomSolver.cpp
                AmesosBottomSolver.cpp
                AmrSmoother.cpp
                AmrMultiGrid.cpp)

# add_executable( testTrilinosGeometry
#                 testTrilinosGeometry.cpp
#                 ../AmrOpal.cpp
#                 ../Solver.cpp
#                 ${F90_source_files}
#                 ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
#                 EllipticDomain.cpp
#                 AmrTrilinos.cpp)

IF ( ${AMREX_DIM} EQUAL 3 )
    add_executable( testCopy testCopy.cpp)

    target_link_libraries (testCopy ${LIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
#     add_executable( testTrilinosSolver
#                 testTrilinosSolver.cpp
#                 ../AmrOpal.cpp
#                 ../Solver.cpp
#                 ${F90_source_files}
#                 ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
#                 EllipticDomain.cpp
#                 AmrMultiGrid.cpp
#                 AmrTrilinos.cpp)
    

# target_link_libraries (testTrilinosSolver ${LIBS}
#                         ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
#                         ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    add_executable( testBitPattern testBitPattern.cpp)
    
    
    target_link_libraries (testBitPattern ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    set( EXECS
         ${EXECS}
         testBitPattern
         testCopy
    )
ENDIF ( ${AMREX_DIM} EQUAL 3 )
                        

IF ( ${AMREX_DIM} EQUAL 2 )
    target_link_libraries (testRestrictionMatrix ${LIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testInterpolationMatrix ${LIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testBoundaryMatrix ${LIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testSmootherMatrix ${LIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testPoissonMatrix ${LIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testNoFine ${LIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testWithFine ${LIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testSpecialPoissonMatrix ${LIBS}
                           ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testLagrange ${LIBS}
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
# target_link_libraries (testTrilinosGeometry ${LIBS}
#                         ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
#                         ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testAmrMultigrid ${LIBS}
                          ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testSolverComparison ${LIBS}
                          ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})


set( EXECS
     ${EXECS}
     testAmrMultigrid
     testSolverComparison
)

install (TARGETS ${EXECS} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/libexec/amr/trilinos")
