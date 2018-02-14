FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test AMR found in ${_relPath}")
#     MESSAGE (STATUS "MPI COMPILER ${MPI_UNDERLYING_COMPILER}")
    
find_package (MPI REQUIRED)
#     message (STATUS "The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
#     message (STATUS "The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
message (STATUS "The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
#     message (STATUS "The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")
    
set (CMAKE_CXX_FLAGS "-O3 ${CMAKE_CXX_FLAGS} -DBL_USE_FORTRAN_MPI")
    
set (CMAKE_CXX_FLAGS
    "${IPPL_CMAKE_CXX_FLAGS} -std=c++11 -DPARALLEL_IO ${CMAKE_CXX_FLAGS}"
)
    
include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/ippl/src
    ${CMAKE_SOURCE_DIR}/src/Classic/
    ${IPPL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${AMREX_INCLUDE_DIR}
    ${HYPRE_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
    
link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${IPPL_LIBRARY_DIR}
    ${AMREX_LIBRARY_DIR}
    ${HYPRE_LIBRARY_DIR}
)
    
# remove digits from hostname: edison03 -> edison
STRING(REGEX MATCH "[^0-9]*" HOSTNAME_BASE "${HOSTNAME}")
    
IF(${HOSTNAME_BASE} MATCHES "edison" OR ${HOSTNAME_BASE} MATCHES "cori" OR
   ${HOSTNAME_BASE} MATCHES "daint")
    SET (LIBS
        ${H5Hut_LIBRARY}
        ${HDF5_LIBRARY}
        ${IPPL_LIBRARY}
        ippl
        m
        z
    )
ELSE()
    SET (LIBS
        ${H5Hut_LIBRARY}
        ${HDF5_LIBRARY}
        ${IPPL_LIBRARY}
        ippl
        ${MPI_LIBRARIES}
        mpi_mpifh
        m
        z
    )
ENDIF(${HOSTNAME_BASE} MATCHES "edison" OR ${HOSTNAME_BASE} MATCHES "cori" OR
      ${HOSTNAME_BASE} MATCHES "daint")
    
set(F90_source_files ../Tagging_nd.f90)
    
#     add_executable (testAmrPartBunch testAmrPartBunch.cpp AmrOpal.cpp AmrPartBunch.cpp ../H5Reader.cpp Distribution.cpp ${F90_source_files} ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/Classic/Physics/Physics.cpp Solver.cpp)
#     add_executable (testSolver testSolver.cpp ../H5Reader.cpp Distribution.cpp Solver.cpp)
add_executable (iterative iterative.cpp)
#     add_executable (testError error.cpp AmrPartBunch.cpp ../H5Reader.cpp Distribution.cpp PoissonProblems.cpp Solver.cpp AmrOpal.cpp ${F90_source_files} ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/Classic/Physics/Physics.cpp)
add_executable (testGaussian
                testGaussian.cpp
                ../AmrOpal.cpp AmrPartBunch.cpp
                ../H5Reader.cpp ../Distribution.cpp
                ${F90_source_files}
                ../Solver.cpp)
    
add_executable (testUniform
                testUniform.cpp
                ../AmrOpal.cpp AmrPartBunch.cpp
                ../H5Reader.cpp ../Distribution.cpp
                ${F90_source_files}
                ../Solver.cpp)
    
add_executable (testUnifSphere
                testUnifSphere.cpp
                ../AmrOpal.cpp AmrPartBunch.cpp
                ${F90_source_files} 
                ../Solver.cpp)# HypreABecLap.H HypreABecLap.cpp)
    
add_executable (testGridSolve
                testGridSolve.cpp
                ${F90_source_files}
                ../Solver.cpp)
    
add_executable (testUnifSphereGrid
                testUnifSphereGrid.cpp
                ../AmrOpal.cpp AmrPartBunch.cpp
                ${F90_source_files}
                ../Solver.cpp)
    
    
add_executable(testReal
                testReal.cpp
                ../H5Reader.cpp ../Distribution.cpp
                ../AmrOpal.cpp AmrPartBunch.cpp
                ${F90_source_files}
                ../Solver.cpp)
    
add_executable(testMultiBunch
                testMultiBunch.cpp
                ../H5Reader.cpp ../Distribution.cpp
                ../AmrOpal.cpp AmrPartBunch.cpp
                ${F90_source_files}
                ../Solver.cpp)
    
add_executable (testH5Read testH5Read.cpp ../H5Reader.cpp)

add_executable (testHdf5 testHdf5.cpp)
    
add_executable(toFile toFile.cpp ../H5Reader.cpp ../Distribution.cpp)
    
add_executable (testTracker
                testTracker.cpp
                ../AmrOpal.cpp AmrPartBunch.cpp
                ../H5Reader.cpp ../Distribution.cpp
                ${F90_source_files}
                ../Solver.cpp)
    
add_executable (testOldPlasma
                testOldPlasma.cpp
                ../AmrOpal.cpp AmrPartBunch.cpp
                ../H5Reader.cpp ../Distribution.cpp
                ${F90_source_files}
                ../Solver.cpp)
    
add_executable (testPerfectLoadBalancing
                testPerfectLoadBalancing.cpp
                ../AmrOpal.cpp AmrPartBunch.cpp
                ../H5Reader.cpp ../Distribution.cpp
                ${F90_source_files}
                ../Solver.cpp)

add_executable (testOldPerformance
                testOldPerformance.cpp
                ../AmrOpal.cpp AmrPartBunch.cpp
                ../H5Reader.cpp ../Distribution.cpp
                ${F90_source_files}
                ../Solver.cpp)
    
    # Boxlib has a circular dependency between -lcboxlib and -lfboxlib --> resolve by: -lcboxlib -lfboxlib -lcboxlib
#     target_link_libraries (testAmrPartBunch ${LIBS} -lgfortran ${AMREX_LIBRARIES} -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
#     target_link_libraries (testSolver ${LIBS} -lgfortran ${AMREX_LIBRARIES} -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
#     target_link_libraries (testError ${LIBS} -lgfortran ${AMREX_LIBRARIES} -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (iterative)

target_link_libraries (testH5Read ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testHdf5 ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
target_link_libraries (testGaussian ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
target_link_libraries (testUniform ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
target_link_libraries (testUnifSphere ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
                        
target_link_libraries (testGridSolve ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testUnifSphereGrid ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testReal ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
target_link_libraries (testMultiBunch ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
target_link_libraries (toFile ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
target_link_libraries (testTracker ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
target_link_libraries (testOldPlasma ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS}
                        -lboost_system -lboost_filesystem)
    
target_link_libraries (testPerfectLoadBalancing ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testOldPerformance ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

