IF (ENABLE_AMR_SOLVER)
    FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    MESSAGE (STATUS "Adding test AMR found in ${_relPath}")
#     MESSAGE (STATUS "MPI COMPILER ${MPI_UNDERLYING_COMPILER}")
    
    find_package (MPI REQUIRED)
#     message (STATUS "The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
#     message (STATUS "The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
    message (STATUS "The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
#     message (STATUS "The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")
    
    set (CMAKE_CXX_FLAGS
        "-O3 ${CMAKE_CXX_FLAGS} -DBL_USE_FORTRAN_MPI"
        )
    
    set (CMAKE_CXX_FLAGS
        "${IPPL_CMAKE_CXX_FLAGS} -std=c++11 -DUSEH5FEDV2 -DPARALLEL_IO ${CMAKE_CXX_FLAGS}"
    )
    
    include_directories (
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/ippl/src
        ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/Classic/
        ${IPPL_INCLUDE_DIR}
        ${H5Hut_INCLUDE_DIR}
        ${HDF5_INCLUDE_DIR}
        ${CCSE_INCLUDE_DIR}
        ${HYPRE_INCLUDE_DIR}
    )
    
    link_directories (
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src
        ${IPPL_LIBRARY_DIR}
        ${CCSE_LIBRARY_DIR}
        ${HYPRE_LIBRARY_DIR}
    )
    
    # remove digits from hostname: edison03 -> edison
    STRING(REGEX MATCH "[^0-9]*" HOSTNAME_BASE "${HOSTNAME}")
    
    IF(${HOSTNAME_BASE} MATCHES "edison" OR ${HOSTNAME_BASE} MATCHES "cori")
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
            mpi
            mpi_cxx
            mpi_mpifh
            m
            z
        )
    ENDIF(${HOSTNAME_BASE} MATCHES "edison" OR ${HOSTNAME_BASE} MATCHES "cori")
    
    set(F90_source_files Tagging_nd.f90)
    
#     add_executable (testAmrPartBunch testAmrPartBunch.cpp AmrOpal.cpp AmrPartBunch.cpp H5Reader.cpp Distribution.cpp ${F90_source_files} ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/Classic/Physics/Physics.cpp Solver.cpp)
#     add_executable (testSolver testSolver.cpp H5Reader.cpp Distribution.cpp Solver.cpp)
    add_executable (iterative iterative.cpp)
#     add_executable (testError error.cpp AmrPartBunch.cpp H5Reader.cpp Distribution.cpp PoissonProblems.cpp Solver.cpp AmrOpal.cpp ${F90_source_files} ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/Classic/Physics/Physics.cpp)
    add_executable (testGaussian testGaussian.cpp AmrOpal.cpp AmrPartBunch.cpp H5Reader.cpp Distribution.cpp ${F90_source_files} ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/Classic/Physics/Physics.cpp Solver.cpp)
    
    add_executable (testUnifSphere testUnifSphere.cpp AmrOpal.cpp AmrPartBunch.cpp ${F90_source_files} ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/Classic/Physics/Physics.cpp Solver.cpp)# HypreABecLap.H HypreABecLap.cpp)
    
    add_executable (testGridSolve testGridSolve.cpp ${F90_source_files} ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/Classic/Physics/Physics.cpp Solver.cpp)
    
    add_executable (testUnifSphereGrid testUnifSphereGrid.cpp AmrOpal.cpp AmrPartBunch.cpp ${F90_source_files} ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/Classic/Physics/Physics.cpp Solver.cpp)
    
    
    add_executable(testReal testReal.cpp H5Reader.cpp Distribution.cpp AmrOpal.cpp AmrPartBunch.cpp ${F90_source_files} ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/Classic/Physics/Physics.cpp Solver.cpp)
    
    add_executable (testH5Read testH5Read.cpp H5Reader.cpp)
    
    add_executable(toFile toFile.cpp H5Reader.cpp Distribution.cpp)
    
    # Boxlib has a circular dependency between -lcboxlib and -lfboxlib --> resolve by: -lcboxlib -lfboxlib -lcboxlib
#     target_link_libraries (testAmrPartBunch ${LIBS} -lgfortran ${CCSE_LIBRARIES} -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
#     target_link_libraries (testSolver ${LIBS} -lgfortran ${CCSE_LIBRARIES} -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
#     target_link_libraries (testError ${LIBS} -lgfortran ${CCSE_LIBRARIES} -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testH5Read ${LIBS} -lgfortran ${CCSE_LIBRARIES} -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testGaussian ${LIBS} -lgfortran ${CCSE_LIBRARIES} -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    target_link_libraries (testUnifSphere ${LIBS} -lgfortran ${CCSE_LIBRARIES} -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    target_link_libraries (testGridSolve ${LIBS} -lgfortran ${CCSE_LIBRARIES} -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    target_link_libraries (testUnifSphereGrid ${LIBS} -lgfortran ${CCSE_LIBRARIES} -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    target_link_libraries (testReal ${LIBS} -lgfortran ${CCSE_LIBRARIES} -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (toFile ${LIBS} -lgfortran ${CCSE_LIBRARIES} -lcboxlib ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
ENDIF(ENABLE_AMR_SOLVER)
