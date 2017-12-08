FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test AMR found in ${_relPath}")
#     MESSAGE (STATUS "MPI COMPILER ${MPI_UNDERLYING_COMPILER}")
    
find_package (MPI REQUIRED)
#     message (STATUS "The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
#     message (STATUS "The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
message (STATUS "The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
#     message (STATUS "The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")
    
set (CMAKE_CXX_FLAGS
    "-O3 ${CMAKE_CXX_FLAGS}"
    )
    
set (CMAKE_CXX_FLAGS
    "${IPPL_CMAKE_CXX_FLAGS} -std=c++11 ${AMREX_DEFINES}"
)

include_directories (
    ${AMREX_INCLUDE_DIR}
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
    
IF(${HOSTNAME_BASE} MATCHES "edison" OR ${HOSTNAME_BASE} MATCHES "cori")
    SET (LIBS
        m
        z
    )
ELSE()
    SET (LIBS
        ${MPI_LIBRARIES}
        mpi_mpifh
        m
        z
    )
ENDIF(${HOSTNAME_BASE} MATCHES "edison" OR ${HOSTNAME_BASE} MATCHES "cori")

if ( ${BL_SPACEDIM} EQUAL 3 )
    add_executable( testDeposition
                    testDeposition.cpp)
    
    add_executable( AssignMultiLevelDensity
                    AssignMultiLevelDensity.cpp)

    add_executable( testGeometry
                    testGeometry.cpp)
    
    add_executable( testDepositionFail
                    testDepositionFail.cpp
                    ../Tagging_nd.f90
                    ../AmrOpal_F.h)
    
    add_executable( testSolverFail
                    testSolverFail.cpp
                    ../Tagging_nd.f90
                    ../AmrOpal_F.h
                    Solver.cpp)
    
    # 4. March 2017,
    # http://stackoverflow.com/questions/13429656/cmake-how-to-copy-contents-of-a-directory-into-build-directory-after-make
    add_custom_command(TARGET testDeposition testDepositionFail POST_BUILD
                       COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/inputs $<TARGET_FILE_DIR:testDeposition>)

    target_link_libraries (testDeposition ${LIBS}
                           -lgfortran ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
                        
    target_link_libraries (AssignMultiLevelDensity ${LIBS}
                           -lgfortran ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testGeometry ${LIBS}
                           -lgfortran ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testDepositionFail ${LIBS}
                           -lgfortran ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testSolverFail ${LIBS}
                           -lgfortran ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
endif()
