IF ( ENABLE_AMR )
    # FIXME Why do we need this?
    # remove digits from hostname: edison03 -> edison
    STRING(REGEX MATCH "[^0-9]*" HOSTNAME_BASE "${HOSTNAME}")
        
    IF (${HOSTNAME_BASE} MATCHES "daint")
        set(MPI_CXX_LIBRARIES "-lmpich -lmpichcxx -lfmpich -lmpichf90")
    ENDIF (${HOSTNAME_BASE} MATCHES "daint")
    
#     add_subdirectory(boxlib-amr) still using BoxLib
    add_subdirectory(amrex-only)
    add_subdirectory(ippl-amr)
    add_subdirectory(ippl-bunch)
    IF ( ENABLE_AMR_MG_SOLVER )
    
        add_subdirectory(trilinos)
        
        add_definitions(${AMREX_DEFINES})
        
        set (CMAKE_CXX_FLAGS "-DHAVE_AMR_MG_SOLVER ${CMAKE_CXX_FLAGS}")
        
        include_directories (
            ${CMAKE_SOURCE_DIR}/src
            ${CMAKE_SOURCE_DIR}/ippl/src
            ${CMAKE_SOURCE_DIR}/src/Classic/
            ${CMAKE_SOURCE_DIR}/src/Classic/Utilities/
            ${CMAKE_SOURCE_DIR}/src/Utilities/
            ${IPPL_INCLUDE_DIR}
            ${AMREX_INCLUDE_DIR}
            ${Trilinos_INCLUDE_DIRS}
            ${Trilinos_TPL_INCLUDE_DIRS}
        )
        
        link_directories (
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_SOURCE_DIR}/src
            ${IPPL_LIBRARY_DIR}
            ${AMREX_LIBRARY_DIR}
            ${Trilinos_LIBRARY_DIRS}
            ${Trilinos_TPL_LIBRARY_DIRS}
        )

        add_library(multigrid
                    trilinos/AmrMultiGridCore.h
                    trilinos/BelosBottomSolver.cpp
                    trilinos/AmesosBottomSolver.cpp
                    trilinos/AmrSmoother.cpp
                    trilinos/AmrMultiGrid.cpp
                    ${CMAKE_SOURCE_DIR}/src/Classic/Utilities/ClassicException.cpp
                    ${CMAKE_SOURCE_DIR}/src/Utilities/OpalException.cpp)
        
        target_link_libraries(multigrid
                              ${MPI_CXX_LIBRARIES}
                              ${MPI_Fortran_LIBRARIES}
                              ${Fortran_LIBS}
                              ${Trilinos_LIBRARIES}
                              ${Trilinos_TPL_LIBRARIES}
                              ${AMREX_EXTRA_Fortran_LINK_LINE}
                              ${AMREX_LIBRARIES})
                              
    ENDIF ( ENABLE_AMR_MG_SOLVER )
ENDIF (ENABLE_AMR )
