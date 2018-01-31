IF ( ENABLE_AMR )
#     add_subdirectory(boxlib-amr) still using BoxLib
    add_subdirectory(amrex-only)
    add_subdirectory(ippl-amr)
    add_subdirectory(ippl-bunch)
    IF ( ENABLE_AMR_MG_SOLVER )
    
        add_subdirectory(trilinos)
        
        add_definitions(${AMREX_DEFINES})
        
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
                              ${Trilinos_LIBRARIES}
                              ${Trilinos_TPL_LIBRARIES}
                              ${AMREX_LIBRARIES}
                              ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
                              
    ENDIF ( ENABLE_AMR_MG_SOLVER )
ENDIF (ENABLE_AMR )
