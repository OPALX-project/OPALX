FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test AMR found in ${_relPath}")
#     MESSAGE (STATUS "MPI COMPILER ${MPI_UNDERLYING_COMPILER}")
    
set (CMAKE_CXX_FLAGS
    "-O3 ${CMAKE_CXX_FLAGS} -DUNIQUE_PTR -funroll-all-loops -malign-double -ffast-math"
    )
    
set (CMAKE_CXX_FLAGS
    "${IPPL_CMAKE_CXX_FLAGS} -std=c++11 -DPARALLEL_IO -DIPPL_AMR ${CMAKE_CXX_FLAGS}"
)

add_definitions(${AMREX_DEFINES})

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
    ${Boost_INCLUDE_DIRS}
    ${BOOST_INCLUDE_DIR}
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIR}
)
    
link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${IPPL_LIBRARY_DIR}
    ${AMREX_LIBRARY_DIR}
    ${HYPRE_LIBRARY_DIR}
    ${Boost_LIBRARY_DIRS}
    ${BOOST_LIBRARY_DIR}
    ${Trilinos_LIBRARY_DIRS}
    ${Trilinos_TPL_LIBRARY_DIRS}
    ${GSL_LIBRARY_DIR}
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
    ${AMREX_EXTRA_Fortran_LINK_LINE}
    ${AMREX_LIBRARIES}
    ${GSL_LIBRARY}
    m
    z
)

set(F90_source_files ../Tagging_nd.f90)

add_library(amropal
            ../AmrOpal.cpp
            ../AmrOpal_F.h
            AmrParticleBase1.hpp
            ../Distribution.cpp
            ../H5Reader.cpp
            ${F90_source_files}
            ../Solver.cpp
            ../MGTSolver.cpp)

target_link_libraries(amropal ${LIBS})

set (TESTLIBS amropal)

if ( ENABLE_AMR_MG_SOLVER )
    set (CMAKE_CXX_FLAGS "-DHAVE_AMR_MG_SOLVER ${CMAKE_CXX_FLAGS}")
    set(TESTLIBS ${TESTLIBS} multigrid)
endif ( ENABLE_AMR_MG_SOLVER )

set (EXECS "")

IF ( ${AMREX_DIM} EQUAL 3 )
    # set_property(SOURCE testAmrPartBase.cpp testPerformance.cpp testPlasma.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " -DUNIQUE_PTR ")
    add_executable(testAmrPartBase testAmrPartBase.cpp AmrParticleBase1.hpp)
    
    add_executable( testNewTracker
                    testNewTracker.cpp)
    
    add_executable( testPerformance
                    testPerformance.cpp)
                    
    add_executable( testPlasma
                    PhaseDist.cpp
                    PlasmaPIC.cpp
                    testPlasma.cpp)
                    
    add_executable( testTagging
                    testTagging.cpp)
                    
    add_executable( testInitialBox
                    testInitialBox.cpp)
                    
    add_executable( testDomainTransformSolve
                    testDomainTransformSolve.cpp)
    
    add_executable( testGaussian
                    testGaussian.cpp)
                    
    add_executable( testUnifSphere
                    testUnifSphere.cpp)
    
    add_executable( testThreeGaussian
                    testThreeGaussian.cpp)
        
    add_executable( testFromH5
                    testFromH5.cpp)
    
    add_executable( testMultipleSources
                    testMultipleSources.cpp)
        
    add_executable( testUnifBox
                    testUnifBox.cpp)
    
    add_executable( testLayout
                    testLayout.cpp)

# add_executable( testInitGuessSolver
#                 testInitGuessSolver.cpp
#                 ../AmrOpal.cpp
#                 AmrParticleBase1.hpp
#                 ../Distribution.cpp
#                 ../H5Reader.cpp
#                 ${F90_source_files}
#                 ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
#                 ../Solver.cpp)

    add_executable(testScatterAMReX testScatterAMReX.cpp AmrParticleBase1.hpp)
    
    add_custom_command(TARGET testScatterAMReX POST_BUILD
                       COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/inputs $<TARGET_FILE_DIR:testScatterAMReX>)
    
    
    target_link_libraries (testAmrPartBase ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testNewTracker ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testPerformance ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testPlasma ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testTagging ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

# target_link_libraries (testInitGuessSolver ${TESTLIBS}
#                         ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testInitialBox ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testScatterAMReX ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testDomainTransformSolve ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testGaussian ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testUnifSphere ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testThreeGaussian ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testFromH5 ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testMultipleSources ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testUnifBox ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testLayout ${TESTLIBS}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})


set (EXECS
    testAmrPartBase
    testNewTracker
    testPerformance
    testPlasma
    testTagging
#    testInitGuessSolver
    testInitialBox
    testScatterAMReX
    testDomainTransformSolve
    testGaussian
    testUnifSphere
    testThreeGaussian
    testFromH5
    testMultipleSources
    testUnifBox
    testLayout
)

ENDIF ( ${AMREX_DIM} EQUAL 3 )

IF ( ${AMREX_DIM} EQUAL 2 )
    add_executable( testPlasma
                    PhaseDist.cpp
                    PlasmaPIC.cpp
                    testPlasma.cpp)
    
    target_link_libraries (testPlasma ${TESTLIBS} -lgslcblas -lboost_mpi
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    set (EXECS testPlasma)
ENDIF ( ${AMREX_DIM} EQUAL 2 )
    
install (TARGETS ${EXECS} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/libexec/amr/ippl-amr")
