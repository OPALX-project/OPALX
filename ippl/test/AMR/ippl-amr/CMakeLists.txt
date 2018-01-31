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

# FIXME Why do we need this?
# remove digits from hostname: edison03 -> edison
STRING(REGEX MATCH "[^0-9]*" HOSTNAME_BASE "${HOSTNAME}")

IF (${HOSTNAME_BASE} MATCHES "daint")
    set(MPI_CXX_LIBRARIES "-lmpich -lmpichcxx -lfmpich -lmpichf90")
ENDIF (${HOSTNAME_BASE} MATCHES "daint")

if ( ENABLE_AMR_MG_SOLVER )
    set (CMAKE_CXX_FLAGS "-DHAVE_AMR_MG_SOLVER ${CMAKE_CXX_FLAGS}")
endif ()

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
    ${MPI_CXX_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
    ${Fortran_LIBS}
    ${AMREX_LIBRARIES}
    m
    z
)

set(F90_source_files ../Tagging_nd.f90)

# set_property(SOURCE testAmrPartBase.cpp testPerformance.cpp testPlasma.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " -DUNIQUE_PTR ")
add_executable(testAmrPartBase testAmrPartBase.cpp AmrParticleBase.hpp)

add_executable( testNewTracker
                testNewTracker.cpp
                ../AmrOpal.cpp
                AmrParticleBase.hpp
                ../Distribution.cpp
                ../H5Reader.cpp
                ${F90_source_files}
                ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                ../Solver.cpp)

add_executable( testPerformance
                testPerformance.cpp
                ../AmrOpal.cpp
                AmrParticleBase.hpp
                ../Distribution.cpp
                ../H5Reader.cpp
                ${F90_source_files}
                ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                ../Solver.cpp)

add_executable( testPlasma
                testPlasma.cpp
                ../AmrOpal.cpp
                AmrParticleBase.hpp
                ../Distribution.cpp
                ../H5Reader.cpp
                ${F90_source_files}
                ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                ../Solver.cpp)

add_executable( testTagging
                testTagging.cpp
                ../AmrOpal.cpp
                AmrParticleBase.hpp
                ../Distribution.cpp
                ../H5Reader.cpp
                ${F90_source_files}
                ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                ../Solver.cpp)
                
add_executable( testInitialBox
                testInitialBox.cpp
                ../AmrOpal.cpp
                AmrParticleBase.hpp
                ../Distribution.cpp
                ../H5Reader.cpp
                ${F90_source_files}
                ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                ../Solver.cpp)
                
                
add_executable( testDomainTransformSolve
                testDomainTransformSolve.cpp
                ../AmrOpal.cpp
                AmrParticleBase.hpp
                ../Distribution.cpp
                ../H5Reader.cpp
                ${F90_source_files}
                ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                ../Solver.cpp)

add_executable( testGaussian
                testGaussian.cpp
                ../AmrOpal.cpp
                AmrParticleBase.hpp
                ../Distribution.cpp
                ../H5Reader.cpp
                ${F90_source_files}
                ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                ../Solver.cpp)
                
if ( ENABLE_AMR_MG_SOLVER )
    set(LIBS ${LIBS} multigrid)
    
    add_executable( testUnifSphere
                    testUnifSphere.cpp
                    ../AmrOpal.cpp
                    AmrParticleBase.hpp
                    ../Distribution.cpp
                    ../H5Reader.cpp
                    ${F90_source_files}
                    ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                    ../MGTSolver.cpp
                    ../Solver.cpp)
    
    add_executable( testFromH5
                    testFromH5.cpp
                    ../AmrOpal.cpp
                    AmrParticleBase.hpp
                    ../Distribution.cpp
                    ../H5Reader.cpp
                    ${F90_source_files}
                    ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                    ../MGTSolver.cpp
                    ../Solver.cpp)
    
    add_executable( testMultipleSources
                    testMultipleSources.cpp
                    ../AmrOpal.cpp
                    AmrParticleBase.hpp
                    ../Distribution.cpp
                    ../H5Reader.cpp
                    ${F90_source_files}
                    ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                    ../MGTSolver.cpp
                    ../Solver.cpp)

    add_executable( testUnifBox
                    testUnifBox.cpp
                    ../AmrOpal.cpp
                    AmrParticleBase.hpp
                    ../Distribution.cpp
                    ../H5Reader.cpp
                    ${F90_source_files}
                    ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                    ../MGTSolver.cpp
                    ../Solver.cpp)
else ()
    add_executable( testUnifSphere
                    testUnifSphere.cpp
                    ../AmrOpal.cpp
                    AmrParticleBase.hpp
                    ../Distribution.cpp
                    ../H5Reader.cpp
                    ${F90_source_files}
                    ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                    ../MGTSolver.cpp
                    ../Solver.cpp)
    
    add_executable( testFromH5
                    testFromH5.cpp
                    ../AmrOpal.cpp
                    AmrParticleBase.hpp
                    ../Distribution.cpp
                    ../H5Reader.cpp
                    ${F90_source_files}
                    ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                    ../MGTSolver.cpp
                    ../Solver.cpp)
    
    add_executable( testMultipleSources
                    testMultipleSources.cpp
                    ../AmrOpal.cpp
                    AmrParticleBase.hpp
                    ../Distribution.cpp
                    ../H5Reader.cpp
                    ${F90_source_files}
                    ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                    ../MGTSolver.cpp
                    ../Solver.cpp)
    
    add_executable( testUnifBox
                    testUnifBox.cpp
                    ../AmrOpal.cpp
                    AmrParticleBase.hpp
                    ../Distribution.cpp
                    ../H5Reader.cpp
                    ${F90_source_files}
                    ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                    ../MGTSolver.cpp
                    ../Solver.cpp)

endif ( ENABLE_AMR_MG_SOLVER )

add_executable( testLayout
                testLayout.cpp
                ../AmrOpal.cpp
                AmrParticleBase.hpp
                ${F90_source_files}
                ../Solver.cpp)

# add_executable( testInitGuessSolver
#                 testInitGuessSolver.cpp
#                 ../AmrOpal.cpp
#                 AmrParticleBase.hpp
#                 ../Distribution.cpp
#                 ../H5Reader.cpp
#                 ${F90_source_files}
#                 ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
#                 ../Solver.cpp)

add_executable(testScatterAMReX testScatterAMReX.cpp AmrParticleBase.hpp)

add_custom_command(TARGET testScatterAMReX POST_BUILD
                   COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/inputs $<TARGET_FILE_DIR:testScatterAMReX>)


target_link_libraries (testAmrPartBase ${LIBS}
                       ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testNewTracker ${LIBS}
                       ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testPerformance ${LIBS}
                       ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testPlasma ${LIBS}
                       ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testTagging ${LIBS}
                       ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

# target_link_libraries (testInitGuessSolver ${LIBS}
#                         ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testInitialBox ${LIBS}
                       ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testScatterAMReX ${LIBS}
                       ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testDomainTransformSolve ${LIBS}
                       ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testGaussian ${LIBS}
                       ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testUnifSphere ${LIBS}
                       ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testFromH5 ${LIBS}
                       ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testMultipleSources ${LIBS}
                       ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testUnifBox ${LIBS}
                       ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testLayout ${LIBS}
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
    testFromH5
    testMultipleSources
    testUnifBox
    testLayout
)

install (TARGETS ${EXECS} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/libexec/amr/ippl-amr")
