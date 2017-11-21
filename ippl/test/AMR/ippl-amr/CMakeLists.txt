FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test AMR found in ${_relPath}")
#     MESSAGE (STATUS "MPI COMPILER ${MPI_UNDERLYING_COMPILER}")
    
find_package (MPI REQUIRED)
#     message (STATUS "The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
#     message (STATUS "The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
message (STATUS "The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
#     message (STATUS "The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")
    
set (CMAKE_CXX_FLAGS
    "-O3 ${CMAKE_CXX_FLAGS} -DUNIQUE_PTR -funroll-all-loops -malign-double -ffast-math"
    )
    
set (CMAKE_CXX_FLAGS
    "${IPPL_CMAKE_CXX_FLAGS} -std=c++11 -DUSEH5FEDV2 -DPARALLEL_IO -DIPPL_AMR ${CMAKE_CXX_FLAGS}"
)

add_definitions(${AMREX_DEFINES})

if ( ${ENABLE_AMR_MG_SOLVER} )
    set (CMAKE_CXX_FLAGS "-DHAVE_AMR_MG_SOLVER ${CMAKE_CXX_FLAGS}")
endif ()

set (BOOSTROOT $ENV{BOOST_DIR})
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package (Boost 1.54.0 REQUIRED COMPONENTS filesystem system iostreams
regex)

if (Boost_INCLUDE_DIRS)
    message (STATUS "Found boost include dir: ${Boost_INCLUDE_DIR}")
    message (STATUS "Found boost library dir: ${Boost_LIBRARY_DIR}")
    message (STATUS "Found boost libraries: ${Boost_LIBRARIES}")
    include_directories (${Boost_INCLUDE_DIRS})
endif ()
    
include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/ippl/src
    ${CMAKE_SOURCE_DIR}/src/Classic/
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
    
# remove digits from hostname: edison03 -> edison
STRING(REGEX MATCH "[^0-9]*" HOSTNAME_BASE "${HOSTNAME}")
    
IF(${HOSTNAME_BASE} MATCHES "edison" OR ${HOSTNAME_BASE} MATCHES "cori" OR
    ${HOSTNAME_BASE} MATCHES "daint")
    SET (LIBS
        ${H5Hut_LIBRARY}
        ${HDF5_LIBRARY}
        ${IPPL_LIBRARY}
        ippl
	${Boost_LIBRARIES}
        m
        z
    )
ELSE()
    SET (LIBS
        ${H5Hut_LIBRARY}
        ${HDF5_LIBRARY}
        ${IPPL_LIBRARY}
        ippl
	${Boost_LIBRARIES}
        mpi
        mpi_cxx
        mpi_mpifh
        m
        z
    )
ENDIF(${HOSTNAME_BASE} MATCHES "edison" OR ${HOSTNAME_BASE} MATCHES "cori" OR
    ${HOSTNAME_BASE} MATCHES "daint")


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
                
if ( ${ENABLE_AMR_MG_SOLVER} )
    add_executable( testUnifSphere
                    testUnifSphere.cpp
                    ../AmrOpal.cpp
                    AmrParticleBase.hpp
                    ../Distribution.cpp
                    ../H5Reader.cpp
                    ${F90_source_files}
                    ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                    ../MGTSolver.cpp
                    ../Solver.cpp
                    ../trilinos/AmrMultiGridCore.h
                    ../trilinos/BelosBottomSolver.cpp
                    ../trilinos/AmesosBottomSolver.cpp
                    ../trilinos/AmrSmoother.cpp
                    ../trilinos/AmrMultiGrid.cpp)
    
    add_executable( testReal
                    testReal.cpp
                    ../AmrOpal.cpp
                    AmrParticleBase.hpp
                    ../Distribution.cpp
                    ../H5Reader.cpp
                    ${F90_source_files}
                    ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                    ../MGTSolver.cpp
                    ../Solver.cpp
                    ../trilinos/BelosBottomSolver.cpp
                    ../trilinos/AmesosBottomSolver.cpp
                    ../trilinos/AmrSmoother.cpp
                    ../trilinos/AmrMultiGrid.cpp)
    
    add_executable( testMultipleSources
                    testMultipleSources.cpp
                    ../AmrOpal.cpp
                    AmrParticleBase.hpp
                    ../Distribution.cpp
                    ../H5Reader.cpp
                    ${F90_source_files}
                    ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                    ../MGTSolver.cpp
                    ../Solver.cpp
                    ../trilinos/AmrMultiGridCore.h
                    ../trilinos/BelosBottomSolver.cpp
                    ../trilinos/AmesosBottomSolver.cpp
                    ../trilinos/AmrSmoother.cpp
                    ../trilinos/AmrMultiGrid.cpp)
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
    
    add_executable( testReal
                    testReal.cpp
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
    
endif ( ${ENABLE_AMR_MG_SOLVER} )

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
                        -lgfortran ${AMREX_LIBRARIES}
                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testNewTracker ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testPerformance ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testPlasma ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS}
                        -lboost_system -lboost_filesystem)

target_link_libraries (testTagging ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

# target_link_libraries (testInitGuessSolver ${LIBS}
#                         -lgfortran ${AMREX_LIBRARIES}
#                         ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testInitialBox ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testScatterAMReX ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testDomainTransformSolve ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testGaussian ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testUnifSphere ${LIBS}
                        -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testReal ${LIBS}
                        -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testMultipleSources ${LIBS}
                        -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testLayout ${LIBS}
                        -lgfortran ${AMREX_LIBRARIES}
                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
