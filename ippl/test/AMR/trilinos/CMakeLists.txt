FILE (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE (STATUS "Adding test AMR found in ${_relPath}")
#     MESSAGE (STATUS "MPI COMPILER ${MPI_UNDERLYING_COMPILER}")
    
find_package (MPI REQUIRED)
#     message (STATUS "The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
#     message (STATUS "The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
message (STATUS "The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
#     message (STATUS "The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")
    
set (CMAKE_CXX_FLAGS
    "-O3 ${CMAKE_CXX_FLAGS} -DBL_USE_FORTRAN_MPI -DUNIQUE_PTR"
    )
    
set (CMAKE_CXX_FLAGS
    "${IPPL_CMAKE_CXX_FLAGS} -std=c++11 -DUSEH5FEDV2 -DPARALLEL_IO -DIPPL_AMR ${AMREX_DEFINES} ${CMAKE_CXX_FLAGS}"
)

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


set(F90_source_files ../Tagging_nd.f90 electrostatic_pic_3d.f90)

# add_executable( testBoundaryVector testBoundaryVector.cpp)

# target_link_libraries (testBoundaryVector ${LIBS}
#                        -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
#                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

IF ( ${BL_SPACEDIM} EQUAL 2 )

    add_executable( testRestrictionMatrix testRestrictionMatrix.cpp)
    
    add_executable( testInterpolationMatrix testInterpolationMatrix.cpp)
    
    add_executable( testBoundaryMatrix testBoundaryMatrix.cpp)
    
    add_executable( testSmootherMatrix testSmootherMatrix.cpp)
    
    add_executable( testPoissonMatrix testPoissonMatrix.cpp)
    
    add_executable( testNoFine testNoFine.cpp)
    
    add_executable( testWithFine testWithFine.cpp)
    
    add_executable( testSpecialPoissonMatrix testSpecialPoissonMatrix.cpp)
    
    add_executable( testLagrange testLagrange.cpp)
    
    add_executable( testAmrMultigrid
                testAmrMultigrid.cpp
                ../Solver.cpp
                ${F90_source_files}
                ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                AmrMultiGrid.cpp)
    
ENDIF ( ${BL_SPACEDIM} EQUAL 2 )

# add_executable( testTrilinosGeometry
#                 testTrilinosGeometry.cpp
#                 ../AmrOpal.cpp
#                 ../Solver.cpp
#                 ${F90_source_files}
#                 ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
#                 EllipticDomain.cpp
#                 AmrTrilinos.cpp)

IF ( ${BL_SPACEDIM} EQUAL 3 )
    add_executable( testCopy testCopy.cpp)

    target_link_libraries (testCopy ${LIBS}
                           -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    add_executable( testComplement testComplement.cpp)
    
    target_link_libraries (testComplement ${LIBS}
                          -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                          ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    add_executable( testTrilinosSolver
                testTrilinosSolver.cpp
                ../AmrOpal.cpp
                ../Solver.cpp
                ${F90_source_files}
                ${CMAKE_SOURCE_DIR}/src/Classic/Physics/Physics.cpp
                EllipticDomain.cpp
                AmrMultiGrid.cpp
                AmrTrilinos.cpp)

target_link_libraries (testTrilinosSolver ${LIBS}
                        -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
ENDIF ( ${BL_SPACEDIM} EQUAL 3 )
                        

IF ( ${BL_SPACEDIM} EQUAL 2 )
    target_link_libraries (testRestrictionMatrix ${LIBS}
                           -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    target_link_libraries (testInterpolationMatrix ${LIBS}
                           -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testBoundaryMatrix ${LIBS}
                           -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testSmootherMatrix ${LIBS}
                           -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testPoissonMatrix ${LIBS}
                           -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testNoFine ${LIBS}
                           -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testWithFine ${LIBS}
                           -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testSpecialPoissonMatrix ${LIBS}
                           -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testLagrange ${LIBS}
                           -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
                           
    target_link_libraries (testAmrMultigrid ${LIBS}
                          -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
                           ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

ENDIF ( ${BL_SPACEDIM} EQUAL 2 )
# target_link_libraries (testTrilinosGeometry ${LIBS}
#                         -lgfortran ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${AMREX_LIBRARIES}
#                         ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
