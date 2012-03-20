cmake_minimum_required(VERSION 2.6)

set (CMAKE_CXX_FLAGS_RELEASE "-g -O2")
set (CMAKE_CXX_FLAGS_DEBUG  "-g -O2")
#find_package (MPI REQUIRED)
set (CMAKE_CXX_COMPILER ${MPI_COMPILER})

project(IPPL CXX Fortran)
set (IPPL_VERSION_MAJOR 1)
set (IPPL_VERSION_MINOR 0)

option (LINUXGCC "Use Linux GCC flags" ON)
option (XT3GCC "Use XT3 GCC flags" OFF)
option (LINUXPGI "Use Linux PGI flags" OFF)

ENABLE_LANGUAGE(Fortran)

#set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")
#find_package (HDF5 REQUIRED)

set (PACKAGE \"ippl\")
set (PACKAGE_BUGREPORT \"andreas.adelmann@psi.ch\")
set (PACKAGE_NAME \"IPPL\")
set (PACKAGE_STRING \"IPPL\ 1.0\")
set (PACKAGE_TARNAME \"ippl\")
set (PACKAGE_VERSION \"1.0\")
set (VERSION \"1.0\")

set (IPPLCOMMON
    "-DIPPL_MPI -DMPICH_SKIP_MPICXX -DIPPL_GCC -DIPPL_DONT_POOL -DIPPL_USE_XDIV_RNG -DIPPL_NO_STRINGSTREAM -DPETE_BITWISE_COPY -DIPPL_HAS_TEMPLATED_COMPLEX -DIPPL_USE_STANDARD_HEADERS -DIPPL_USE_PARTIAL_SPECIALIZATION -DIPPL_STDSTL -DIPPL_LONGLONG -DWITH_BRICK -DnoCOMP_GNUOLD -Drestrict=__restrict__ -DNOCTAssert"
)

set (IPPLLINUXGCC
    "-DIPPL_DEBUG -DIPPL_LINUX -DIPPL_STRINGSTREAM"
)

set (IPPLLINUXGCCF
    "-w -DIPPL_DEBUG -DIPPL_LINUX -DIPPL_STRINGSTREAM"
)

set (IPPLXT3GCC
    "-DIPPL_XT3 -DIPPL_DEBUG -DIPPL_LINUX -DIPPL_STRINGSTREAM"
)

set (IPPLXT3PGI
    "-O3 -w -DIPPLProf -DIPPL_GPLUSPLUS -DIPPL_XT3 -DMPICH_IGNORE_CXX_SEEK"
)

set (IPPLLINUXPGI
    "-DIPPLProf -DIPPL_GPLUSPLUS -DMPICH_IGNORE_CXX_SEEK -DIPPL_LINUX "
)

IF (LINUXGCC)
    MESSAGE("Linux GCC flags set")
    set (IPPL_CMAKE_CXX_FLAGS 
        "${IPPLLINUXGCC} ${IPPLCOMMON}"
    ) 
    set (OTHER_CMAKE_CXX_FLAGS
        "-Wall -ftemplate-depth-80 -funroll-loops -fstrict-aliasing "
    )
    set (IPPL_VERSION_NAME \"V1.0 01-2007 Linux GNU\")

    # compiler dependent flags
    STRING(REGEX REPLACE ".*/([A-Za-z]*)$" "\\1" COMPILER_NAME ${CMAKE_CXX_COMPILER})
    MESSAGE("Your compiler is: ${COMPILER_NAME}")

    IF ("${COMPILER_NAME}" STREQUAL "mpicxx")

        SET(MPI_UNDERLYING_COMPILER
            "")

        EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} -show
                        OUTPUT_VARIABLE MPI_COMPILER_OUTPUT
                        ERROR_VARIABLE MPI_COMPILER_ERROR)

        IF("${MPI_COMPILER_ERROR}" STREQUAL "")
            STRING(REGEX REPLACE "([A-Za-z]*) .*$" "\\1" MPI_UNDERLYING_COMPILER ${MPI_COMPILER_OUTPUT})
        ENDIF("${MPI_COMPILER_ERROR}" STREQUAL "")

        IF("${MPI_UNDERLYING_COMPILER}" STREQUAL "icpc")
            # using intel compiler
            MESSAGE("The underlying compiler of ${COMPILER_NAME} is: icpc")
            SET (OTHER_CMAKE_CXX_FLAGS
                 "-diag-disable 383 -diag-disable 981 ${OTHER_CMAKE_CXX_FLAGS}")

        ELSEIF("${MPI_UNDERLYING_COMPILER}" STREQUAL "clang++")
            # using clang compiler
            MESSAGE("The underlying compiler of ${COMPILER_NAME} is: clang++")
            SET (OTHER_CMAKE_CXX_FLAGS
                 "-Wsign-compare -Wunused-variable -Warray-bounds -DIPPL_RESTRICT_BUG ${OTHER_CMAKE_CXX_FLAGS}")

        ELSEIF("${MPI_UNDERLYING_COMPILER}" STREQUAL "g++")
            # using gnu compiler
            MESSAGE("The underlying compiler of ${COMPILER_NAME} is: g++")

        ELSE("${MPI_UNDERLYING_COMPILER}" STREQUAL "icpc")
            # using unknown compiler
	    MESSAGE("${COMPILER_NAME} is using a compiler (${MPI_UNDERLYING_COMPILER}) we were not thinking of! 
                     Please use the gnu compiler or the intel compiler if you are having problems.")
        ENDIF("${MPI_UNDERLYING_COMPILER}" STREQUAL "icpc")

    ELSE ("${COMPILER_NAME}" STREQUAL "mpicxx")
        # using unknown mpi implementation
        MESSAGE("You are using an unsupported MPI compiler: ${COMPILER_NAME}!
                 Please use the OpenMPI or MPICH if you are having problems.")
    ENDIF ("${COMPILER_NAME}" STREQUAL "mpicxx")
    
ENDIF (LINUXGCC)

IF (XT3GCC)
    MESSAGE("XT3 GCC flags set")
    set (IPPL_CMAKE_CXX_FLAGS 
        "${IPPLXT3GCC} ${IPPLCOMMON}"
    ) 
    set (OTHER_CMAKE_CXX_FLAGS
        "-Wall -Wno-deprecated -ftemplate-depth-80 -funroll-loops -fstrict-aliasing -g -O3 "
    )
    set (IPPL_VERSION_NAME \"V1.0 01-2007 XT3 GNU\")
ENDIF (XT3GCC)

IF (LINUXPGI)
    MESSAGE("Linux PGI flags set")
    set (IPPL_CMAKE_CXX_FLAGS 
        "${IPPLLINUXPGI} ${IPPLCOMMON}"
    ) 
    set (OTHER_CMAKE_CXX_FLAGS
        "-O3 -w "
    )
    set (IPPL_VERSION_NAME \"V1.0 01-2007 Linux PGI\")
ENDIF (LINUXPGI)

set (CMAKE_CXX_FLAGS
    "${OTHER_CMAKE_CXX_FLAGS} ${IPPL_CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS}"
)
    
add_subdirectory (src)

IF (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   add_subdirectory (test)
ENDIF (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake )

INSTALL(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake
  DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/${PROJECT_NAME}"
  RENAME ${PROJECT_NAME}Config.cmake
)
