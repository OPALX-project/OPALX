CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (IPPL CXX Fortran)
SET (IPPL_VERSION_MAJOR 1)
SET (IPPL_VERSION_MINOR 0)
SET (PACKAGE \"ippl\")
SET (PACKAGE_BUGREPORT \"andreas.adelmann@psi.ch\")
SET (PACKAGE_NAME \"IPPL\")
SET (PACKAGE_STRING \"IPPL\ 1.0\")
SET (PACKAGE_TARNAME \"ippl\")
SET (PACKAGE_VERSION \"1.0\")
SET (VERSION \"1.0\")

ENABLE_LANGUAGE (Fortran)

OPTION (LINUXGCC "Use Linux GCC flags" ON)
OPTION (XT3GCC "Use XT3 GCC flags" OFF)
OPTION (LINUXPGI "Use Linux PGI flags" OFF)

SET (CMAKE_CXX_COMPILER ${MPI_COMPILER})

SET (IPPL_CXX_FLAGS "")
SET (OTHER_CXX_FLAGS "")

IF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET (CMAKE_BUILD_TYPE RelWithDebInfo)
ENDIF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
SET (CMAKE_CXX_FLAGS_RELEASE "-O2")
SET (CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

### CXX_FLAGS specific to build type ########################
IF (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    SET (OTHER_CXX_FLAGS
        "-DNOPAssert ${OTHER_CXX_FLAGS}"
    )
ELSE (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    SET (OTHER_CXX_FLAGS
        "-DIPPL_DEBUG ${OTHER_CXX_FLAGS}"
    )
ENDIF (${CMAKE_BUILD_TYPE} STREQUAL "Release")


### CXX_FLAGS common to all compilers and platforms which ###
### may be written to the configure file                  ###
SET (IPPL_CXX_FLAGS
    "-std=c++1y -DIPPL_MPI -DMPICH_SKIP_MPICXX -DIPPL_GCC -DIPPL_DONT_POOL -DIPPL_USE_XDIV_RNG -DIPPL_NO_STRINGSTREAM -DPETE_BITWISE_COPY -DIPPL_HAS_TEMPLATED_COMPLEX -DIPPL_USE_STANDARD_HEADERS -DIPPL_USE_PARTIAL_SPECIALIZATION -DIPPL_STDSTL -DIPPL_LONGLONG -DWITH_BRICK -DnoCOMP_GNUOLD -Drestrict=__restrict__ -DNOCTAssert ${IPPL_CXX_FLAGS}"
)


### CXX_FLAGS common to all compilers and platforms which ###
### may *not* be written to the configure file            ###
SET (OTHER_CXX_FLAGS
    "-Wall ${OTHER_CXX_FLAGS}"
)


### CXX_FLAGS specific to compiler and platform #############
IF (LINUXGCC)
    MESSAGE (STATUS "Linux GCC flags SET")
    SET (IPPL_CXX_FLAGS
        "-DIPPL_LINUX -DIPPL_STRINGSTREAM ${IPPL_CXX_FLAGS}"
    )
    SET (OTHER_CXX_FLAGS
        "-ftemplate-depth-80 -funroll-loops -fstrict-aliasing ${OTHER_CXX_FLAGS}"
    )

    # compiler dependent flags
    STRING (REGEX REPLACE ".*/([A-Za-z]*)$" "\\1" COMPILER_NAME ${CMAKE_CXX_COMPILER})
    MESSAGE (STATUS "Your compiler is: ${COMPILER_NAME}")

    IF ("${COMPILER_NAME}" STREQUAL "mpicxx")

        SET (MPI_UNDERLYING_COMPILER
             "")

        EXECUTE_PROCESS (COMMAND ${CMAKE_CXX_COMPILER} -show
                         OUTPUT_VARIABLE MPI_COMPILER_OUTPUT
                         ERROR_VARIABLE MPI_COMPILER_ERROR)

        IF ("${MPI_COMPILER_ERROR}" STREQUAL "")
            STRING (REGEX REPLACE "([A-Za-z]*) .*$" "\\1" MPI_UNDERLYING_COMPILER ${MPI_COMPILER_OUTPUT})
        ENDIF ("${MPI_COMPILER_ERROR}" STREQUAL "")

        IF ("${MPI_UNDERLYING_COMPILER}" STREQUAL "icpc")
            # using intel compiler
            MESSAGE (STATUS "The underlying compiler of ${COMPILER_NAME} is: icpc")
            SET (OTHER_CXX_FLAGS
                 "-diag-disable 383,981 ${OTHER_CXX_FLAGS}")

        ELSEIF ("${MPI_UNDERLYING_COMPILER}" STREQUAL "clang++")
            # using clang compiler
            MESSAGE (STATUS "The underlying compiler of ${COMPILER_NAME} is: clang++")
            SET (OTHER_CXX_FLAGS
                 "-Wall -Wsign-compare -Wunused-variable -Warray-bounds -DIPPL_RESTRICT_BUG ${OTHER_CXX_FLAGS}")

        ELSEIF ("${MPI_UNDERLYING_COMPILER}" STREQUAL "g++")
            # using gnu compiler
            MESSAGE (STATUS "The underlying compiler of ${COMPILER_NAME} is: g++")
            EXECUTE_PROCESS (COMMAND ${MPI_UNDERLYING_COMPILER} --version
                             OUTPUT_VARIABLE GCC_VERSION_OUTPUT
                             ERROR_VARIABLE GCC_VERSION_ERROR)
            STRING (REGEX REPLACE ".*([0-9]\\.[0-9]\\.[0-9]).*" "\\1" GCC_VERSION ${GCC_VERSION_OUTPUT})
            MESSAGE (STATUS "Your gcc version is: ${GCC_VERSION}")

            IF (${GCC_VERSION} VERSION_LESS "4.8.0")
                SET (OTHER_CXX_FLAGS
    	            "-Wall -std=c++0x")

            ELSEIF (${GCC_VERSION} VERSION_GREATER "4.7")
                SET (OTHER_CXX_FLAGS
                    "-Wall -Wno-reorder -std=c++11")
            ENDIF ()
        ELSE ("${MPI_UNDERLYING_COMPILER}" STREQUAL "icpc")
            # using unknown compiler
	    MESSAGE (STATUS "${COMPILER_NAME} is using a compiler (${MPI_UNDERLYING_COMPILER}) we were not thinking of!
                     Please use the gnu compiler or the intel compiler if you are having problems.")
        ENDIF ("${MPI_UNDERLYING_COMPILER}" STREQUAL "icpc")

    ELSE ("${COMPILER_NAME}" STREQUAL "mpicxx")
        # using unknown mpi implementation
        MESSAGE (STATUS "You are using an unsupported MPI compiler: ${COMPILER_NAME}!
                 Please use the OpenMPI or MPICH if you are having problems.")
    ENDIF ("${COMPILER_NAME}" STREQUAL "mpicxx")

    SET (IPPL_VERSION_NAME \"V1.0 01-2007 Linux GNU\")
ENDIF (LINUXGCC)

IF (XT3GCC)
    MESSAGE (STATUS "XT3 GCC flags SET")
    SET (IPPL_CXX_FLAGS
        "-DIPPL_XT3 -DIPPL_DEBUG -DIPPL_LINUX -DIPPL_STRINGSTREAM ${IPPL_CXX_FLAGS}"
    )
    SET (OTHER_CXX_FLAGS
        "-Wno-deprecated -ftemplate-depth-80 -funroll-loops -fstrict-aliasing"
    )
    SET (IPPL_VERSION_NAME \"V1.0 01-2007 XT3 GNU\")
ENDIF (XT3GCC)

IF (LINUXPGI)
    MESSAGE(STATUS "Linux PGI flags SET")
    SET (IPPL_CXX_FLAGS
        "-DIPPLProf -DIPPL_GPLUSPLUS -DMPICH_IGNORE_CXX_SEEK -DIPPL_LINUX ${IPPL_CXX_FLAGS}"
    )
    SET (IPPL_VERSION_NAME \"V1.0 01-2007 Linux PGI\")
ENDIF (LINUXPGI)

SET (IPPLLINUXGCCF
    "-DIPPL_LINUX -DIPPL_STRINGSTREAM"
)

SET (IPPLXT3PGI
    "-DIPPLProf -DIPPL_GPLUSPLUS -DIPPL_XT3 -DMPICH_IGNORE_CXX_SEEK"
)


SET (CMAKE_CXX_FLAGS
    "${OTHER_CXX_FLAGS} ${IPPL_CXX_FLAGS} ${CMAKE_CXX_FLAGS} "
)

### add test directories ####################################
ADD_SUBDIRECTORY (src)

IF (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   ADD_SUBDIRECTORY (test)
ENDIF (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")


### write configure file ####################################
CONFIGURE_FILE ( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake )


### install files ###########################################
INSTALL (
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake
    DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/${PROJECT_NAME}"
    RENAME ${PROJECT_NAME}Config.cmake
)
