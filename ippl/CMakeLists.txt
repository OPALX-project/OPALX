CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (IPPL CXX Fortran)
SET (IPPL_VERSION_MAJOR 1)
SET (IPPL_VERSION_MINOR 1.4)
SET (PACKAGE \"ippl\")
SET (PACKAGE_BUGREPORT \"andreas.adelmann@psi.ch\")
SET (PACKAGE_NAME \"IPPL\")
SET (PACKAGE_STRING \"IPPL\ 1.1.4\")
SET (PACKAGE_TARNAME \"ippl\")
SET (PACKAGE_VERSION \"1.1.4\")
SET (VERSION \"1.1.4\")

set (IPPL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

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

SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
SET (CMAKE_CXX_FLAGS_RELEASE "-O3")
SET (CMAKE_CXX_FLAGS_DEBUG  "-O3 -g")

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
    "-DIPPL_MPI -DMPICH_SKIP_MPICXX -DIPPL_GCC -DIPPL_DONT_POOL -DIPPL_USE_XDIV_RNG -DIPPL_NO_STRINGSTREAM -DPETE_BITWISE_COPY -DIPPL_HAS_TEMPLATED_COMPLEX -DIPPL_USE_STANDARD_HEADERS -DIPPL_USE_PARTIAL_SPECIALIZATION -DIPPL_STDSTL -DIPPL_LONGLONG -DWITH_BRICK -DnoCOMP_GNUOLD -Drestrict=__restrict__ -DNOCTAssert ${IPPL_CXX_FLAGS}"
)


### CXX_FLAGS common to all compilers and platforms which ###
### may *not* be written to the configure file            ###
SET (OTHER_CXX_FLAGS
    "-Wall ${OTHER_CXX_FLAGS}"
)

### Uldis Locans 2014.09.05 ###
OPTION (ENABLE_OPENCL "Enable OpenCL" ON)
OPTION (ENABLE_CUDA "Enable CUDA" ON)
OPTION (ENABLE_OPENMP "Enable OpenMP + offload" OFF)

OPTION (USE_OPENCL "Use OpenCL" OFF)
OPTION (USE_CUDA "Use CUDA" ON)
OPTION (USE_MIC "Use MIC" OFF)


### Load DKS package ###
IF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_OPENMP)
   FIND_PACKAGE(DKS REQUIRED HINTS $ENV{DKS_PREFIX} $ENV{DKS_DIR} $ENV{DKS})
   MESSAGE (STATUS "Found DKS library: ${DKS_LIBRARY_DIR}")
   MESSAGE (STATUS "Found DKS include dir: ${DKS_INCLUDE_DIR}")
   SET (IPPL_CXX_FLAGS "${IPPL_CXX_FLAGS} -DDKS_MPI")
ENDIF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_OPENMP)

### OpenCL compiler flags ###
IF (ENABLE_OPENCL)
	SET (IPPL_CXX_FLAGS "${IPPL_CXX_FLAGS} -lOpenCL -lpthread -DDKS_OPENCL")
ENDIF (ENABLE_OPENCL)
### End OpenCL ###

### CUDA compiler flags ###
IF (ENABLE_CUDA)
    SET (IPPL_CXX_FLAGS "${IPPL_CXX_FLAGS} -lcudart -lcufft -lcublas -lnvToolsExt -DDKS_CUDA")
ENDIF (ENABLE_CUDA)
### End Cuda ###

### MIC compiler flags ###
IF (ENABLE_OPENMP)
    SET (IPPL_CXX_FLAGS "${IPPL_CXX_FLAGS} -openmp -mkl -offload -DDKS_MIC")

    FIND_PROGRAM(XIAR xiar)	
    IF(XIAR)
	SET(CMAKE_AR "${XIAR}")
    ENDIF(XIAR)
    MARK_AS_ADVANCED(XIAR)
    SET(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> rcs -qoffload-build <TARGET> <LINK_FLAGS> <OBJECTS>")
    SET(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> rcs -qoffload-build <TARGET> <LINK_FLAGS> <OBJECTS>")

    FIND_PROGRAM(XILD xild)
    IF(XILD)
	SET(CMAKE_LINKER "${XILD}")
    ENDIF(XILD)
    MARK_AS_ADVANCED(XILD)

ENDIF (ENABLE_OPENMP)
### END MIC ###

### if any accelerator enabled set flag to use DKS ###
IF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_OPENMP)
	SET (IPPL_CXX_FLAGS "${IPPL_CXX_FLAGS} -DIPPL_DKS")
ENDIF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_OPENMP)

### select API for DKS (temporary - will be moved to runtime) ###
IF (USE_OPENCL)
	SET (IPPL_CXX_FLAGS "${IPPL_CXX_FLAGS} -DIPPL_DKS_OPENCL")
ENDIF (USE_OPENCL)
IF (USE_CUDA)
	SET (IPPL_CXX_FLAGS "${IPPL_CXX_FLAGS} -DIPPL_DKS_CUDA")
ENDIF (USE_CUDA)
IF (USE_MIC)
	SET (IPPL_CXX_FLAGS "${IPPL_CXX_FLAGS} -DIPPL_DKS_MIC")
ENDIF (USE_MIC)

### CXX_FLAGS specific to compiler and platform #############
find_package (MPI REQUIRED)
message (STATUS "The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
message (STATUS "The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
message (STATUS "The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
message (STATUS "The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")

IF (LINUXGCC)
    MESSAGE (STATUS "Linux GCC flags SET")
    SET (IPPL_CXX_FLAGS
        "-DIPPL_LINUX -DIPPL_STRINGSTREAM ${IPPL_CXX_FLAGS}"
    )
    SET (OTHER_CXX_FLAGS
        "-ftemplate-depth-80 -funroll-loops -fstrict-aliasing ${OTHER_CXX_FLAGS}"
    )

    IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        SET (OTHER_CXX_FLAGS
             "-std=c++0x -diag-disable 383,981 ${OTHER_CXX_FLAGS}")

    ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        SET (OTHER_CXX_FLAGS
             "-Wall -Wsign-compare -Wunused-variable -Warray-bounds -DIPPL_RESTRICT_BUG ${OTHER_CXX_FLAGS}")

    ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        IF (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.8.0")
            SET (OTHER_CXX_FLAGS
	            "-Wall -std=c++0x")

        ELSEIF (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.7")
            SET (OTHER_CXX_FLAGS
                "-Wall -Wno-reorder -std=c++1y")
        ENDIF ()
    ELSE ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        # using unknown compiler
        MESSAGE (STATUS "${COMPILER_NAME} is using a compiler (${CMAKE_CXX_COMPILER_ID}) we were not thinking of!
                 Please use the gnu compiler or the intel compiler if you are having problems.")
    ENDIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")

    SET (IPPL_VERSION_NAME "\"V${IPPL_VERSION_MAJOR}.${IPPL_VERSION_MINOR} Linux GNU\"")
ENDIF (LINUXGCC)

IF (XT3GCC)
    MESSAGE (STATUS "XT3 GCC flags SET")
    SET (IPPL_CXX_FLAGS
        "-DIPPL_XT3 -DIPPL_DEBUG -DIPPL_LINUX -DIPPL_STRINGSTREAM ${IPPL_CXX_FLAGS}"
    )
    SET (OTHER_CXX_FLAGS
        "-Wno-deprecated -ftemplate-depth-80 -funroll-loops -fstrict-aliasing"
    )
    SET (IPPL_VERSION_NAME "\"V${IPPL_VERSION_MAJOR}.${IPPL_VERSION_MINOR} XT3 GNU\"")
ENDIF (XT3GCC)

IF (LINUXPGI)
    MESSAGE(STATUS "Linux PGI flags SET")
    SET (IPPL_CXX_FLAGS
        "-DIPPLProf -DIPPL_GPLUSPLUS -DMPICH_IGNORE_CXX_SEEK -DIPPL_LINUX ${IPPL_CXX_FLAGS}"
    )
    SET (IPPL_VERSION_NAME "\"V${IPPL_VERSION_MAJOR}.${IPPL_VERSION_MINOR} Linux PGI\"")
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
IF (ENABLE_IPPLTESTS)
    ADD_SUBDIRECTORY (test)
ENDIF (ENABLE_IPPLTESTS)

### write configure file ####################################
CONFIGURE_FILE ( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake )

### install files ###########################################
INSTALL (
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake
    DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/${PROJECT_NAME}"
    RENAME ${PROJECT_NAME}Config.cmake
)
