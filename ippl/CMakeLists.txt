IF (NOT IPPL_USED_FROM_OPAL)
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.10)
PROJECT (IPPL)
SET (IPPL_VERSION_MAJOR 1)
SET (IPPL_VERSION_MINOR 5.00.1)

IF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
	)
ENDIF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

MESSAGE (STATUS "Build type is: " ${CMAKE_BUILD_TYPE})

# Select flags.
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ")
SET (CMAKE_CXX_FLAGS_RELEASE "-O3")
SET (CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

# Resolve all library dependencies
SET (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

find_package (MPI REQUIRED)
message (STATUS "The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
message (STATUS "The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
message (STATUS "The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
message (STATUS "The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")

# compiler dependent flags
if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -diag-disable 383 -diag-disable 981 -std=c++11")

elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wsign-compare -Wunused-variable -Warray-bounds -DIPPL_RESTRICT_BUG")

elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    if (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.7.3")
        message (FATAL_ERROR "To build OPAL you need GCC version 4.7.3 or greater")
    elseif (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.7.2")
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-unused-local-typedefs -std=c++11")
    endif ()
else ()
    message (STATUS "Unknown C++ compiler. Please use the GNU or Intel compiler, IF you are having problems.")
endif ()

FIND_PACKAGE (H5Hut REQUIRED)
FIND_PACKAGE (HDF5 REQUIRED)
FIND_PACKAGE (GSL REQUIRED)

set (BOOSTROOT $ENV{BOOST_DIR})
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package (Boost 1.55.0 REQUIRED COMPONENTS  regex filesystem system)
if (Boost_INCLUDE_DIRS)
    message (STATUS "Found boost include dir: ${Boost_INCLUDE_DIR}")
    message (STATUS "Found boost library dir: ${Boost_LIBRARY_DIR}")
    message (STATUS "Found boost libraries: ${Boost_LIBRARIES}")
    include_directories (${Boost_INCLUDE_DIRS})
endif ()


### Uldis Locans 2014.09.05 ###
#OPTION (ENABLE_DKS "Enable DKS" OFF)

### Load DKS package ###
IF (ENABLE_DKS)
   #FIND_PACKAGE(DKS 1.1.0 REQUIRED HINTS $ENV{DKS_PREFIX} $ENV{DKS_DIR} $ENV{DKS})
   #MESSAGE (STATUS "Found DKS library: ${DKS_LIBRARY}")
   #MESSAGE (STATUS "Found DKS include dir: ${DKS_INCLUDE_DIR}")
   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDKS_MPI")

   ### OpenCL compiler flags ###
   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lOpenCL -pthread -DDKS_OPENCL")

   ### CUDA compiler flags ###
   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lcudart -lcufft -lcublas -lnvToolsExt -DDKS_CUDA")

   ### if any accelerator enabled set flag to use DKS ###
   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPPL_DKS -DIPPL_DKS_CUDA -DOPAL_DKS")
ENDIF (ENABLE_DKS)
### END DKS part ###

# Handle options
OPTION (ENABLE_AMR_SOLVER "Enable BoxLib based AMR solver" OFF)

OPTION (ENABLE_SAAMG_SOLVER "Enable iteartive SA-AMG-PCG self field solver" OFF)

OPTION (DBG_SCALARFIELD "Enable dump of scalar field rho_m" OFF)
OPTION (DBG_STENCIL "Dump stencil of iterative solver in matlab format" OFF)
OPTION (BUILD_LIBOPAL "Enable building the OPAL library" OFF)

OPTION (NO_FIELD_ASSIGN_OPTIMIZATION "Disable compiler optimization of IPPL field assignment" OFF)

IF (ENABLE_AMR_SOLVER)
    # set specific parameters for BoxLib used in OPAL
    ENABLE_LANGUAGE (Fortran)

    SET ( CCSE_DIR $ENV{BOXLIB_PREFIX} )
    SET ( CCSE_LIBRARY_DIR $ENV{BOXLIB_LIBRARY_DIR} )

    FIND_PACKAGE (CCSE REQUIRED)

    IF (CCSE_FOUND)
        MESSAGE (STATUS "Found BoxLib include dir: ${CCSE_INCLUDE_DIR}")
        MESSAGE (STATUS "Found BoxLib library dir: ${CCSE_LIBRARY_DIR}")
        MESSAGE (STATUS "Found BoxLib perl dir: ${CCSE_PERL_DIR}")
        INCLUDE_DIRECTORIES (${CCSE_INCLUDE_DIR})
    ENDIF (CCSE_FOUND)

    SET (BL_SPACEDIM 3 CACHE INT "Dimension of BoxLib build")
    SET (ENABLE_MPI 1 CACHE INT "Enable build with MPI")
    SET (ENABLE_OpenMP 0 CACHE INT "Enable build with OpenMP")
    SET (BL_PRECISION "DOUBLE" CACHE INT "Precision of BoxLib build")
    SET (BL_USE_PARTICLES 1 CACHE INT "Include Particles classes in BoxLib build")
    SET (ENABLE_PROFILING 0 CACHE INT "Include profiling information in BoxLib build")
    SET (ENABLE_BACKTRACE 0 CACHE INT "Include backtrace information in BoxLib build")

    # disable due to BoxLib
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable")

    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBL_PRECISION=${BL_PRECISION}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_MPI=${ENABLE_MPI}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_OpenMP=${ENABLE_OpenMP}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_PARTICLES=${BL_USE_PARTICLES}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBL_DEBUG=${BL_DEBUG}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_BACKTRACE=${ENABLE_BACKTRACE}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_PROFILING=${ENABLE_PROFILING}")

    MESSAGE (STATUS "Use following BoxLib settings:")
    MESSAGE (STATUS "  BL_SPACEDIM = ${BL_SPACEDIM} (INT: 1, 2, 3)")
    MESSAGE (STATUS "  BL_MACHINE = ${BL_MACHINE} (STRING: <ARCH>)")
    MESSAGE (STATUS "  BL_PRECISION = ${BL_PRECISION} (STRING: \"FLOAT\", \"DOUBLE\")")
    MESSAGE (STATUS "  ENABLE_MPI = ${ENABLE_MPI} (INT: 0, 1)")
    MESSAGE (STATUS "  ENABLE_OpenMP = ${ENABLE_OpenMP} (INT: 0, 1)")
    MESSAGE (STATUS "  BL_DEBUG = ${BL_DEBUG} (INT: 0, 1)")
    MESSAGE (STATUS "  BL_USE_PARTICLES = ${BL_USE_PARTICLES} (INT: 0, 1)")
    MESSAGE (STATUS "  ENABLE_BACKTRACE = ${ENABLE_BACKTRACE} (INT: 0, 1)")
    MESSAGE (STATUS "  ENABLE_PROFILING = ${ENABLE_PROFILING} (INT: 0, 1)")
ENDIF (ENABLE_AMR_SOLVER)

if (ENABLE_SAAMG_SOLVER OR ENABLE_AMR_SOLVER)
    find_package (Trilinos REQUIRED HINTS $ENV{TRILINOS_PREFIX} $ENV{TRILINOS_DIR} $ENV{TRILINOS})

    message (STATUS "Found Trilinos: ${Trilinos_DIR}")
    message (STATUS "      Trilinos version: ${Trilinos_VERSION}")
    message (STATUS "      Trilinos package list: ${Trilinos_PACKAGE_LIST}")
    message (STATUS "      Trilinos libraries: ${Trilinos_LIBRARIES}")
    message (STATUS "      Trilinos TPL libraries: ${Trilinos_TPL_LIBRARIES}")

    if (NOT ${MPI_CXX_COMPILER} STREQUAL ${Trilinos_CXX_COMPILER} )
        message (STATUS "Compiler mismatch:")
        message (STATUS "        Trilinos was compiled with:")
        message (STATUS "            ${Trilinos_C_COMPILER}")
        message (STATUS "            ${Trilinos_CXX_COMPILER}")
        message (STATUS "            ${Trilinos_Fortran_COMPILER}")
        message (STATUS "        You are using:")
        message (STATUS "            ${CMAKE_C_COMPILER}")
        message (STATUS "            ${CMAKE_CXX_COMPILER}")
        message (STATUS "            ${CMAKE_Fortran_COMPILER}")
        message (FATAL  "")
    endif ()
endif ()

IF (DBG_SCALARFIELD)
    MESSAGE (STATUS "\nWrite scalar rho_m field is enabled ")
    SET (CMAKE_CXX_FLAGS "-DDBG_SCALARFIELD ${CMAKE_CXX_FLAGS}")
ENDIF (DBG_SCALARFIELD)

IF (DBG_STENCIL)
    MESSAGE (STATUS "\nWrite stencil is enabled ")
    SET (CMAKE_CXX_FLAGS "-DDBG_STENCIL ${CMAKE_CXX_FLAGS}")
ENDIF (DBG_STENCIL)

IF (BUILD_LIBOPAL)
    MESSAGE (STATUS "\nBuilding OPAL library")
ENDIF (BUILD_LIBOPAL)

IF (NO_FIELD_ASSIGN_OPTIMIZATION)
    MESSAGE (STATUS "\nCompiler will not try to optimize field assignment")
    SET (CMAKE_CXX_FLAGS
        "-DdontOPTIMIZE_FIELD_ASSIGNMENT ${CMAKE_CXX_FLAGS}"
	)
ENDIF (NO_FIELD_ASSIGN_OPTIMIZATION)

INCLUDE(FindGit)
IF(GIT_FOUND)
    SET (OPAL_WC_REVISION 0)
    EXECUTE_PROCESS(
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMAND git rev-parse HEAD
        OUTPUT_VARIABLE GIT_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    SET (OPAL_WC_REVISION \"${GIT_OUT}\")
ELSE()
    MESSAGE (STATUS "could not find git")
    SET (OPAL_WC_REVISION \"\")
ENDIF()

SET (IPPL_USED_FROM_OPAL "TRUE")
ADD_SUBDIRECTORY (ippl)
#ADD_SUBDIRECTORY (classic/5.0/src)
#ADD_SUBDIRECTORY (src)

IF (BUILD_OPAL_UNIT_TESTS)
     ADD_SUBDIRECTORY(tests)
ENDIF (BUILD_OPAL_UNIT_TESTS)

ELSE(NOT IPPL_USED_FROM_OPAL)
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
    "-DIPPL_MPI -DMPICH_SKIP_MPICXX -DIPPL_DONT_POOL -DIPPL_USE_XDIV_RNG -DPETE_BITWISE_COPY -DIPPL_HAS_TEMPLATED_COMPLEX -DIPPL_USE_PARTIAL_SPECIALIZATION -Drestrict=__restrict__ -DNOCTAssert ${IPPL_CXX_FLAGS}"
    )

### CXX_FLAGS common to all compilers and platforms which ###
### may *not* be written to the configure file            ###
SET (OTHER_CXX_FLAGS
    "-Wall ${OTHER_CXX_FLAGS}"
    )

### Uldis Locans 2014.09.05 ###
#OPTION (ENABLE_DKS "Enable DKS" OFF)

### Load DKS package ###
IF (ENABLE_DKS)
    #FIND_PACKAGE(DKS 1.1.0 REQUIRED HINTS $ENV{DKS_PREFIX} $ENV{DKS_DIR} $ENV{DKS})
    #MESSAGE (STATUS "Found DKS library: ${DKS_LIBRARY_DIR}")
    #MESSAGE (STATUS "Found DKS include dir: ${DKS_INCLUDE_DIR}")
    SET (IPPL_CXX_FLAGS "${IPPL_CXX_FLAGS} -DDKS_MPI")

    ### OpenCL compiler flags ###
    SET (IPPL_CXX_FLAGS "${IPPL_CXX_FLAGS} -lOpenCL -lpthread -DDKS_OPENCL")
    ### End OpenCL ###

    ### CUDA compiler flags ###
    SET (IPPL_CXX_FLAGS "${IPPL_CXX_FLAGS} -lcudart -lcufft -lcublas -lnvToolsExt -DDKS_CUDA")
    ### End Cuda ###

    ### if any accelerator enabled set flag to use DKS ###
    SET (IPPL_CXX_FLAGS "${IPPL_CXX_FLAGS} -DIPPL_DKS")
ENDIF (ENABLE_DKS)

### CXX_FLAGS specific to compiler and platform #############
#find_package (MPI REQUIRED)
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

    IF (NOT IPPL_USED_FROM_OPAL)
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
     ENDIF (NOT IPPL_USED_FROM_OPAL)
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

 SET (IPPL_CMAKE_CXX_FLAGS "${IPPL_CXX_FLAGS}" PARENT_SCOPE)
 SET (IPPL_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src" PARENT_SCOPE)
 SET (IPPL_LIBRARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src" PARENT_SCOPE)
 SET (IPPL_LIBRARY "ippl" PARENT_SCOPE)

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

ENDIF (NOT IPPL_USED_FROM_OPAL)