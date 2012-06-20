cmake_minimum_required (VERSION 2.6)

project (MC4Initializer)
set (MC4Initializer_VERSION_MAJOR 0)
set (MC4Initializer_VERSION_MINOR 1)

SET (CMAKE_CXX_COMPILER ${MPI_COMPILER})

IF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
    )
ENDIF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

MESSAGE (STATUS "Build type is: " ${CMAKE_BUILD_TYPE})

SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
SET (CMAKE_CXX_FLAGS_RELEASE "-O2")
SET (CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

SET (MC4Initializer_SRC
    Cosmology.cpp
    DataBase.cpp
    Initializer.cpp
    InputParser.cpp
    MT_Random.cpp
    Output.cpp
    Parallelization.cpp
    PerfMon.c
)

SET (CMAKE_CXX_FLAGS
    "-Wall -Werror -std=c++0x -DUSENAMESPACE ${CMAKE_CXX_FLAGS}")

INCLUDE_DIRECTORIES (
    ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_LIBRARY (MC4Initializer ${MC4Initializer_SRC})
