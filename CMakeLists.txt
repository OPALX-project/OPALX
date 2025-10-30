cmake_minimum_required(VERSION 3.20)

# --------------------------------------------------
# Global policy to avoid old subproject issues
# --------------------------------------------------
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE INTERNAL "")

project(OPALX VERSION 2024.1.0 LANGUAGES C CXX)

include(FetchContent)

# --------------------------------------------------
# General settings
# --------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

message(STATUS "=== Configuring OPALX ${PROJECT_VERSION} ===")

# --------------------------------------------------
# Helper: fetch dependencies if missing
# --------------------------------------------------
function(fetch_dep pkg_name repo branch)
    FetchContent_Declare(${pkg_name}
        GIT_REPOSITORY ${repo}
        GIT_TAG ${branch}
        GIT_SHALLOW TRUE
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
endfunction()

# --------------------------------------------------
# Declare dependencies
# --------------------------------------------------
fetch_dep(ippl https://github.com/IPPL-framework/ippl.git master)
fetch_dep(gsl https://github.com/ampl/gsl.git master)          # use master branch
fetch_dep(boost https://github.com/boostorg/boost.git boost-1.82.0)
fetch_dep(hdf5 https://github.com/HDFGroup/hdf5.git hdf5-1_10_8)

option(BUILD_OPAL_UNIT_TESTS "Build unit tests" OFF)
if(BUILD_OPAL_UNIT_TESTS)
    fetch_dep(gtest https://github.com/google/googletest.git v1.13.0)
endif()

# --------------------------------------------------
# Custom target: dependencies
# --------------------------------------------------
add_custom_target(dependencies
    COMMENT "Fetching and building all external dependencies..."
)

# --------------------------------------------------
# Helper macro to populate and add subdirectories
# --------------------------------------------------
macro(populate_dep dep)
    FetchContent_GetProperties(${dep})
    if(NOT ${dep}_POPULATED)
        FetchContent_Populate(${dep})
    endif()

    # Patch old cmake_minimum_required in subproject (like GSL)
    set(cmake_file "${${dep}_SOURCE_DIR}/CMakeLists.txt")
    if(EXISTS "${cmake_file}")
        file(READ "${cmake_file}" contents)
        string(REGEX REPLACE "cmake_minimum_required\\(VERSION [0-9.]+\\)"
                             "cmake_minimum_required(VERSION 3.20)" 
                             contents "${contents}")
        file(WRITE "${cmake_file}" "${contents}")
    endif()

    # Add subdirectory if it exists
    if(EXISTS "${${dep}_SOURCE_DIR}/CMakeLists.txt")
        add_subdirectory(${${dep}_SOURCE_DIR} ${${dep}_BINARY_DIR})
    endif()

    # Only add dependency if target exists
    if(TARGET ${dep})
        add_dependencies(dependencies ${dep})
    endif()
endmacro()

# --------------------------------------------------
# Populate all dependencies
# --------------------------------------------------
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(Boost_USE_ICU OFF)

foreach(dep ippl gsl boost hdf5)
    populate_dep(${dep})
endforeach()

if(BUILD_OPAL_UNIT_TESTS)
    populate_dep(gtest)
endif()

# --------------------------------------------------
# Add OPALX sources
# --------------------------------------------------
add_subdirectory(src)
include_directories(${CMAKE_SOURCE_DIR}/src)

# --------------------------------------------------
# Make OPALX depend on all dependencies
# --------------------------------------------------
if(TARGET OPALX)
    target_link_libraries(OPALX
        PRIVATE
            $<$<TARGET_EXISTS:ippl>:ippl>
            $<$<TARGET_EXISTS:gsl>:gsl>
            $<$<TARGET_EXISTS:boost>:boost>
            $<$<TARGET_EXISTS:hdf5>:hdf5>
            $<$<TARGET_EXISTS:gtest>:gtest>
    )
    add_dependencies(OPALX dependencies)
endif()

message(STATUS "✅ Setup complete. Run 'make dependencies' first, then 'make' to build OPALX.")
