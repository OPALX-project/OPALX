# -----------------------------------------------------------------------------
# CompilerOptions.cmake
# cmake-format: off
#
# Sets compiler flags that affect how all OPALX targets are built.
#
# Responsibilities:
#   - Warning flags (-Wall, -Wextra, etc.)
#   - Debug sanitizers (ASan/UBSan) when using Debug build type
#   - Compiler-specific warning suppressions
#
# Not responsible for:
#   - Enabling CUDA/OpenMP/Serial                 → Platforms.cmake
#   - Selecting platform specific compiler flags  → Platforms.cmake 
#
# This file is only concerned with general correctness and development-time safety.
#
# cmake-format: on
# -----------------------------------------------------------------------------

# === Basic warnings (apply to OPALX's C and C++ languages) ===
# OPALX enables Fortran only to build the fetched reference LAPACK. Do not pass
# the project's warning policy into that third-party implementation.
add_compile_options(
  $<$<COMPILE_LANGUAGE:C,CXX>:-Wall>
  $<$<COMPILE_LANGUAGE:C,CXX>:-Wextra>
  $<$<COMPILE_LANGUAGE:C,CXX>:-Wno-deprecated-declarations>)

# === Use modified variant implementation ===
if(OPALX_USE_ALTERNATIVE_VARIANT)
  add_definitions(-DOPALX_USE_ALTERNATIVE_VARIANT)
endif()

# === Code coverage options ===
if(OPALX_ENABLE_COVERAGE AND (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang"))
  message(STATUS "${ColorYellow}Code coverage enabled.${ColorReset}")
  add_compile_options(-fprofile-arcs -ftest-coverage -g)
  add_link_options(-fprofile-arcs -ftest-coverage)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
endif()

# === Compiler-specific warning suppressions ===
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  add_compile_options(
    $<$<NOT:$<COMPILE_LANGUAGE:Fortran>>:-Wno-deprecated-copy>
    $<$<NOT:$<COMPILE_LANGUAGE:Fortran>>:-Wno-sign-compare>)
endif()

# GCC 12+ false positives for buffer overflows, restrict, etc.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12)
  add_compile_options(
    $<$<NOT:$<COMPILE_LANGUAGE:Fortran>>:-Wno-stringop-overflow>
    $<$<NOT:$<COMPILE_LANGUAGE:Fortran>>:-Wno-array-bounds>
    $<$<NOT:$<COMPILE_LANGUAGE:Fortran>>:-Wno-restrict>)
endif()

# === Debug-specific sanitizers ===
if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU"
   AND OPALX_ENABLE_SANITIZER)
  message(STATUS "✅ Enabling AddressSanitizer and UBSan for Debug build")
  #  add_compile_options(-fsanitize=address)
  # add_link_options(-fsanitize=address)
endif()

# === Position Independent Code (PIC) for shared libraries ===
if(BUILD_SHARED_LIBS)
  message(STATUS "✅ Enabling Position Independent Code (PIC) for shared libraries")
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  set(CMAKE_LINK_DEPENDS_NO_SHARED true)
endif()

message(STATUS "✅ Compiler options configured")
