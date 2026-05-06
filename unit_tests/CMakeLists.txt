# -----------------------------------------------------------------------------
# unit_test/CMakeLists.txt
#
# Configures and adds unit tests. This file assumes OPALX_ENABLE_UNIT_TESTS is ON.
# -----------------------------------------------------------------------------

if(NOT TARGET GTest::gtest)
  message(FATAL_ERROR "❌ GoogleTest not found. Check Dependencies.cmake.")
endif()

include(GoogleTest)

if(NOT TARGET opalx_test_support)
  add_library(opalx_test_support INTERFACE)
  target_include_directories(opalx_test_support
                             INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
  add_library(opalx::test_support ALIAS opalx_test_support)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/cmake/AddOPALXTest.cmake")

add_subdirectory(AbsBeamline)
add_subdirectory(BeamlineCore)
add_subdirectory(BeamlineGeometry)
add_subdirectory(Distribution)
add_subdirectory(Fields)
add_subdirectory(Physics)
add_subdirectory(Utilities)
add_subdirectory(Structure)
add_subdirectory(Algorithms)
add_subdirectory(PartBunch)
add_subdirectory(BasicActions)
add_subdirectory(Elements)
add_subdirectory(Processes)
