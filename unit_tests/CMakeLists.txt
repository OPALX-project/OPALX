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

add_subdirectory(OpalParser)

add_subdirectory(AbsBeamline)
add_subdirectory(BeamlineGeometry)
add_subdirectory(Distribution)
add_subdirectory(Fields)
add_subdirectory(Physics)
add_subdirectory(Utilities)
add_subdirectory(Structure)
add_subdirectory(Algorithms)
add_subdirectory(PartBunch)
add_subdirectory(SpaceCharge)
add_subdirectory(BasicActions)
add_subdirectory(Elements)
add_subdirectory(Moments)
add_subdirectory(Processes)
add_subdirectory(Steppers)
add_subdirectory(Track)

# Run-side sources read no OPTION globals and no OpalData and include no input-layer header.
if(BUILD_TESTING)
  add_test(NAME TestRunSideBoundary
           COMMAND ${CMAKE_COMMAND} -DOPALX_SOURCE_DIR=${PROJECT_SOURCE_DIR}/src -P
                   ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CheckRunSideBoundary.cmake)
  set_tests_properties(TestRunSideBoundary PROPERTIES LABELS unit TIMEOUT 120)
endif()

# Links every object of opalx_run without opalx, so a run-side reference to an input-side or
# bridge symbol is an undefined symbol at link time. A static opalx_run is linked whole, since the
# linker otherwise takes only the archive members it needs. A shared opalx_run holds all its
# objects anyway: Apple ld rejects its undefined symbols when it links the library, GNU ld when it
# links this executable.
add_executable(TestRunSideLinksAlone TestRunSideLinksAlone.cpp)
get_target_property(_opalx_run_type opalx_run TYPE)
if(_opalx_run_type STREQUAL "STATIC_LIBRARY")
  target_link_libraries(TestRunSideLinksAlone PRIVATE "$<LINK_LIBRARY:WHOLE_ARCHIVE,opalx_run>")
else()
  target_link_libraries(TestRunSideLinksAlone PRIVATE opalx_run)
endif()
if(BUILD_TESTING)
  add_test(NAME TestRunSideLinksAlone COMMAND TestRunSideLinksAlone)
  set_tests_properties(TestRunSideLinksAlone PROPERTIES LABELS unit TIMEOUT 60)
endif()
