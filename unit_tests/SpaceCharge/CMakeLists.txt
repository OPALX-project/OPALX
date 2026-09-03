file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestSelfFieldConfig USE_GTEST_MAIN)
add_opalx_test(TestPic2d5Solver USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestSelfFieldSystem USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestPic3DPolicies USE_GTEST_MAIN)
add_opalx_test(TestBinningPlan USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestCartesianPicInfrastructure USE_GTEST_MAIN RUN_SERIAL)

set(_pic2d5_regression_dir "${CMAKE_CURRENT_BINARY_DIR}/Pic2d5EndToEnd")
file(MAKE_DIRECTORY "${_pic2d5_regression_dir}")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/data/Pic2d5EndToEnd.in"
  "${_pic2d5_regression_dir}/Pic2d5EndToEnd.in" COPYONLY)

if(BUILD_TESTING)
  add_test(
    NAME Pic2d5EndToEnd
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:opalx_exe> Pic2d5EndToEnd.in --info 2)
  set_tests_properties(
    Pic2d5EndToEnd
    PROPERTIES WORKING_DIRECTORY "${_pic2d5_regression_dir}"
               LABELS "unit;regression;space-charge;fft2d5"
               RUN_SERIAL TRUE
               TIMEOUT 120
               PASS_REGULAR_EXPRESSION "Total FieldSolver calls: 16")
endif()
