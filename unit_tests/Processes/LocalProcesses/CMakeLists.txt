file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestLocalProcessAttachment USE_GTEST_MAIN)
add_opalx_test(TestTestCounterProcess USE_GTEST_MAIN)
add_opalx_test(TestOpalLocalProcess USE_GTEST_MAIN)
