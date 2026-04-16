file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestBCHandler USE_GTEST_MAIN)
add_opalx_test(TestBinnedFieldSolver USE_GTEST_MAIN)
add_opalx_test(TestMultiContainerPartBunch USE_GTEST_MAIN)
add_opalx_test(TestParticleContainer USE_GTEST_MAIN)

add_subdirectory(Binning)

