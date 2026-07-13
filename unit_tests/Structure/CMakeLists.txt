file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestBinningCmd USE_GTEST_MAIN)
add_opalx_test(TestBinConfigWriter USE_GTEST_MAIN)
add_opalx_test(TestDirichletPlaneWriter USE_GTEST_MAIN)
add_opalx_test(TestLaser USE_GTEST_MAIN)
add_opalx_test(TestBeam USE_GTEST_MAIN)
add_opalx_test(TestLossDataSink USE_GTEST_MAIN)
