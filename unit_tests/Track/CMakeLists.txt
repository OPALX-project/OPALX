file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestDecayChain USE_GTEST_MAIN)

add_opalx_test(TestMatchCmd USE_GTEST_MAIN LABELS tps)
add_opalx_test(TestMatchCmdMPI SOURCES TestMatchCmd.cpp NUM_PROCS 2 USE_GTEST_MAIN LABELS tps)
