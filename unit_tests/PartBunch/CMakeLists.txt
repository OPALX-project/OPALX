file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestBCHandler USE_GTEST_MAIN)
add_opalx_test(TestBeamBeam USE_GTEST_MAIN)
add_opalx_test(TestBunchStateHandler USE_GTEST_MAIN)
add_opalx_test(TestMultiContainerPartBunch USE_GTEST_MAIN)
add_opalx_test(TestWitnessDeletionMPI2
  SOURCES TestMultiContainerPartBunch.cpp NUM_PROCS 2 TIMEOUT 60 USE_GTEST_MAIN
  ARGS --gtest_filter=MultiContainerPartBunchTest.WitnessDeletionPreservesIdentityAndPhaseSpaceInH5)
add_opalx_test(TestParticleContainer USE_GTEST_MAIN)

add_subdirectory(Binning)
