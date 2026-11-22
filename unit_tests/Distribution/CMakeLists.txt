file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestGaussian USE_GTEST_MAIN)
add_opalx_test(TestMultiVariateGaussian USE_GTEST_MAIN)
add_opalx_test(TestFlatTop USE_GTEST_MAIN)
add_opalx_test(TestOpalFlatTop USE_GTEST_MAIN)
add_opalx_test(TestFromFile USE_GTEST_MAIN)
add_opalx_test(TestEmittedFromFile USE_GTEST_MAIN)
add_opalx_test(TestEmittedFromFileMPI2
  SOURCES TestEmittedFromFile.cpp NUM_PROCS 2 TIMEOUT 30 USE_GTEST_MAIN)
if(BUILD_TESTING AND TARGET TestEmittedFromFile)
  set_tests_properties(TestEmittedFromFile PROPERTIES
    TIMEOUT 30 RUN_SERIAL TRUE LABELS "unit;beambeam;regression"
    ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  set_tests_properties(TestEmittedFromFileMPI2 PROPERTIES
    RUN_SERIAL TRUE LABELS "unit;beambeam;regression;mpi"
    ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  # Exercise both legacy headers and explicit timed births on one and two ranks.
endif()
add_opalx_test(TestUniform USE_GTEST_MAIN)
