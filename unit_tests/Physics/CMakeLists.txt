file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestLinearCompton USE_GTEST_MAIN)

add_opalx_test(TestLinearComptonSpectrum USE_GTEST_MAIN NO_MPI)

add_executable(LinearComptonSpectrumBenchmark LinearComptonSpectrumBenchmark.cpp)
target_link_libraries(LinearComptonSpectrumBenchmark PRIVATE opalx)
target_include_directories(LinearComptonSpectrumBenchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(TestLinearComptonSpectrum PRIVATE
    OPALX_TEST_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
    OPALX_LINEAR_COMPTON_BENCHMARK="$<TARGET_FILE:LinearComptonSpectrumBenchmark>"
)
