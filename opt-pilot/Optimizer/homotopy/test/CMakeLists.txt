include_directories(../src)
include_directories(../dependencies)
include_directories(../dependencies/lib)
#include_directories(../dependencies/gtest-1.6.0/include)
include_directories( ${gtest_SOURCE_DIR}/include )
include_directories(../dependencies/nlopt-2.3/api)
include_directories(../dependencies/alglib-3.6.0/cpp/src)

enable_testing()

add_executable(testSandbox
    test_sandbox.cpp
)
target_link_libraries(testSandbox gtest gtest_main nlopt Problems )
add_test(
    NAME testSandbox
    COMMAND testSandbox
)

add_executable(testNlopt
    test_nlopt.cpp
)
target_link_libraries(testNlopt gtest gtest_main nlopt )
add_test(
    NAME testNlopt
    COMMAND testNlopt
)

add_executable(testAlgLib
    test_alglib.cpp
)
target_link_libraries(testAlgLib gtest gtest_main alglib )
add_test(
    NAME testAlgLib
    COMMAND testAlgLib
)

add_executable(testInterpolator
    test_interpolator.cpp
)
target_link_libraries(testInterpolator gtest gtest_main alglib )
add_test(
    NAME testInterpolator
    COMMAND testInterpolator
)

file( COPY data/ DESTINATION data/)
add_executable(testHypervolume
    test_hypervolume.cpp
)
target_link_libraries(testHypervolume gtest gtest_main nlopt wfghypervolume)
add_test(
    NAME testHypervolume
    COMMAND testHypervolume
)

add_executable(testWFGPROBLEMS
    test_wfgproblems.cpp
)
target_link_libraries(testWFGPROBLEMS gtest gtest_main nlopt wfgProblems)
add_test(
    NAME testWFGPROBLEMS
    COMMAND testWFGPROBLEMS
)

add_executable(testProblems
    test_problems.cpp
)
target_link_libraries(testProblems gtest gtest_main nlopt Problems)
add_test(
    NAME testProblems
    COMMAND testProblems
)

add_executable(testMesh
    test_mesh.cpp
)
target_link_libraries(testMesh gtest gtest_main )
add_test(
    NAME testMesh
    COMMAND testMesh
)

add_executable(testMeshSet
    test_meshset.cpp
)
target_link_libraries(testMeshSet gtest gtest_main )
add_test(
    NAME testMeshSet
    COMMAND testMeshSet
)

add_executable(testEvaluator
    test_evaluator.cpp
)
target_link_libraries(testEvaluator gtest gtest_main Problems )
add_test(
    NAME testEvaluator
    COMMAND testEvaluator
)

add_executable(testOptimizer
    test_optimizer.cpp
)
target_link_libraries(testOptimizer gtest gtest_main nlopt Problems )
add_test(
    NAME testOptimizer
    COMMAND testOptimizer
)

add_executable(testMarching
    test_marching.cpp
)
target_link_libraries(testMarching gtest gtest_main nlopt Problems )
add_test(
    NAME testMarching
    COMMAND testMarching
)

add_executable(testJobQueue
    test_jobqueue.cpp
)
target_link_libraries(testJobQueue gtest gtest_main nlopt Problems )
add_test(
    NAME testJobQueue
    COMMAND testJobQueue
)

add_executable(testHomotopy
    test_homotopy.cpp
)
target_link_libraries(testHomotopy gtest gtest_main nlopt Problems )
add_test(
    NAME testHomotopy
    COMMAND testHomotopy
)

add_executable(testOptimizerDetail
    test_optimizer_detail.cpp
)
target_link_libraries(testOptimizerDetail gtest gtest_main nlopt Problems )
add_test(
    NAME testOptimizerDetail
    COMMAND testOptimizerDetail
)

add_executable(testBenchmarks
    test_benchmarks.cpp
)
target_link_libraries(testBenchmarks gtest gtest_main nlopt Problems )
add_test(
    NAME testBenchmarks
    COMMAND testBenchmarks
)
