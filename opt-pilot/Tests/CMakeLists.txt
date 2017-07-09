# gtest and python libraries are required for the tests
find_package(GTest)
find_package(PythonLibs)

if(NOT PYTHONLIBS_FOUND OR NOT GTEST_FOUND)
    if (NOT PYTHONLIBS_FOUND)
        message (STATUS "Python libraries not found.")
    endif()
    if (NOT GTEST_FOUND)
        message (STATUS "gtest not found.")
    endif()
    message (STATUS "No testing possible.")
    return()
endif()

include_directories (
    $ENV{GTEST_PREFIX}/include
    ${Boost_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
    ${PYTHON_INCLUDE_DIRS}
)

link_directories (
    $ENV{GTEST_PREFIX}/lib
    ${PYTHON_LIBRARIES}
)

set (OPT_PILOT_SOURCE_DIR ${CMAKE_SOURCE_DIR})
IF (OPTP_USED_FROM_OPAL)
    set (OPT_PILOT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/opt-pilot)
ENDIF()

set (ExpressionTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSReader.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/expression.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/evaluator.cpp
)

set (SumErrSqExpressionTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Expression/SumErrSq.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSReader.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/expression.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/evaluator.cpp
)

set (FromFileExpressionTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Expression/FromFile.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/expression.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/evaluator.cpp
)

set (PythonExprTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/expression.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/evaluator.cpp
)


set (HashNameGeneratorTest_SRC
)

set (CmdArgumentsTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Util/CmdArguments.cpp
)

set (MPIHelperTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Util/MPIHelper.cpp
)

set (ManagedIDsTest_SRC
)

set (IndividualTest_SRC
)

set (PopulationTest_SRC
)

set (SDDSParserTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSReader.cpp
)

set (OpalInputFileParserTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Expression/SumErrSq.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/FromFile.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSReader.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/CommentAnnotatedInputFileParser.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/expression.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/evaluator.cpp
)

set (PILOT_TESTS
    MPIHelperTest
    OpalInputFileParserTest
    SDDSParserTest
    ManagedIDsTest
    IndividualTest
    PopulationTest
    ExpressionTest
    SumErrSqExpressionTest
    FromFileExpressionTest
    PythonExprTest
    CmdArgumentsTest
    HashNameGeneratorTest
)

set (TEST_LIBS
    gtest
    pthread
    boost_system
    boost_serialization
    mpi
    mpi_cxx
)

set (PYTHON_LIBS
    ${PYTHON_LIBRARIES}
    boost_python
)

FOREACH(test ${PILOT_TESTS})
    SET( SRC ${${test}_SRC} )
    ADD_EXECUTABLE( ${test}.exe ${test}.cpp ${SRC} )
    TARGET_LINK_LIBRARIES( ${test}.exe ${TEST_LIBS} ${PYTHON_LIBS})
    #target_link_libraries(foo ${GTEST_BOTH_LIBRARIES})
    ADD_TEST(${test} ${test}.exe)
ENDFOREACH(test ${PILOT_TESTS})

#HACK: link resources dir to test dir
if (NOT EXISTS ${PROJECT_BINARY_DIR}/Tests/resources)
    EXECUTE_PROCESS(
        COMMAND ln -s ${OPT_PILOT_SOURCE_DIR}/Tests/resources ${PROJECT_BINARY_DIR}/Tests/
        OUTPUT_VARIABLE LN_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
        )
endif()