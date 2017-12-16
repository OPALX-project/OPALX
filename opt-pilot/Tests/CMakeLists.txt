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

IF (OPTP_USED_FROM_OPAL)
    STRING (REGEX REPLACE "/[^/]+$" "/" GTEST_LIBRARY_DIRS ${GTEST_LIBRARY})
    SET (RELATIVE_TO_OPAL "opt-pilot/")
ENDIF()

include_directories (
    ${GTEST_INCLUDE_DIRS}
    $ENV{GTEST_PREFIX}/include
    ${Boost_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
    ${PYTHON_INCLUDE_DIRS}
)

link_directories (
    $ENV{GTEST_PREFIX}/lib
    ${GTEST_LIBRARY_DIRS}
    ${PYTHON_LIBRARIES}
)

set (OPT_PILOT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/opt-pilot)
set (SDDSPARSER_SRC
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/array.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/associate.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/ast.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/column.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/data.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/description.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/file.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/include.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/parameter.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/version.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser.cpp
)

set (ExpressionTest_SRC
    ${SDDSPARSER_SRC}
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/expression.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/evaluator.cpp
)

set (SumErrSqExpressionTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Expression/SumErrSq.cpp
    ${SDDSPARSER_SRC}
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
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/expression.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/evaluator.cpp
)

set (PopulationTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/expression.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/evaluator.cpp
)

set (SDDSParserTest_SRC
    ${SDDSPARSER_SRC}
)

set (OpalInputFileParserTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Expression/SumErrSq.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/FromFile.cpp
    ${SDDSPARSER_SRC}
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
#    PythonExprTest
    CmdArgumentsTest
    HashNameGeneratorTest
)

set (TEST_LIBS
    gtest
    pthread
    ${Boost_LIBRARIES}
    ${MPI_LIBRARIES}
)

set (PYTHON_LIBS
#    ${PYTHON_LIBRARIES}
#    boost_python
)

FOREACH(test ${PILOT_TESTS})
    SET( SRC ${${test}_SRC} )
    ADD_EXECUTABLE( ${test}.exe ${test}.cpp ${SRC} )
    TARGET_LINK_LIBRARIES( ${test}.exe ${TEST_LIBS} ${PYTHON_LIBS})
    #target_link_libraries(foo ${GTEST_BOTH_LIBRARIES})
    ADD_TEST(${test} ${test}.exe)
ENDFOREACH(test ${PILOT_TESTS})

#HACK: link resources dir to test dir
if (NOT EXISTS ${PROJECT_BINARY_DIR}/${RELATIVE_TO_OPAL}Tests/resources)
    EXECUTE_PROCESS(
        COMMAND ln -s ${OPT_PILOT_SOURCE_DIR}/Tests/resources ${PROJECT_BINARY_DIR}/${RELATIVE_TO_OPAL}Tests/
        OUTPUT_VARIABLE LN_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
        )
endif()