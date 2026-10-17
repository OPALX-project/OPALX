find_package(Python3 COMPONENTS Interpreter QUIET)
if(BUILD_TESTING AND Python3_Interpreter_FOUND)
  add_test(NAME ParseOnlyCLI
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_parse_only.py"
            $<TARGET_FILE:opalx_exe>)
  set_tests_properties(ParseOnlyCLI PROPERTIES TIMEOUT 180 LABELS "integration;parser")
endif()
