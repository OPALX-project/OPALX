# Preserve the numerical regression executables and their process-isolated cache tests.
function(add_tps_regression name)
  add_executable(tps_${name} ${name}.cpp ${ARGN})
  target_link_libraries(tps_${name} PRIVATE TPS::tps)
  add_test(NAME tps.${name} COMMAND tps_${name})
  set_tests_properties(tps.${name} PROPERTIES LABELS "unit;tps" TIMEOUT 60)
endfunction()

foreach(name numerical truncation cache_bounds coefficient_indexing array_indexing
             column_iterators modern_storage workflow_algebra map_defects map_context
             kokkos_linear_algebra runtime_lifecycle public_api)
  add_tps_regression(${name})
endforeach()
add_tps_regression(prepared_map device_map_header.cpp)
foreach(regression constant_multiplication linear_substitution)
  add_test(NAME tps.${regression} COMMAND tps_numerical ${regression})
  set_tests_properties(tps.${regression} PROPERTIES LABELS "unit;tps" TIMEOUT 60)
endforeach()
add_executable(tps_zero_order zero_order.cpp)
target_link_libraries(tps_zero_order PRIVATE TPS::tps)
foreach(variables 1 2 6)
  foreach(mode monomial series)
    if(mode STREQUAL series)
      set(args series_first)
    else()
      set(args)
    endif()
    add_test(NAME tps.zero_order_${mode}_${variables}
             COMMAND tps_zero_order ${variables} ${args})
    set_tests_properties(tps.zero_order_${mode}_${variables}
                         PROPERTIES LABELS "unit;tps" TIMEOUT 60)
  endforeach()
endforeach()
foreach(example maps context_maps prepared_evaluation)
  add_executable(tps_example_${example} examples/${example}.cpp)
  target_link_libraries(tps_example_${example} PRIVATE TPS::tps)
  add_test(NAME tps.example_${example} COMMAND tps_example_${example})
  set_tests_properties(tps.example_${example} PROPERTIES LABELS "unit;tps" TIMEOUT 60)
endforeach()
foreach(order forward reverse)
  add_executable(tps_opalx_${order} opalx_headers.cpp)
  target_link_libraries(tps_opalx_${order} PRIVATE OPALX::opalx)
  if(order STREQUAL reverse)
    target_compile_definitions(tps_opalx_${order} PRIVATE TPS_REVERSE_INCLUDE_ORDER)
  endif()
  add_test(NAME tps.opalx_headers_${order} COMMAND tps_opalx_${order})
  set_tests_properties(tps.opalx_headers_${order} PROPERTIES LABELS "unit;tps" TIMEOUT 60)
endforeach()

# Check that no public header relies on inclusion through tps.hpp.
file(GLOB_RECURSE tps_headers CONFIGURE_DEPENDS
     RELATIVE "${PROJECT_SOURCE_DIR}/src" "${PROJECT_SOURCE_DIR}/src/tps/*.h"
     "${PROJECT_SOURCE_DIR}/src/tps/*.hpp")
set(header_checks)
foreach(header IN LISTS tps_headers)
  string(MAKE_C_IDENTIFIER "${header}" header_id)
  set(check "${CMAKE_CURRENT_BINARY_DIR}/header-checks/${header_id}.cpp")
  file(WRITE "${check}" "#include <${header}>\n")
  list(APPEND header_checks "${check}")
endforeach()
add_library(tps_header_checks OBJECT ${header_checks})
target_link_libraries(tps_header_checks PRIVATE TPS::tps)

# Use OPALX's real IPPL/MPI lifecycle and its normal unit-test registration helper.
add_opalx_test(TestTPSIntegration USE_GTEST_MAIN LABELS tps NUM_PROCS 1)
add_opalx_test(TestTPSIntegrationMPI SOURCES TestTPSIntegration.cpp
               USE_GTEST_MAIN LABELS tps NUM_PROCS 2)
