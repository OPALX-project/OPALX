# Use the same registration as the other OPALX tests, including custom-main tests.
# MPIEXEC also supplies the CSCS srun/uenv environment needed by the CUDA runtime.
function(add_tps_regression name)
  add_opalx_test(tps_${name} SOURCES ${name}.cpp ${ARGN} NUM_PROCS 1 TIMEOUT 60
                 LABELS tps CTEST_NAME tps.${name})
endfunction()

foreach(name numerical truncation cache_bounds coefficient_indexing array_indexing
             column_iterators modern_storage workflow_algebra map_defects map_context
             kokkos_linear_algebra linear_solve runtime_lifecycle public_api)
  add_tps_regression(${name})
endforeach()
add_tps_regression(prepared_map device_map_header.cpp)
foreach(regression constant_multiplication linear_substitution)
  add_opalx_test(tps_${regression} SOURCES numerical.cpp ARGS ${regression}
                 NUM_PROCS 1 TIMEOUT 60 LABELS tps CTEST_NAME tps.${regression})
endforeach()
foreach(variables 1 2 6)
  foreach(mode monomial series)
    if(mode STREQUAL series)
      set(args series_first)
    else()
      set(args)
    endif()
    add_opalx_test(tps_zero_order_${mode}_${variables} SOURCES zero_order.cpp
                   ARGS ${variables} ${args} NUM_PROCS 1 TIMEOUT 60 LABELS tps
                   CTEST_NAME tps.zero_order_${mode}_${variables})
  endforeach()
endforeach()
foreach(order forward reverse)
  # OPALX-linked custom-main tests need the normal launcher/uenv for MPI libraries.
  # Keep their existing executable and CTest names while sharing launch policy.
  add_opalx_test(tps_opalx_${order} SOURCES opalx_headers.cpp NUM_PROCS 1 TIMEOUT 60
                 LABELS tps CTEST_NAME tps.opalx_headers_${order})
  if(order STREQUAL reverse)
    target_compile_definitions(tps_opalx_${order} PRIVATE TPS_REVERSE_INCLUDE_ORDER)
  endif()
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
