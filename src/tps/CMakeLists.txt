# TPS shares OPALX's selected Kokkos backend; no second runtime is configured.
include(${PROJECT_SOURCE_DIR}/cmake/TPSDependencies.cmake)

add_library(tps Algebra/TpsData.cpp Algebra/TpsMonomial.cpp)
add_library(TPS::tps ALIAS tps)
target_compile_features(tps PUBLIC cxx_std_20)
set_target_properties(tps PROPERTIES CXX_EXTENSIONS OFF POSITION_INDEPENDENT_CODE ON)
target_include_directories(tps PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/opalx>)
target_link_libraries(tps PUBLIC Kokkos::kokkos Kokkos::kokkoskernels)

install(TARGETS tps EXPORT OPALXTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES Copyright.readme DESTINATION ${CMAKE_INSTALL_DATADIR}/opalx/tps)
