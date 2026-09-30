# Prefer the headers bundled with the Kokkos built by IPPL. A broad imported
# include directory (for example a uenv MPI prefix) can also contain an older
# Kokkos/Desul installation and otherwise precede these SYSTEM directories.
# Keep upstream targets and link order unchanged. PUBLIC covers OPALX consumers;
# BUILD_INTERFACE prevents dependency checkout paths leaking into installed exports.
function(opalx_prefer_fetched_kokkos_headers target)
  if(NOT TARGET Kokkos::kokkoscore)
    return()
  endif()

  get_target_property(_kokkos_imported Kokkos::kokkoscore IMPORTED)
  if(_kokkos_imported)
    return()
  endif()

  get_target_property(_kokkos_system_includes Kokkos::kokkoscore INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
  if(_kokkos_system_includes)
    target_include_directories(
      ${target} SYSTEM BEFORE PUBLIC "$<BUILD_INTERFACE:${_kokkos_system_includes}>"
    )
  endif()
endfunction()
