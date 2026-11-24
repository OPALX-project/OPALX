# IPPL owns Kokkos backend selection. TPS needs dense Kokkos Kernels but does not
# enable IPPL's optional eigenanalysis providers or fetch another Kokkos runtime.
if(NOT TARGET Kokkos::kokkos)
  message(FATAL_ERROR "TPS must be added after OPALX/IPPL configures Kokkos")
endif()
if(TARGET Kokkos::kokkoskernels)
  return()
endif()

set(KokkosKernels_VERSION "5.2.0" CACHE STRING "Kokkos Kernels release or git.tag/sha")
if(NOT KokkosKernels_VERSION MATCHES "^git\\.")
  # An in-tree Kokkos target is sufficient; the package config checks this flag.
  set(_tps_kokkos_found "${Kokkos_FOUND}")
  set(Kokkos_FOUND TRUE)
  find_package(KokkosKernels ${KokkosKernels_VERSION} CONFIG QUIET)
  set(Kokkos_FOUND "${_tps_kokkos_found}")
endif()
if(NOT TARGET Kokkos::kokkoskernels)
  include(FetchContent)
  set(KokkosKernels_ADD_DEFAULT_ETI OFF CACHE BOOL "Instantiate kernels used by TPS")
  set(KokkosKernels_ENABLE_COMPONENT_BLAS ON CACHE BOOL "TPS dense BLAS")
  set(KokkosKernels_ENABLE_COMPONENT_BATCHED ON CACHE BOOL "TPS small pivoted solves")
  set(KokkosKernels_ENABLE_SUPERNODAL_SPTRSV OFF CACHE BOOL "No sparse solver ETI")
  string(REGEX REPLACE "^git\\." "" _tps_kernels_ref "${KokkosKernels_VERSION}")
  FetchContent_Declare(KokkosKernels
    GIT_REPOSITORY https://github.com/kokkos/kokkos-kernels.git
    GIT_TAG "${_tps_kernels_ref}")
  FetchContent_MakeAvailable(KokkosKernels)
endif()
