# A source-only cleanup for the fetched Kokkos Kernels 5.2.0 Getrf header.
# The unused device-implementation A0 view triggers HIP -Wunused-variable.
# Header hashes restrict this edit to the reviewed upstream file: installed
# packages, modified checkouts and other versions are left untouched. Remove
# this workaround when the selected upstream release includes the cleanup.
include_guard(GLOBAL)

function(opalx_patch_kokkos_getrf source_dir)
  if(NOT source_dir)
    return()
  endif()
  set(_header "${source_dir}/batched/dense/impl/KokkosBatched_Getrf_Serial_Internal.hpp")
  if(NOT EXISTS "${_header}")
    return()
  endif()
  set(_original_sha "c690e58f602d2aa944fd24cd5abbd627b5ccaa6c32030c6ebec04aea1d4b9e47")
  set(_patched_sha "651292ba7aa5ed59f8d187539ec2b7947b2468bfce4e10dcd2dbebb62630759f")
  file(SHA256 "${_header}" _actual_sha)
  if(_actual_sha STREQUAL _patched_sha)
    return()
  elseif(NOT _actual_sha STREQUAL _original_sha)
    message(STATUS "Kokkos Getrf header differs from the reviewed 5.2.0 file; leaving it unchanged")
    return()
  endif()

  file(READ "${_header}" _original)
  # A_current identifies the device implementation; the host implementation's
  # A0 is required for recursion and must remain in place.
  string(REPLACE
    "    auto A0    = Kokkos::subview(A_current, Kokkos::ALL, Kokkos::pair<int, int>(0, n1));\n"
    "" _patched "${_original}")
  string(SHA256 _result_sha "${_patched}")
  if(NOT _result_sha STREQUAL _patched_sha)
    message(FATAL_ERROR "Kokkos Getrf cleanup did not produce the reviewed header")
  endif()
  file(WRITE "${_header}" "${_patched}")
  message(STATUS "Applied Kokkos Kernels 5.2.0 unused device A0 cleanup")
endfunction()
