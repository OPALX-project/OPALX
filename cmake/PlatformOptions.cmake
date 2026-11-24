# -----------------------------------------------------------------------------
# PlatformOptions.cmake
#
# Determines target-dependent compile-time flags (e.g., fencing behavior).
# -----------------------------------------------------------------------------

# Set default
set(TimerFences true)

# Determine from Kokkos device list if available
if(Kokkos_DEVICES)
  string(REPLACE ";SERIAL" "" KokkosAccelerators "${Kokkos_DEVICES}")
  string(REPLACE "SERIAL;" "" KokkosAccelerators "${KokkosAccelerators}")

  if("${KokkosAccelerators}" MATCHES ".+;.+")
    set(TimerFences false)
    message(
      STATUS "⏱️  Multiple Kokkos devices detected (${Kokkos_DEVICES}) — disabling timer fences.")
  else()
    message(STATUS "⏱️  Single Kokkos device detected (${Kokkos_DEVICES}) — enabling timer fences.")
  endif()
else()
  message(STATUS "⏱️  Kokkos_DEVICES not set — defaulting to TimerFences = true.")
endif()

# Define macro for use in source code, in both OPALX libraries (OPALX_LIBRARY_TARGETS)
foreach(_opalx_library IN LISTS OPALX_LIBRARY_TARGETS)
  target_compile_definitions(${_opalx_library} PUBLIC OPALX_ENABLE_TIMER_FENCES=${TimerFences})
endforeach()

# Field debug option - apply compile definition to the libraries
if(OPALX_FIELD_DEBUG)
  message(STATUS "⚠️  OPALX_FIELD_DEBUG enabled — field dumps will be emitted during simulation")
  foreach(_opalx_library IN LISTS OPALX_LIBRARY_TARGETS)
    target_compile_definitions(${_opalx_library} PUBLIC OPALX_FIELD_DEBUG)
  endforeach()
endif()