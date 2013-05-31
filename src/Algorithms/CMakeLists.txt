set (_SRCS
  BeamBeam3D.cpp
  LieMapper.cpp
  MPSplitIntegrator.cpp
  ThickMapper.cpp
  ThickTracker.cpp
  ParallelTTracker.cpp
  ParallelCyclotronTracker.cpp
  ParallelSliceTracker.cpp
  TransportMapper.cpp
  bet/EnvelopeSlice.cpp
  bet/EnvelopeBunch.cpp
  bet/profile.cpp
  bet/math/sort.cpp
  bet/math/integrate.cpp
  bet/math/interpol.cpp
  bet/math/rk.cpp
  bet/math/functions.cpp
  bet/math/linfit.cpp
  bet/math/root.cpp
  bet/math/savgol.cpp
  bet/math/svdfit.cpp
  bet/BetError.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
