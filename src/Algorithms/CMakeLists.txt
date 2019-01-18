set (_SRCS
  Hamiltonian.cpp
  IndexMap.cpp
  OrbitThreader.cpp
  BeamBeam3D.cpp
  LieMapper.cpp
  MapAnalyser.cpp
  NilTracker.cpp
  MPSplitIntegrator.cpp
  ThickMapper.cpp
  ThickTracker.cpp
  # AutophaseTracker.cpp
  CavityAutophaser.cpp
  Ctunes.cpp
  lomb.cpp
  ParallelTTracker.cpp
  ParallelCyclotronTracker.cpp
  ParallelSliceTracker.cpp
  StatisticalErrors.cpp
  StatisticalErrorsUtilities.cpp
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

add_opal_sources(${_SRCS})

set (HDRS
    #AutophaseTracker.h
    BeamBeam3D.h
    CavityAutophaser.h
    Ctunes.h
    Hamiltonian.h
    IndexMap.h
    LieMapper.h
    lomb.h
    MapAnalyser.h
    MPSplitIntegrator.h
    NilTracker.h
    OrbitThreader.h
    ParallelCyclotronTracker.h
    ParallelSliceTracker.h
    ParallelTTracker.h
    StatisticalErrors.h
    StatisticalErrorsUtilities.h
    ThickMapper.h
    ThickTracker.h
    TransportMapper.h
    bet/BetError.h
    bet/EnvelopeBunch.h
    bet/EnvelopeSlice.h
    bet/error.h
    bet/profile.h
    bet/math/functions.h
    bet/math/integrate.h
    bet/math/interpol.h
    bet/math/linfit.h
    bet/math/rk.h
    bet/math/root.h
    bet/math/savgol.h
    bet/math/sort.h
    bet/math/svdfit.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")