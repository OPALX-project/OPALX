set (_SRCS
  BeamBeamRep.cpp
  CollimatorRep.cpp
  CorrectorRep.cpp
  CyclotronRep.cpp
  DegraderRep.cpp
  DriftRep.cpp
  MarkerRep.cpp
  MonitorRep.cpp
  MultipoleRep.cpp
  Octupole.cpp
  PatchRep.cpp
  ProbeRep.cpp
  Quadrupole.cpp
  RBendRep.cpp
  RFCavityRep.cpp
  TravelingWaveRep.cpp
  SBendRep.cpp
  SeparatorRep.cpp
  SeptumRep.cpp
  Sextupole.cpp
  SkewOctupole.cpp
  SkewQuadrupole.cpp
  SkewSextupole.cpp
  SolenoidRep.cpp
  XCorrectorRep.cpp
  XMonitorRep.cpp
  YCorrectorRep.cpp
  YMonitorRep.cpp
  ParallelPlateRep.cpp
  CyclotronValleyRep.cpp
  StripperRep.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_cl_sources(${_SRCS})
