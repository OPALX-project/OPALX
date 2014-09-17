set (_SRCS
  AlignWrapper.cpp
  AttributeSet.cpp
  BeamBeam.cpp
  BeamlineVisitor.cpp
  Collimator.cpp
  Component.cpp
  Corrector.cpp
  Cyclotron.cpp
  Degrader.cpp
  Diagnostic.cpp
  Drift.cpp
  ElementBase.cpp
  ElementImage.cpp
  Integrator.cpp
  Lambertson.cpp
  Marker.cpp
  Monitor.cpp
  Multipole.cpp
  Offset.cpp
  Patch.cpp
  Probe.cpp
  RBend.cpp
  RFCavity.cpp
  VariableRFCavity.cpp
  TravelingWave.cpp
  RFQuadrupole.cpp
  SBend.cpp
  SBend3D.cpp
  Separator.cpp
  Septum.cpp
  Solenoid.cpp
  ParallelPlate.cpp
  CyclotronValley.cpp
  Stripper.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_cl_sources(${_SRCS})
