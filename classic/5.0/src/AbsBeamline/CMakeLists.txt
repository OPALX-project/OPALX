set (_SRCS
  AlignWrapper.cpp
  AttributeSet.cpp
  BeamBeam.cpp
  BeamlineVisitor.cpp
  Collimator.cpp
  Component.cpp
  Corrector.cpp
  Cyclotron.cpp
  CyclotronValley.cpp
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
  ParallelPlate.cpp
  Patch.cpp
  Probe.cpp
  RBend.cpp
  RFCavity.cpp
  RFQuadrupole.cpp
  Ring.cpp
  SBend.cpp
  SBend3D.cpp
  Separator.cpp
  Septum.cpp
  Solenoid.cpp
  Stripper.cpp
  TravelingWave.cpp
  VariableRFCavity.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_cl_sources(${_SRCS})
