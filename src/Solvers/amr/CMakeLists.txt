IF (ENABLE_AMR_SOLVER)
set (_SRCS
  Accel.cpp
  AccelBld.cpp
  AccelParticles.cpp
  Accel_advance.cpp
  Accel_error.cpp
  Accel_output.cpp
  Accel_setup.cpp
  Electrostatic.cpp
  ParticleDerive.cpp  
  Src_3d/Derive_3d.f90
  Src_3d/Electrostatic_3d.f90
  Src_3d/Nyx_3d.f90
  Src_3d/ProbInit.f90
  Src_3d/Tagging_3d.f90
  Src_3d/bc_fill_3d.f90
  TrilinosSolver.cpp
  init_trilinos.cpp
)
ENDIF (ENABLE_AMR_SOLVER)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
