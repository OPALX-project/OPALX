set (_SRCS
  Beamline.cpp
  BeamlineGeometry.cpp
  ElmPtr.cpp
  FlaggedBeamline.cpp
  FlaggedElmPtr.cpp
  SimpleBeamline.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_cl_sources(${_SRCS})