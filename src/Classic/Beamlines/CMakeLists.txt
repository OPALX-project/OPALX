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

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    BeamlineGeometry.h
    Beamline.h
    ElmPtr.h
    FlaggedBeamline.h
    FlaggedElmPtr.h
    SimpleBeamline.h
    TBeamline.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Beamlines")