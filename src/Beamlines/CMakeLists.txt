set (_SRCS
    Beamline.cpp
    BeamlineGeometry.cpp
    ElmPtr.cpp
    FlaggedBeamline.cpp
    FlaggedElmPtr.cpp
    )

set (HDRS
    BeamlineGeometry.h
    Beamline.h
    ElmPtr.h
    FlaggedBeamline.h
    FlaggedElmPtr.h
    TBeamline.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Beamlines")