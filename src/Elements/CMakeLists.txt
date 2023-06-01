set (_SRCS
    OpalCavity.cpp
    OpalTravelingWave.cpp
    OpalDrift.cpp
    OpalElement.cpp
    OpalMarker.cpp
    OpalMonitor.cpp
    OpalMultipole.cpp
    OpalMultipoleT.cpp
    OpalVacuum.cpp
    OpalBeamline.cpp
    )

set (HDRS
    OpalBeamline.h
    OpalCavity.h
    OpalDrift.h
    OpalElement.h
    OpalMarker.h
    OpalMonitor.h
    OpalMultipole.h
    OpalMultipoleT.h
    OpalSource.h
    OpalTravelingWave.h
    OpalVacuum.h
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Elements")