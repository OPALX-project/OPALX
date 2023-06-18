set (_SRCS
    OpalCavity.cpp
    OpalTravelingWave.cpp
    OpalDrift.cpp
    OpalElement.cpp
    OpalMarker.cpp
    OpalMultipole.cpp
    OpalMultipoleT.cpp
    OpalBeamline.cpp
    )

set (HDRS
    OpalBeamline.h
    OpalCavity.h
    OpalDrift.h
    OpalElement.h
    OpalMarker.h
    OpalMultipole.h
    OpalMultipoleT.h
    OpalSource.h
    OpalTravelingWave.h
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Elements")