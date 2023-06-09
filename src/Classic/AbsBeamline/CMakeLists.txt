add_subdirectory(EndFieldModel)
add_subdirectory(MultipoleTFunctions)

set (_SRCS
    AttributeSet.cpp
    BeamlineVisitor.cpp
    Component.cpp
    Drift.cpp
    ElementBase.cpp
    Marker.cpp
    Monitor.cpp
    Multipole.cpp
    MultipoleT.cpp
    MultipoleTBase.cpp
    MultipoleTStraight.cpp
    MultipoleTCurvedConstRadius.cpp
    MultipoleTCurvedVarRadius.cpp
    RFCavity.cpp
    TravelingWave.cpp
    )

set (HDRS
    AttributeSet.h
    BeamlineVisitor.h
    Component.h
    Drift.h
    ElementBase.h
    Marker.h
    Monitor.h
    Multipole.h
    MultipoleT.h
    MultipoleTBase.h
    MultipoleTStraight.h
    MultipoleTCurvedConstRadius.h
    MultipoleTCurvedVarRadius.h
    RFCavity.h
    TravelingWave.h
    SpecificElementVisitor.h
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/AbsBeamline")
