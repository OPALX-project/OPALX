# AbsBeamline/CMakeLists.txt

# Include subdirectories
add_subdirectory(EndFieldModel)
add_subdirectory(MultipoleTFunctions)

# List local source and header files
set(_SRCS
    AttributeSet.cpp
    BeamlineVisitor.cpp
    Component.cpp
    Drift.cpp
    ElementBase.cpp
    Marker.cpp
    Multipole.cpp
    MultipoleT.cpp
    MultipoleTBase.cpp
    MultipoleTStraight.cpp
    MultipoleTCurvedConstRadius.cpp
    MultipoleTCurvedVarRadius.cpp
    RFCavity.cpp
    Ring.cpp
    ScalingFFAMagnet.cpp        
    Solenoid.cpp
    TravelingWave.cpp
    Offset.cpp
    VerticalFFAMagnet.cpp
    Probe.cpp
    PluginElement.cpp
)

set(_HDRS
    AttributeSet.h
    BeamlineVisitor.h
    Component.h
    Drift.h
    ElementBase.h
    Marker.h
    Multipole.h
    MultipoleT.h
    MultipoleTBase.h
    MultipoleTStraight.h
    MultipoleTCurvedConstRadius.h
    MultipoleTCurvedVarRadius.h
    RFCavity.h
    TravelingWave.h
    SpecificElementVisitor.h
    ScalingFFAMagnet.h
    Ring.h
    Solenoid.h
    Offset.h
    VerticalFFAMagnet.h
    Probe.h
    PluginElement.h
)

# Add local include directory
target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Append sources to the main executable target
target_sources(opalx PRIVATE ${_SRCS})

# Optionally install headers
install(FILES ${_HDRS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/AbsBeamline)
