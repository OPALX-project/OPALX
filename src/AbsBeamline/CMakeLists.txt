# AbsBeamline/CMakeLists.txt

# Include subdirectories

# List local source and header files
set(_SRCS
    Box.cpp
    Collimator.cpp
    ConstantEFieldCavity.cpp
    ConstantFocusing.cpp
    Drift.cpp
    FieldmapElement.cpp
    ElementBase.cpp
    Laser.cpp
    Marker.cpp
    Monitor.cpp
    Multipole.cpp
    MultipoleT.cpp
    MultipoleTBase.cpp
    MultipoleTStraight.cpp
    MultipoleTCurvedConstRadius.cpp
    RFCavity.cpp
    Solenoid.cpp
    TravelingWave.cpp
    RBend.cpp
    VariableRFCavity.cpp
    SBend.cpp
)

set(_HDRS
    Box.h
    Collimator.h
    ConstantEFieldCavity.h
    ConstantFocusing.h
    Drift.h
    FieldmapElement.h
    ElementBase.h
    Laser.h
    Marker.h
    Multipole.h
    MultipoleT.h
    MultipoleTBase.h
    MultipoleTStraight.h
    MultipoleTCurvedConstRadius.h
    RFCavity.h
    TravelingWave.h
    Solenoid.h
    RBend.h
    SBend.h
    VariableRFCavity.h
    MultipoleTConfig.h
)

# Add local include directory
opalx_include_source_directory()

# Append sources to the main executable target
target_sources(opalx_run PRIVATE ${_SRCS})

# Optionally install headers
install(FILES ${_HDRS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/AbsBeamline)
