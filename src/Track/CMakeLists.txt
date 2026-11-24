# Input side: the TRACK, RUN and ENDTRACK commands and freeze().
set (_INPUT_SRCS
    FreezeRun.cpp
    TrackCmd.cpp
    TrackEnd.cpp
    TrackRun.cpp
    )

# Build bridge, built into opalx: turns a frozen run description into elements and tracking
# objects through value-only static functions of the element, EMISSIONSOURCE and FIELDSOLVER
# commands.
set (_BRIDGE_SRCS
    BuildElements.cpp
    BuildTracking.cpp
    )

# Run side: the run description, its checks and its text dump.
set (_RUN_SRCS
    DecayChain.cpp
    ObjectValues.cpp
    RunDescriptionWriter.cpp
    ValidateRun.cpp
    )

set (HDRS
    BuildElements.h
    BuildTracking.h
    DecayChain.h
    FreezeRun.h
    ObjectValues.h
    RunDescription.h
    RunDescriptionWriter.h
    TrackCmd.h
    TrackEnd.h
    TrackRun.h
    ValidateRun.h
    )

opalx_include_source_directory()

target_sources(opalx PRIVATE ${_INPUT_SRCS} ${_BRIDGE_SRCS})
target_sources(opalx_run PRIVATE ${_RUN_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Track")
