set (_SRCS
    BuildElements.cpp
    BuildTracking.cpp
    DecayChain.cpp
    FreezeRun.cpp
    ObjectValues.cpp
    RunDescriptionWriter.cpp
    TrackCmd.cpp
    TrackEnd.cpp
    TrackRun.cpp
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

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Track")
