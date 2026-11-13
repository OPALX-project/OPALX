set (_SRCS
    CofCmd.cpp
    BuildElements.cpp
    DecayChain.cpp
    FreezeRun.cpp
    RunDescriptionWriter.cpp
    Track.cpp
    TrackCmd.cpp
    TrackEnd.cpp
    TrackParser.cpp
    TrackRun.cpp
    ValidateRun.cpp
    )

set (HDRS
    CofCmd.h
    BuildElements.h
    DecayChain.h
    FreezeRun.h
    RunDescription.h
    RunDescriptionWriter.h
    TrackCmd.h
    TrackEnd.h
    Track.h
    TrackParser.h
    TrackRun.h
    ValidateRun.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Track")
