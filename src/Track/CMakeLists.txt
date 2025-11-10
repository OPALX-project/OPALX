set (_SRCS
    Track.cpp
    TrackCmd.cpp
    TrackEnd.cpp
    TrackParser.cpp
    TrackRun.cpp
    )

set (HDRS
    TrackCmd.h
    TrackEnd.h
    Track.h
    TrackParser.h
    TrackRun.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Track")