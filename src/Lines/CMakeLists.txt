set (_SRCS
    EmissionSourceList.cpp
    Line.cpp
    LineTemplate.cpp
    Ring.cpp
    )

set (HDRS
    EmissionSourceList.h
    Line.h
    LineTemplate.h
    Ring.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Lines")
