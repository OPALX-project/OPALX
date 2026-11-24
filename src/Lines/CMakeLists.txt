set (_SRCS
    EmissionSourceList.cpp
    Line.cpp
    LineTemplate.cpp
    )

set (HDRS
    EmissionSourceList.h
    Line.h
    LineTemplate.h
    )

opalx_include_source_directory()

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Lines")
