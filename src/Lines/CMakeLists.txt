set (_SRCS
    Line.cpp
    LineTemplate.cpp
    Replacer.cpp
    Sequence.cpp
    SequenceMember.cpp
    SequenceParser.cpp
    SequenceTemplate.cpp
    )

set (HDRS
    Line.h
    LineTemplate.h
    Replacer.h
    Sequence.h
    SequenceMember.h
    SequenceParser.h
    SequenceTemplate.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Lines")
