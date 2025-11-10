set (_SRCS
    Call.cpp
    DumpFields.cpp
    DumpEMFields.cpp
    Echo.cpp
    Help.cpp
    Option.cpp
    Select.cpp
    Stop.cpp
    Quit.cpp
    PSystem.cpp
    System.cpp
    Title.cpp
    Value.cpp
    )

set (HDRS
    Call.h
    DumpFields.h
    DumpEMFields.h
    Echo.h
    Help.h
    Option.h
    PSystem.h
    Quit.h
    Select.h
    Stop.h
    System.h
    Title.h
    Value.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/BasicActions")