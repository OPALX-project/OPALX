set (_SRCS
    Call.cpp
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

opalx_include_source_directory()

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/BasicActions")