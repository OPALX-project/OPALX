set (_SRCS
    Configure.cpp
    )

set (HDRS
    Configure.h
    )

opalx_include_source_directory()

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpalConfigure")