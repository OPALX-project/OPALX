set (_SRCS
    Channel.cpp
    ConstChannel.cpp
    )

set (HDRS
    Channel.h
    ConstChannel.h
    DirectChannel.h
    IndexedChannel.h
    IndirectChannel.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Channels")