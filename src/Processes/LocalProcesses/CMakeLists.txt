set(_SRCS
    LocalProcess.cpp
    OpalLocalProcess.cpp
    OpalTestProcess.cpp
    TestCounterProcess.cpp
)

set(HDRS
    LocalProcess.h
    OpalLocalProcess.h
    OpalTestProcess.h
    TestCounterProcess.h
)

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(opalx PRIVATE ${_SRCS})

install(FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Processes/LocalProcesses")
