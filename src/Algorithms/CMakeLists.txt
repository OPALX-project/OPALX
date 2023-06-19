set (_SRCS
    CavityAutophaser.cpp
    IndexMap.cpp
    OrbitThreader.cpp
    ParallelTTracker.cpp
    StepSizeConfig.cpp
    Tracker.cpp
)   

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    CavityAutophaser.h
    IndexMap.h
    OrbitThreader.h
    ParallelTTracker.h
    StepSizeConfig.h
    Tracker.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")
