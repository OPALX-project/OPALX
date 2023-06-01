set (_SRCS
    CavityAutophaser.cpp
#    Ctunes.cpp
    IndexMap.cpp
#    Hamiltonian.cpp
#    lomb.cpp
#    MapAnalyser.cpp
    OrbitThreader.cpp
    ParallelTTracker.cpp
    StepSizeConfig.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    CavityAutophaser.h
#    Ctunes.h
#    Hamiltonian.h
    IndexMap.h
#    lomb.h
#    MapAnalyser.h
    OrbitThreader.h
    ParallelTTracker.h
    StepSizeConfig.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")
