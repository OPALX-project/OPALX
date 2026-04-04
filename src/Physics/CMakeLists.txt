set (_SRCS
    LinearCompton.cpp
    LinearBreitWheeler.cpp
    Material.cpp
    ParticleProperties.cpp
    )

set (HDRS
    LinearCompton.h
    LinearBreitWheeler.h
    Material.h
    ParticleProperties.h
    Physics.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Physics")