set (_SRCS
    BeamStrippingPhysics.cpp
    CSRIGFWakeFunction.cpp
    CSRWakeFunction.cpp
    GreenWakeFunction.cpp
    ScatteringPhysics.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    BeamStrippingPhysics.hh
    CSRIGFWakeFunction.hh
    CSRWakeFunction.hh
    GreenWakeFunction.hh
    ParticleMatterInteractionHandler.hh
    ScatteringPhysics.hh
    WakeFunction.hh
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Solvers")
