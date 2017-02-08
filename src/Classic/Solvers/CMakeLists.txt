set (_SRCS
  CollimatorPhysics.cpp
  CSRWakeFunction.cpp
  CSRIGFWakeFunction.cpp
  GreenWakeFunction.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    GreenWakeFunction.hh
    CollimatorPhysics.hh
    CSRWakeFunction.hh
    SurfacePhysicsHandler.hh
    WakeFunction.hh
    CSRIGFWakeFunction.hh
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Solvers")