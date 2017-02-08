set (_SRCS
  Dynamic.cpp
  MakeSequence.cpp
  SetIntegrator.cpp
  Static.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_opal_sources(${_SRCS})

set (HDRS
    Dynamic.h
    MakeSequence.h
    SetIntegrator.h
    Static.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/PhysicsActions")