set (_SRCS
  AmrBoxLib.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_opal_sources(${_SRCS})

set (HDRS
    AmrObject.h
    AmrBoxLib.h
    BoxLibLayout.h
    BoxLibLayout.hpp
    BoxLibParticle.h
    BoxLibParticle.hpp
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Amr")