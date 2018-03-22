set (_SRCS
    AmrBoxLib.cpp
    AmrObject.cpp
    AmrYtWriter.cpp
    Tagging_nd.f90
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_opal_sources(${_SRCS})

set (HDRS
    AbstractAmrWriter.h
    AmrBoxLib.h
    AmrBoxLib_F.h
    AmrObject.h
    AmrYtWriter.h
    BoxLibLayout.h
    BoxLibLayout.hpp
    BoxLibParticle.h
    BoxLibParticle.hpp
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Amr")