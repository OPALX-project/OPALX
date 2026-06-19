set (_SRCS
    PartBunch.cpp
    BunchStateHandler.cpp
)

set (HDRS
    PartBunch.h
    PartBunchFwd.h
    BunchStateHandler.h
    FieldContainer.hpp
    FieldStorage.hpp
    LoadBalancer.hpp
    ParticleContainer.hpp
)

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/PartBunch")
