set (_SRCS
    PartBunch.cpp
    FieldSolver.cpp
)

set (HDRS
    PartBunch.h
    FieldSolver.hpp
    FieldContainer.hpp
    LoadBalancer.hpp
    ParticleContainer.hpp
    datatypes.h
)

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/PartBunch")

add_subdirectory (Binning)