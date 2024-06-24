set (_SRCS
    PartBunch.cpp
    FieldSolver.cpp
)

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    PartBunch.hpp    
    FieldSolver.hpp
    FieldContainer.hpp
    LoadBalancer.hpp
    ParticleContainer.hpp
    datatypes.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/PartBunch")