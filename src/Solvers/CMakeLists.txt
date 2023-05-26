set (_SRCS
    FFTBoxPoissonSolver.cpp
    FFTPoissonSolver.cpp
    )

set (HDRS
    FFTBoxPoissonSolver.h
    FFTPoissonSolver.h
    PoissonSolver.h
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Solvers")
