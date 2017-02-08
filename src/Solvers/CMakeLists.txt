set (_SRCS
  P3MPoissonSolver.cpp
  FFTPoissonSolver.cpp
  FFTBoxPoissonSolver.cpp
  MGPoissonSolver.cpp
  RectangularDomain.cpp
  EllipticDomain.cpp
  ArbitraryDomain.cpp
  BoxCornerDomain.cpp
  TaperDomain.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_opal_sources(${_SRCS})

set (HDRS
    ArbitraryDomain.h
    BoxCornerDomain.h
    EllipticDomain.h
    FFTBoxPoissonSolver.h
    FFTPoissonSolver.h
    IrregularDomain.h
    MGPoissonSolver.h
    P3MPoissonSolver.h
    PoissonSolver.h
    RectangularDomain.h
    TaperDomain.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Solvers")