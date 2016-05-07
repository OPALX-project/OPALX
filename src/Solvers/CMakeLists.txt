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

add_sources(${_SRCS})
