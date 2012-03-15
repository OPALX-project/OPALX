set (_SRCS
  AbstractFun.cpp
  AbstractVar.cpp
  ConcreteFun.cpp
  ConcreteVar.cpp
  ConstraintCmd.cpp
  LMDif.cpp
  Match.cpp
  MatchCmd.cpp
  MatchEnd.cpp
  MatchOption.cpp
  MatchParser.cpp
  Migrad.cpp
  Simplex.cpp
  VaryCmd.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
