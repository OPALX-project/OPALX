SET (_SRCS
  )

SET (_HDRS
  AntiSymTenzor.h
  AppTypeTraits.h
  dcomplex.h
  GenVektor.h
  SymTenzor.h
  Tenzor.h
  TSVMetaAssign.h
  TSVMetaBinary.h
  TSVMetaCompare.h
  TSVMetaCross.h
  TSVMetaDotDot.h
  TSVMetaDot.h
  TSVMeta.h
  TSVMetaUnary.h
  Vektor.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_IPPL_SOURCES (${_SRCS})
ADD_IPPL_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/AppTypes)
