set (_SRCS
  )

set (_HDRS
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

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
add_headers(${_HDRS})

install(FILES ${_HDRS} DESTINATION include/AppTypes)
