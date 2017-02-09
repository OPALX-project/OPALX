SET (_SRCS
  )

SET (_HDRS
  IpplExpressions.h
  IpplTypeComputations.h
  PETE.h
  TypeComputations.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_IPPL_SOURCES (${_SRCS})
ADD_IPPL_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/PETE)
