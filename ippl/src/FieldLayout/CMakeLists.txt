SET (_SRCS
  FieldLayoutUser.cpp
  MultiBalancer.cpp
  VRB.cpp
  )

SET (_HDRS
  BinaryBalancer.hpp
  BinaryBalancer.h
  CenteredFieldLayout.hpp
  CenteredFieldLayout.h
  ConejoBalancer.hpp
  ConejoBalancer.h
  FieldLayout.hpp
  FieldLayout.h
  FieldLayoutUser.h
  MultiBalancer.h
  Vnode.h
  VnodeMultiBalancer.hpp
  VnodeMultiBalancer.h
  VRB.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_IPPL_SOURCES (${_SRCS})
ADD_IPPL_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/FieldLayout)
