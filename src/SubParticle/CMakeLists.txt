SET (_SRCS
  )

SET (_HDRS
  SubParticleAssign.cpp
  SubParticleAssignDefs.h
  SubParticleAssign.h
  SubParticleAttrib.cpp
  SubParticleAttrib.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/SubParticle)
