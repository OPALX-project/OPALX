SET (_SRCS
  )

SET (_HDRS
  SubParticleAssign.hpp
  SubParticleAssignDefs.h
  SubParticleAssign.h
  SubParticleAttrib.hpp
  SubParticleAttrib.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/SubParticle)
