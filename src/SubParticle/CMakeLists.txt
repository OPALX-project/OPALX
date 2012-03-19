set (_SRCS
  )

set (_HDRS
  SubParticleAssign.cpp
  SubParticleAssignDefs.h
  SubParticleAssign.h
  SubParticleAttrib.cpp
  SubParticleAttrib.h
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
add_headers(${_HDRS})

install(FILES ${_HDRS} DESTINATION include/SubParticle)
