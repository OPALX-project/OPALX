SET (_SRCS
  )

SET (_HDRS
  BoxParticleCachingPolicy.h
  CellParticleCachingPolicy.h
  GenArrayParticle.h
  GenParticle.h
  IntCIC.h
  Interpolator.h
  IntNGP.h
  IntSUDS.h
  NoParticleCachingPolicy.h
  PairBuilder/BasicPairBuilder.h
  PairBuilder/HashPairBuilder.h
  PairBuilder/PairConditions.h
  PairBuilder/SortingPairBuilder.h
  ParticleAttribBase.h
  ParticleAttrib.cpp
  ParticleAttribElem.h
  ParticleAttrib.h
  ParticleBalancer.cpp
  ParticleBalancer.h
  ParticleBase.cpp
  ParticleBase.h
  ParticleBConds.h
  ParticleCashedLayout.cpp
  ParticleCashedLayout.h
  ParticleInteractAttrib.cpp
  ParticleInteractAttrib.h
  ParticleInteractLayout.cpp
  ParticleInteractLayout.h
  ParticleLayout.cpp
  ParticleLayoutFromGrid.h
  ParticleLayout.h
  ParticleSpatialLayout.cpp
  ParticleSpatialLayout.h
  ParticleUniformLayout.cpp
  ParticleUniformLayout.h
  PAssign.cpp
  PAssignDefs.h
  PAssign.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/Particle)

ADD_SUBDIRECTORY (PairBuilder)
