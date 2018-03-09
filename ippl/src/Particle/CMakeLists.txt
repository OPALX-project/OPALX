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
  IntTSC.h
  IntSUDS.h
  NoParticleCachingPolicy.h
  PairBuilder/BasicPairBuilder.h
  PairBuilder/HashPairBuilder.h
  PairBuilder/HashPairBuilderPeriodic.h
  PairBuilder/PairConditions.h
  PairBuilder/SortingPairBuilder.h
  ParticleAttribBase.h
  ParticleAttrib.hpp
  ParticleAttribElem.h
  ParticleAttrib.h
  ParticleBalancer.hpp
  ParticleBalancer.h
  IpplParticleBase.hpp
  IpplParticleBase.h
  ParticleBConds.h
  ParticleCashedLayout.hpp
  ParticleCashedLayout.h
  ParticleInteractAttrib.hpp
  ParticleInteractAttrib.h
  ParticleInteractLayout.hpp
  ParticleInteractLayout.h
  ParticleLayout.hpp
  ParticleLayoutFromGrid.h
  ParticleLayout.h
  ParticleSpatialLayout.hpp
  ParticleSpatialLayout.h
  ParticleUniformLayout.hpp
  ParticleUniformLayout.h
  PAssign.hpp
  PAssignDefs.h
  PAssign.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_IPPL_SOURCES (${_SRCS})
ADD_IPPL_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/Particle)

ADD_SUBDIRECTORY (PairBuilder)
