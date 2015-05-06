SET (_SRCS
  DiscBuffer.cpp
  DiscConfig.cpp
  DiscMeta.cpp
  DiscParticleFunctions.cpp
  FieldDebugFunctions.cpp
  ParticleDebugFunctions.cpp
  Inform.cpp
  Pool.cpp
  IpplCounter.cpp
  IpplStats.cpp
  IpplTimings.cpp
  IpplMessageCounter.cpp
  IpplMemoryUsage.cpp
  IpplInfo.cpp
  RNGSimple.cpp
  RNGXDiv.cpp
  RNGXCI.cpp
  RandomNumberGen.cpp
  Timer.cpp
  Unique.cpp
  User.cpp
  UserList.cpp
  PAssert.cpp
  )

SET (_HDRS
  bstring.h
  DiscBuffer.h
  DiscConfig.h
  DiscField.hpp
  DiscField.h
  DiscMeta.h
  DiscParticle.hpp
  DiscParticle.h
  DiscType.h
  FieldBlock.hpp
  FieldBlock.h
  FieldDebug.hpp
  FieldDebug.h
  FieldDebugFunctions.h
  FieldDebugPrint.hpp
  FieldDebugPrint.h
  FieldPrint.hpp
  FieldPrint.h
  FieldView.hpp
  FieldView.h
  Inform.h
  IplCounter.h
  IpplCounter.h
  IpplInfo.h
  IpplMemoryUsage.h
  IpplMessageCounter.h
  IpplStats.h
  IpplTimings.h
  my_auto_ptr.h
  NamedObj.h
  ParticleDebug.hpp
  ParticleDebug.h
  ParticleDebugFunctions.h
  PAssert.h
  Pooled.hpp
  Pooled.h
  Pool.h
  Pstring.h
  RandomNumberGen.h
  RefBlock.h
  RefCounted.h
  RNGAssignDefs.h
  RNGBitReverse.h
  RNGLattice.h
  RNGRand.h
  RNGSimple.h
  RNGXCI.h
  RNGXDiv.h
  SequenceGen.h
  Timer.h
  Unique.h
  User.h
  UserList.h
  Vec.h
  vmap.hpp
  vmap.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/Utility)
