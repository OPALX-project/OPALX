# Input side: the object model of the input language and OpalData.
set (_INPUT_SRCS
  Action.cpp
  Attribute.cpp
  AttributeBase.cpp
  AttributeHandler.cpp
  BeamSequence.cpp
  Definition.cpp
  Directory.cpp
  Element.cpp
  Invalidator.cpp
  OpalData.cpp
  Object.cpp
  ObjectFunction.cpp
  PlaceRep.cpp
  RangeRep.cpp
  ValueDefinition.cpp
  )

# Run side: OpalParticle, the runtime record of one particle.
set (_RUN_SRCS
  OpalParticle.cpp
  )

set (HDRS
    Action.h
    AttributeBase.h
    Attribute.h
    AttributeHandler.h
    BeamSequence.h
    Definition.h
    Directory.h
    Element.h
    Expressions.h
    Invalidator.h
    ObjectFunction.h
    Object.h
    OpalData.h
    OpalParticle.h
    PlaceRep.h
    RangeRep.h
    ValueDefinition.h
)

opalx_include_source_directory()

target_sources(opalx PRIVATE ${_INPUT_SRCS})
target_sources(opalx_run PRIVATE ${_RUN_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/AbstractObjects")