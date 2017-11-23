set (_SRCS
  Beam.cpp
  OpalWake.cpp
  ParticleMatterInteraction.cpp
  SecondaryEmissionPhysics.cpp
  PriEmissionPhysics.cpp
  BoundaryGeometry.cpp
  FieldSolver.cpp
  DataSink.cpp
  H5PartWrapper.cpp
  H5PartWrapperForPC.cpp
  H5PartWrapperForPS.cpp
  H5PartWrapperForPT.cpp
  SDDSParser.cpp
  SDDSParser/array.cpp
  SDDSParser/associate.cpp
  SDDSParser/ast.cpp
  SDDSParser/column.cpp
  SDDSParser/data.cpp
  SDDSParser/description.cpp
  SDDSParser/file.cpp
  SDDSParser/include.cpp
  SDDSParser/parameter.cpp
  SDDSParser/version.cpp
  OpalInputInterpreter.cpp
  IpplInfoWrapper.cpp
  TracerParticles.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_opal_sources(${_SRCS})

set (HDRS
    Beam.h
    BoundaryGeometry.h
    DataSink.h
    FieldSolver.h
    H5PartWrapperForPC.h
    H5PartWrapperForPS.h
    H5PartWrapperForPT.h
    H5PartWrapper.h
    IpplInfoWrapper.h
    OpalInputInterpreter.h
    OpalWake.h
    PriEmissionPhysics.h
    SDDSParser.h
    SecondaryEmissionPhysics.h
    ParticleMatterInteraction.h
    SDDSParser/array_def.hpp
    SDDSParser/array.hpp
    SDDSParser/associate_def.hpp
    SDDSParser/associate.hpp
    SDDSParser/ast_def.hpp
    SDDSParser/ast.hpp
    SDDSParser/column_def.hpp
    SDDSParser/column.hpp
    SDDSParser/data_def.hpp
    SDDSParser/data.hpp
    SDDSParser/description_def.hpp
    SDDSParser/description.hpp
    SDDSParser/error_handler.hpp
    SDDSParser/file_def.hpp
    SDDSParser/file.hpp
    SDDSParser/include_def.hpp
    SDDSParser/include.hpp
    SDDSParser/parameter_def.hpp
    SDDSParser/parameter.hpp
    SDDSParser/skipper.hpp
    SDDSParser/version_def.hpp
    SDDSParser/version.hpp
    TracerParticles.hpp
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Structure")