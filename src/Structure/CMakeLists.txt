# Input side: the BEAM, EMISSIONSOURCE, FIELDSOLVER and BINNING commands.
set (_INPUT_SRCS
  Beam.cpp
  EmissionSource.cpp
  FieldSolverCmd.cpp
  BinningCmd.cpp
)

# Run side: output writers, checkpoints, bounding box and mesh generator.
set (_RUN_SRCS
  BoundingBox.cpp
  CheckpointFile.cpp
  BinConfigWriter.cpp
  DirichletPlaneWriter.cpp
  DataSink.cpp
  EMFieldDumpWriter.cpp
  H5FileSpace.cpp
  H5PartWrapper.cpp
  H5PartWrapperForPT.cpp
  H5Writer.cpp
  MeshGenerator.cpp       
  LossDataSink.cpp        
  StatBaseWriter.cpp
  StatWriter.cpp
  SDDSColumn.cpp
  SDDSColumnSet.cpp
  SDDSWriter.cpp
)

set (HDRS
    Beam.h
    EmissionSource.h
    CheckpointFile.h
    BoundingBox.h
    BinConfigWriter.h
    DirichletPlaneWriter.h
    DataSink.h
    DataSinkOptions.h
    EMFieldDumpWriter.h
    FieldSolverCmd.h
    BinningCmd.h
    H5PartWrapperForPT.h
    H5PartWrapper.h
    H5Writer.h
    MeshGenerator.h
    LossDataSink.h
    StatBaseWriter.h
    StatWriter.h
    SDDSColumn.h
    SDDSColumnSet.h
    SDDSWriter.h
)

opalx_include_source_directory()

target_sources(opalx PRIVATE ${_INPUT_SRCS})
target_sources(opalx_run PRIVATE ${_RUN_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Structure")
