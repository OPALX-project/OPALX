set (_SRCS
  Beam.cpp
  EmissionSource.cpp
  BoundingBox.cpp
  CheckpointFile.cpp
  BinConfigWriter.cpp
  DirichletPlaneWriter.cpp
  DataSink.cpp
  EMFieldDumpWriter.cpp
  FieldSolverCmd.cpp
  BinningCmd.cpp
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

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Structure")
