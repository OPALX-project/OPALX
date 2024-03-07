set (_SRCS
  Beam.cpp
  BoundingBox.cpp
  BoundaryGeometry.cpp
  DataSink.cpp
  ElementPositionWriter.cpp
  FieldsolverCmd.cpp
  H5PartWrapper.cpp
  H5PartWrapperForPT.cpp
  MonitorStatisticsWriter.cpp
  IpplInfoWrapper.cpp
  H5Writer.cpp
  MeshGenerator.cpp       
  LBalWriter.cpp
  LossDataSink.cpp        
  StatBaseWriter.cpp
  StatWriter.cpp
  SDDSColumn.cpp
  SDDSColumnSet.cpp
  SDDSWriter.cpp
  PeakFinder.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

set (HDRS
    Beam.h
    BoundaryGeometry.h
    BoundingBox.h
    DataSink.h
    ElementPositionWriter.h
    FieldSolverCmd.h
    H5PartWrapperForPT.h
    H5PartWrapper.h
    IpplInfoWrapper.h
    MonitorStatisticsWriter.h
    H5Writer.h
    MeshGenerator.h
    LBalWriter.h
    LossDataSink.h
    StatBaseWriter.h
    StatWriter.h
    SDDSColumn.h
    SDDSColumnSet.h
    SDDSWriter.h
    PeakFinder.h
)

if (UNIX AND NOT APPLE)
    # for Linux, BSD, Solaris, Minix
    list (APPEND _SRCS MemoryProfiler.cpp)
    list (APPEND HDRS MemoryProfiler.h)
endif()

add_opal_sources (${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Structure")
