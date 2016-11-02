SET (_SRCS
  DataSource.cpp
  DataConnect.cpp
  DataConnectCreator.cpp
  DataSource.cpp
  PtclAttribDataSource.cpp
  PtclBaseDataSource.cpp
  )

SET (_HDRS
  DataConnectCreator.h
  DataConnect.h
  DataSource.h
  DataSourceObject.h
  FieldDataSource.hpp
  FieldDataSource.h
  FileDataConnect.h
  FileFieldDataSource.hpp
  FileFieldDataSource.h
  FilePtclAttribDataSource.hpp
  FilePtclAttribDataSource.h
  FilePtclBaseDataSource.hpp
  FilePtclBaseDataSource.h
  MakeDataSource.hpp
  MakeDataSource.h
  PtclAttribDataSource.h
  PtclBaseDataSource.h
  ScalarDataSource.hpp
  ScalarDataSource.h
  StringDataSource.hpp
  StringDataSource.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/DataSource)
