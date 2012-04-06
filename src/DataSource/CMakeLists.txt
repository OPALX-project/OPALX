SET (_SRCS
  DataSource.cpp
  DataConnect.cpp
  DataConnectCreator.cpp
  DataSource.cpp
  PtclAttribDataSource.cpp
  PtclBaseDataSource.cpp
  )

SET (_HDRS
  ACLVISDataConnect.h
  ACLVISFieldDataSource.cpp
  ACLVISFieldDataSource.h
  ACLVISOperations.h
  ACLVISPtclAttribDataSource.cpp
  ACLVISPtclAttribDataSource.h
  ACLVISPtclBaseDataSource.cpp
  ACLVISPtclBaseDataSource.h
  DataConnectCreator.h
  DataConnect.h
  DataSource.h
  DataSourceObject.h
  FieldDataSource.cpp
  FieldDataSource.h
  FileDataConnect.h
  FileFieldDataSource.cpp
  FileFieldDataSource.h
  FilePtclAttribDataSource.cpp
  FilePtclAttribDataSource.h
  FilePtclBaseDataSource.cpp
  FilePtclBaseDataSource.h
  MakeDataSource.cpp
  MakeDataSource.h
  PtclAttribDataSource.h
  PtclBaseDataSource.h
  ScalarDataSource.cpp
  ScalarDataSource.h
  StringDataSource.cpp
  StringDataSource.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/DataSource)
