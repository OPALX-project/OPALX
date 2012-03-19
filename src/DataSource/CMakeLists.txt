set (_SRCS
  DataSource.cpp
  DataConnect.cpp
  DataConnectCreator.cpp
  DataSource.cpp
  PtclAttribDataSource.cpp
  PtclBaseDataSource.cpp
  )

set (_HDRS
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

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
add_headers(${_HDRS})

install(FILES ${_HDRS} DESTINATION include/DataSource)
