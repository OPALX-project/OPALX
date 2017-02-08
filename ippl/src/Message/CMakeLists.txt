SET (_SRCS
  CRC.cpp
  Communicate.cpp
  CommCreator.cpp
  CommMPI.cpp
  MessageFunctions.cpp
  Formatter.cpp
  )

SET (_HDRS
  CommCreator.h
  CommMPI.h
  CommPM.h
  CommSHMEMPI.h
  Communicate.h
  CRC.h
  Formatter.h
  GlobalComm.hpp
  GlobalComm.h
  Message.hpp
  Message.h
  TagMaker.h
  Tags.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_IPPL_SOURCES (${_SRCS})
ADD_IPPL_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/Message)
