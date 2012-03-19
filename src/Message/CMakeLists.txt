set (_SRCS
  CRC.cpp
  Communicate.cpp
  CommCreator.cpp
  CommMPI.cpp
  MessageFunctions.cpp
  Formatter.cpp
  )

set (_HDRS
  CommCreator.h
  CommMPI.h
  CommPM.h
  CommSHMEMPI.h
  Communicate.h
  CRC.h
  Formatter.h
  GlobalComm.cpp
  GlobalComm.h
  Message.cpp
  Message.h
  TagMaker.h
  Tags.h
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
add_headers(${_HDRS})

install(FILES ${_HDRS} DESTINATION include/Message)
