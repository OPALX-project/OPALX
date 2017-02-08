set (_SRCS
  AlignBase.cpp
  AlignHandler.cpp
  AlignReader.cpp
  AlignRemover.cpp
  AlignWriter.cpp
  Error.cpp
  ErrorAlign.cpp
  ErrorCmd.cpp
  ErrorComp.cpp
  ErrorEnd.cpp
  ErrorField.cpp
  ErrorParser.cpp
  ErrorPrint.cpp
  ErrorSave.cpp
  ErrorSelect.cpp
  MPBase.cpp
  MPHandler.cpp
  MPReader.cpp
  MPRemover.cpp
  MPWriter.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_opal_sources(${_SRCS})

set (HDRS
    AlignBase.h
    AlignHandler.h
    AlignReader.h
    AlignRemover.h
    AlignWriter.h
    ErrorAlign.h
    ErrorCmd.h
    ErrorComp.h
    ErrorEnd.h
    ErrorField.h
    Error.h
    ErrorParser.h
    ErrorPrint.h
    ErrorSave.h
    ErrorSelect.h
    MPBase.h
    MPHandler.h
    MPReader.h
    MPRemover.h
    MPWriter.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Errors")