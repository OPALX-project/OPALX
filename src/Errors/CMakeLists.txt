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

add_sources(${_SRCS})
