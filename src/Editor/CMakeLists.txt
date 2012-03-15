set (_SRCS
  Edit.cpp
  EditCmd.cpp
  EditCycle.cpp
  EditEnd.cpp
  EditFlatten.cpp
  EditInstall.cpp
  EditMove.cpp
  EditParser.cpp
  EditReflect.cpp
  EditRemove.cpp
  EditReplace.cpp
  EditSelect.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
