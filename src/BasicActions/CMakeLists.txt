set (_SRCS
  Call.cpp
  Dump.cpp
  Echo.cpp
  Help.cpp
  Option.cpp
  Save.cpp
  Select.cpp
  Show.cpp
  Stop.cpp
  Quit.cpp
  System.cpp
  Title.cpp
  Value.cpp
  What.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
