set (_SRCS
  Dynamic.cpp
  MakeSequence.cpp
  SetIntegrator.cpp
  Static.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
