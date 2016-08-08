set (_SRCS
    LossDataSink.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_cl_sources(${_SRCS})
