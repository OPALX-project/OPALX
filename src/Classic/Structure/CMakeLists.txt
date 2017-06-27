set (_SRCS
    LossDataSink.cpp
    PeakFinder.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
