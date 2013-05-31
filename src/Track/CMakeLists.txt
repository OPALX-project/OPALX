set (_SRCS
  Track.cpp
  TrackCmd.cpp
  TrackEnd.cpp
  TrackParser.cpp
  TrackRun.cpp
  TrackSave.cpp
  TrackStart.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
