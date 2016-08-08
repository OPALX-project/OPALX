set (_SRCS
  NumToStr.cpp
  OpalException.cpp
  OpalFilter.cpp
  OpalSection.cpp
  RegularExpression.cpp
  Round.cpp
  Timer.cpp
  TpsWerrf.cpp
  Truncate.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
