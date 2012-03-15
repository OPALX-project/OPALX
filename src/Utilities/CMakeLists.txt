set (_SRCS
  OpalException.cpp
  OpalFilter.cpp
  NumToStr.cpp
  RegularExpression.cpp
  Round.cpp
  Timer.cpp
  TpsWerrf.cpp
  Truncate.cpp
  OpalSection.cpp
  OpalField.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
