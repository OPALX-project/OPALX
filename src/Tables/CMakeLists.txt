set (_SRCS
  AttList.cpp
  AttWriter.cpp
  CorrectionBase.cpp
  Eigen.cpp
  Envelope.cpp
  Insertion.cpp
  List.cpp
  MatrixCmd.cpp
  Micado.cpp
  Period.cpp
  RangeSelector.cpp
  Selector.cpp
  Survey.cpp
  TableTester.cpp
  ThreadAll.cpp
  ThreadBpm.cpp
  Twiss.cpp
  Twiss3.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
