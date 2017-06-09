set (_SRCS
    OffsetTest.cpp
    RingTest.cpp
    SBend3DTest.cpp
    SpiralSectorTest.cpp
    VariableRFCavityTest.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
