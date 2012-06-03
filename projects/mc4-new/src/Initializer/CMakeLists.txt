set (_SRCS
Cosmology.cpp
DataBase.cpp
IPPLInitializer.cpp
Initializer.cpp
InputParser.cpp
MT_Random.cpp
Output.cpp
Parallelization.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})

