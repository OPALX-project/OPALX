set (_SRCS
    BarnesHutAlgorithm.cpp
    )

set (_PRIVATE_HDRS
    BarnesHutAlgorithm.h
    )

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})
