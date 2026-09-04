set (_SRCS
    PoissonSolver.cpp
    )

set (_PRIVATE_HDRS
    PoissonSolver.h
    )

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})
