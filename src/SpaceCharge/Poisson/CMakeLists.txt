set (_SRCS
    PoissonSolver.cpp
    )

set (_PRIVATE_HDRS
    PoissonSolver.h
    PoissonSolverFactory.h
    NullPoissonAdapter.h
    OpenPoissonAdapter.h
    PeriodicPoissonAdapter.h
    P3MAdapters.h
    )

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})
