set (_SRCS
    PoissonBackendRegistry.cpp
)

set (HDRS
    PoissonBackend.hpp
    PoissonBackendRegistry.hpp
    PoissonBackends.hpp
    SolveRequest.hpp
    SolverCapabilities.hpp
)

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge/Solvers")
