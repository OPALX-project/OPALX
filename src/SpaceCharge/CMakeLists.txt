set (_SRCS
    ParticleFieldSet.cpp
    SpaceChargeConfig.cpp
    SpaceChargeConfigBuilder.cpp
    SpaceChargeSolverFactory.cpp
    SpaceChargeRequestSchedule.cpp
    SpaceChargeSolver.cpp
    SpaceChargeSolveContext.cpp
    )

set (HDRS
    ParticleFieldSet.h
    SpaceChargeAlgorithm.h
    SpaceChargeConfig.h
    SpaceChargeDiagnostics.h
    SpaceChargeRequestSchedule.h
    SpaceChargeSolver.h
    SpaceChargeSolveContext.h
    SpaceChargeCapabilities.h
    )

set (_PRIVATE_HDRS
    SpaceChargeFrameGuard.h
    SpaceChargeConfigBuilder.h
    SpaceChargeSolverFactory.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})

add_subdirectory (Poisson)
add_subdirectory (FFT2D5)
add_subdirectory (CartesianPIC)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge")
