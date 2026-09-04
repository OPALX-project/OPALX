set (_SRCS
    SpaceChargeConfig.cpp
    SpaceChargeConfigBuilder.cpp
    SpaceChargeFactory.cpp
    SpaceChargeSolver.cpp
    SpaceChargeSolveContext.cpp
    )

set (HDRS
    SpaceChargeConfig.h
    SpaceChargeSolveContext.h
    )

set (_PRIVATE_HDRS
    SpaceChargeAlgorithm.h
    SpaceChargeSolver.h
    SpaceChargeFrames.h
    SpaceChargeConfigBuilder.h
    SpaceChargeFactory.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})

add_subdirectory (Poisson)
add_subdirectory (FFT2D5)
add_subdirectory (CartesianPIC)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge")
