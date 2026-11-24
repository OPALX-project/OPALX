# Build bridge, built into opalx: turns the RUN description into the space-charge configuration
# through value-only functions of the BINNING, EMISSIONSOURCE and FIELDSOLVER commands.
set (_BRIDGE_SRCS
    SpaceChargeConfigBuilder.cpp
    SpaceChargeConfigBuilder.h
    )

# Run side.
set (_RUN_SRCS
    SpaceChargeConfig.cpp
    SpaceChargeFactory.cpp
    SpaceChargeSolveContext.cpp
    SpaceChargeSolver.cpp
    )

set (HDRS
    SpaceChargeConfig.h
    SpaceChargeSolveContext.h
    )

set (_PRIVATE_HDRS
    SpaceChargeAlgorithm.h
    SpaceChargeSolver.h
    SpaceChargeFactory.h
    )

opalx_include_source_directory()

target_sources(opalx PRIVATE ${_BRIDGE_SRCS})
target_sources(opalx_run PRIVATE ${_RUN_SRCS} ${_PRIVATE_HDRS})

add_subdirectory (Poisson)
add_subdirectory (CartesianPIC3D)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge")
