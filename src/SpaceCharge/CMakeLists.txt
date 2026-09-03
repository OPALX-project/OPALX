set (_SRCS
    ParticleSetView.cpp
    SelfFieldConfig.cpp
    SelfFieldConfigBuilder.cpp
    SelfFieldFactory.cpp
    SelfFieldRequestPolicy.cpp
    SelfFieldSystem.cpp
    SolveContext.cpp
    )

set (HDRS
    ParticleSetView.h
    ParticleStorageConfig.h
    SelfFieldAlgorithm.h
    SelfFieldConfig.h
    SelfFieldDiagnostics.h
    SelfFieldRequestPolicy.h
    SelfFieldSystem.h
    SolveContext.h
    SolverCapabilities.h
    )

set (_PRIVATE_HDRS
    ParticleFrameGuard.h
    SelfFieldConfigBuilder.h
    SelfFieldFactory.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})

add_subdirectory (Ippl)
add_subdirectory (Pic2d5)
add_subdirectory (Pic3D)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge")
