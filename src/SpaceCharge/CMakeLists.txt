set (_SRCS
    LegacyPic3DAlgorithm.cpp
    ParticleSetView.cpp
    SelfFieldConfig.cpp
    SelfFieldConfigBuilder.cpp
    SelfFieldDiagnostics.cpp
    SelfFieldFactory.cpp
    SelfFieldSystem.cpp
    SolveContext.cpp
    )

set (HDRS
    ParticleSetView.h
    SelfFieldAlgorithm.h
    SelfFieldConfig.h
    SelfFieldConfigBuilder.h
    SelfFieldDiagnostics.h
    SelfFieldFactory.h
    SelfFieldSystem.h
    SolveContext.h
    SolverCapabilities.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge")
