set (_SRCS
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

add_subdirectory (Ippl)
add_subdirectory (Pic)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge")
