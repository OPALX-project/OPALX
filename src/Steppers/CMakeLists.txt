set (_SRCS
    NumericalIntegrator.cpp
    RungeKuttaIntegrator.cpp
    )

set (HDRS
    BorisIntegrators.h
    BorisPusher.h
    NumericalIntegrator.h
    PhaseSpaceRate.h
    RungeKuttaIntegrator.h
    RungeKuttaTableau.h
    )

opalx_include_source_directory()

target_sources(opalx_run PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Steppers")