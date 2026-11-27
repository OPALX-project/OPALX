set (_SRCS
    NumericalIntegrator.cpp
    )

set (HDRS
    BorisIntegrators.h
    BorisPusher.h
    NumericalIntegrator.h
    )

opalx_include_source_directory()

target_sources(opalx_run PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Steppers")