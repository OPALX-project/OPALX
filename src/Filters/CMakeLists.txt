set (_SRCS
    SavitzkyGolay.cpp
    Stencil.cpp
    FixedFFTLowPass.cpp
    RelativeFFTLowPass.cpp
    )

set (HDRS
    Filter.h
    Filters.h
    FixedFFTLowPass.h
    RelativeFFTLowPass.h
    SavitzkyGolay.h
    Stencil.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Filters")
