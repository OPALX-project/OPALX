set (_SRCS
    )

set (HDRS
    BorisPusher.h
    LF2.h
    LF2.hpp
    RK4.h
    RK4.hpp
    Stepper.h
    Steppers.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Steppers")