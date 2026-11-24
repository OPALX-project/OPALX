set (_SRCS
    )

set (HDRS
    BorisPusher.h
    )

opalx_include_source_directory()

target_sources(opalx_run PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Steppers")