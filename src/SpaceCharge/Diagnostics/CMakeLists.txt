set (_SRCS
    FieldDiagnostics.cpp
)

set (HDRS
    FieldDiagnostics.hpp
)

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge/Diagnostics")
