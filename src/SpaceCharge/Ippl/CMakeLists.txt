set (_SRCS
    IpplPoissonAdapter.cpp
    )

set (HDRS
    IpplPoissonAdapter.h
    IpplPoissonTypes.h
    )

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge/Ippl")
