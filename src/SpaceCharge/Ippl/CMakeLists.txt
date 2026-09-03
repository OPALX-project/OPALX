set (_SRCS
    IpplPoissonAdapter.cpp
    )

set (_PRIVATE_HDRS
    IpplPoissonAdapter.h
    IpplPoissonTypes.h
    )

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})
