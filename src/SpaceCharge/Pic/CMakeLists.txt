set (_SRCS
    PicDomainManager.cpp
    PicParticleDomainAdapter.cpp
    PicWorkspace.cpp
    )

set (_PRIVATE_HDRS
    PicDomainManager.h
    PicParticleDomainAdapter.h
    )

set (HDRS
    PicWorkspace.h
    PicWorkspace.tpp
    )

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge/Pic")
