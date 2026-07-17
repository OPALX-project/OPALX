set (_SRCS
    FieldComposer.cpp
    PicDomainManager.cpp
    PicParticleDomainAdapter.cpp
    PicScatterGather.cpp
    PicWorkspace.cpp
    )

set (_PRIVATE_HDRS
    PicDomainManager.h
    PicParticleDomainAdapter.h
    )

set (HDRS
    FieldComposer.h
    FieldComposer.tpp
    PicScatterGather.h
    PicScatterGather.tpp
    PicWorkspace.h
    PicWorkspace.tpp
    )

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge/Pic")
