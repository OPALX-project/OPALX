set (_SRCS
    BinningPlan.cpp
    CorrectionPlan.cpp
    FieldComposer.cpp
    IterationPlanFactory.cpp
    PicDomainManager.cpp
    PicParticleDomainAdapter.cpp
    PicScatterGather.cpp
    PicWorkspace.cpp
    WholeBunchPlan.cpp
    )

set (_PRIVATE_HDRS
    BinningPlan.h
    BinningPlan.tpp
    IterationPlanFactory.h
    PicDomainManager.h
    PicParticleDomainAdapter.h
    WholeBunchPlan.h
    )

set (HDRS
    CorrectionPlan.h
    CorrectionPlan.tpp
    FieldComposer.h
    FieldComposer.tpp
    IterationPlan.h
    PicScatterGather.h
    PicScatterGather.tpp
    PicWorkspace.h
    PicWorkspace.tpp
    SolvePass.h
    )

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge/Pic")
