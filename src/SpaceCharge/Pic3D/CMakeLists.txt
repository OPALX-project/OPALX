set (_SRCS
    BinningPlan.cpp
    CorrectionPlan.cpp
    FieldFramePolicy.cpp
    FieldComposer.cpp
    IterationPlanFactory.cpp
    Pic3DSolver.cpp
    PicDomainManager.cpp
    PicParticleDomainAdapter.cpp
    PicScatterGather.cpp
    PicWorkspace.cpp
    WholeBunchPlan.cpp
    )

set (_PRIVATE_HDRS
    BinningPlan.h
    BinningPlan.tpp
    FieldFramePolicy.h
    IterationPlanFactory.h
    Pic3DSolver.h
    PicDomainManager.h
    PicParticleDomainAdapter.h
    WholeBunchPlan.h
    )

set (HDRS
    CorrectionPlan.h
    CorrectionPlan.tpp
    FieldComposer.h
    FieldComposer.tpp
    FieldMirror.hpp
    IterationPlan.h
    PicScatterGather.h
    PicScatterGather.tpp
    PicWorkspace.h
    PicWorkspace.tpp
    SolvePass.h
    )

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge/Pic3D")
