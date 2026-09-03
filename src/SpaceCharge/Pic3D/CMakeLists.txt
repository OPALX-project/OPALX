set (_SRCS
    BinningPlan.cpp
    CorrectionPlan.cpp
    FieldComposer.cpp
    Pic3DSolver.cpp
    PicDomainManager.cpp
    PicParticleDomainAdapter.cpp
    PicScatterGather.cpp
    PicWorkspace.cpp
    )

set (_PRIVATE_HDRS
    BinningPlan.h
    BinningPlan.tpp
    P3MShortRangeInteraction.h
    Pic3DSolver.h
    PicDomainManager.h
    PicParticleDomainAdapter.h
    CorrectionPlan.h
    FieldComposer.h
    FieldComposer.tpp
    FieldMirror.hpp
    PicScatterGather.h
    PicScatterGather.tpp
    PicWorkspace.h
    PicWorkspace.tpp
    SolvePass.h
    )

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})
