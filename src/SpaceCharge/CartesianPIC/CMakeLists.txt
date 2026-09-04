set (_SRCS
    ParticleBinTraversal.cpp
    CorrectionPassSchedule.cpp
    RelativisticFieldComposer.cpp
    CartesianPICAlgorithm.cpp
    CartesianDomainUpdater.cpp
    ParticleDomainOperations.cpp
    ParticleMeshFieldTransfer.cpp
    CartesianPICFieldStorage.cpp
    )

set (_PRIVATE_HDRS
    ParticleBinTraversal.h
    ParticleBinTraversal.tpp
    P3MShortRangeInteraction.h
    CartesianPICAlgorithm.h
    CartesianDomainUpdater.h
    ParticleDomainOperations.h
    CorrectionPassSchedule.h
    RelativisticFieldComposer.h
    RelativisticFieldComposer.tpp
    FieldMirror.hpp
    ParticleMeshFieldTransfer.h
    ParticleMeshFieldTransfer.tpp
    CartesianPICFieldStorage.h
    CartesianPICFieldStorage.tpp
    CartesianPICPass.h
    )

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})
