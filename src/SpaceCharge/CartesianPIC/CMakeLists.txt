set (_SRCS
    ParticleBinTraversal.cpp
    RelativisticFieldComposer.cpp
    CartesianPICAlgorithm.cpp
    CartesianDomainUpdater.cpp
    ParticleMeshFieldTransfer.cpp
    CartesianPICFieldStorage.cpp
    )

set (_PRIVATE_HDRS
    ParticleBinTraversal.h
    P3MShortRangeInteraction.h
    CartesianPICAlgorithm.h
    CartesianDomainUpdater.h
    RelativisticFieldComposer.h
    FieldMirror.hpp
    ParticleMeshFieldTransfer.h
    CartesianPICFieldStorage.h
    CartesianPICFieldStorage.tpp
    )

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})
