set (_SRCS
    ParticleBinTraversal.cpp
    RelativisticFieldComposer.cpp
    CartesianPIC3DAlgorithm.cpp
    CartesianDomainUpdater.cpp
    ParticleMeshFieldTransfer.cpp
    CartesianPIC3DFieldStorage.cpp
    )

set (_PRIVATE_HDRS
    ParticleBinTraversal.h
    CartesianPIC3DAlgorithm.h
    CartesianDomainUpdater.h
    RelativisticFieldComposer.h
    FieldMirror.hpp
    ParticleMeshFieldTransfer.h
    CartesianPIC3DFieldStorage.h
    CartesianPIC3DFieldStorage.tpp
    )

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})
