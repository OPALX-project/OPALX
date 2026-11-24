set (_SRCS
    PartBunch.cpp
    BunchStateHandler.cpp
    CartesianDomain.cpp
)

set (HDRS
    PartBunch.h
    PartBunchFwd.h
    BunchStateHandler.h
    CartesianDomainConfig.h
    CartesianDomain.h
    CartesianDomain.tpp
    ParticleContainerTypes.h
    ParticleContainer.hpp
    ParticleFrame.h
    BCHandler.hpp
)

opalx_include_source_directory()

target_sources(opalx_run PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/PartBunch")

add_subdirectory (Binning)
