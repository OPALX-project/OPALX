set (_SRCS
    PartBunch.cpp
    FieldSolver.cpp
    BinnedFieldSolver.cpp
    BunchStateHandler.cpp
)

set (HDRS
    PartBunch.h
    PartBunchFwd.h
    BinContext.hpp
    BunchStateHandler.h
    FieldSolver.hpp
    FieldAccumulator.hpp
    FieldContainer.hpp
    FieldMirror.hpp
    FieldOps.hpp
    LoadBalancer.hpp
    ParticleContainer.hpp
    BCHandler.hpp
    BinnedFieldSolver.h
    Corrections/BoundaryCorrection.hpp
    Corrections/ImageChargeCorrection.hpp
    Corrections/ShiftedGreensCorrection.hpp
    Drivers/BinnedLorentzSelfFieldDriver.hpp
    Drivers/MonolithicSelfFieldDriver.hpp
    Drivers/SelfFieldDriver.hpp
)

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/PartBunch")

add_subdirectory (Binning)
add_subdirectory (Diagnostics)
add_subdirectory (Solvers)
