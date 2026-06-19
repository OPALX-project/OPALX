set (_SRCS
    FieldSolver.cpp
    BinnedFieldSolver.cpp
)

set (HDRS
    BCHandler.hpp
    BinContext.hpp
    BinnedFieldSolver.h
    BinnedFieldSolver.tpp
    FieldAccumulator.hpp
    FieldMirror.hpp
    FieldOps.hpp
    FieldSolver.hpp
    FieldSolverConfigBuilder.hpp
    FieldSolverConfig.hpp
    ImageChargeScatterController.h
    ImageChargeScatterController.tpp
)

set (CORRECTION_HDRS
    Corrections/BoundaryCorrection.hpp
    Corrections/ImageChargeCorrection.hpp
    Corrections/ShiftedGreensCorrection.hpp
)

set (DRIVER_HDRS
    Drivers/BinnedLorentzSelfFieldDriver.hpp
    Drivers/BinnedLorentzSelfFieldDriver.tpp
    Drivers/MonolithicSelfFieldDriver.hpp
    Drivers/MonolithicSelfFieldDriver.tpp
    Drivers/SelfFieldDriver.hpp
)

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge")
install (FILES ${CORRECTION_HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge/Corrections")
install (FILES ${DRIVER_HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge/Drivers")

add_subdirectory (Binning)
add_subdirectory (Diagnostics)
add_subdirectory (Solvers)
