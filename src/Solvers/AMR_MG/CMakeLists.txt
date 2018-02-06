set (_SRCS
    AmesosBottomSolver.cpp
    AmrMultiGrid.cpp
    AmrSmoother.cpp
    BelosBottomSolver.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_opal_sources(${_SRCS})

set (HDRS
    AmesosBottomSolver.h
    AmrBoundary.h
    AmrDirichletBoundary.h
    AmrInterpolater.h
    AmrLagrangeInterpolater.h
    AmrLagrangeInterpolater.hpp
    AmrMultiGridCore.h
    AmrMultiGrid.h
    AmrMultiGridLevel.h
    AmrMultiGridLevel.hpp
    AmrOpenBoundary.h
    AmrPCInterpolater.h
    AmrPCInterpolater.hpp
    AmrPeriodicBoundary.h
    AmrSmoother.h
    AmrTrilinearInterpolater.h
    AmrTrilinearInterpolater.hpp
    BelosBottomSolver.h
    BottomSolver.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Solvers/AMR_MG/")