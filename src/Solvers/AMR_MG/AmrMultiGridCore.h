#ifndef AMR_MULTI_GRID_CORE_H
#define AMR_MULTI_GRID_CORE_H

// boundary handlers
#include "AmrDirichletBoundary.h"
#include "AmrOpenBoundary.h"
#include "AmrPeriodicBoundary.h"

// interpolaters
#include "AmrTrilinearInterpolater.h"
#include "AmrLagrangeInterpolater.h"
#include "AmrPCInterpolater.h"

// base level solvers
#include "BottomSolver.h"
#include "BelosBottomSolver.h"
#include "Amesos2BottomSolver.h"
#include "MueLuBottomSolver.h"

#include "AmrSmoother.h"

// preconditioners
#include "Ifpack2Preconditioner.h"
#include "MueLuPreconditioner.h"

#endif
