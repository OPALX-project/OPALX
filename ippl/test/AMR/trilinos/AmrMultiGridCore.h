#ifndef AMR_MULTI_GRID_CORE_H
#define AMR_MULTI_GRID_CORE

// boundary handlers
#include "AmrDirichletBoundary.h"
#include "AmrOpenBoundary.h"

// interpolaters
#include "AmrTrilinearInterpolater.h"
#include "AmrLagrangeInterpolater.h"

// linear solvers
#include "LinearSolver.h"
#include "BlockCGSolMgr.h"

// Trilinos headers
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
// #include <Teuchos_DefaultMpiComm.hpp> // wrapper for our communicator

namespace amr {
    typedef Epetra_CrsMatrix    matrix_t;
    typedef Epetra_Vector       vector_t;
    typedef Epetra_Map          dmap_t;
    typedef Epetra_MpiComm      comm_t;
}

#endif
