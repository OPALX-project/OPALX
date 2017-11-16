#ifndef AMR_MULTI_GRID_CORE_H
#define AMR_MULTI_GRID_CORE

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

#include "AmrSmoother.h"

// Trilinos headers
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DefaultMpiComm.hpp> // wrapper for our communicator

#include <Kokkos_DefaultNode.hpp>

namespace amr {
    // All Tpetra
    typedef double scalar_t;
    typedef int local_ordinal_t;
    typedef int global_ordinal_t;
    typedef KokkosClassic::DefaultNode::DefaultNodeType node_t;
    
    typedef Tpetra::CrsMatrix<scalar_t,
                              local_ordinal_t,
                              global_ordinal_t,
                              node_t
            > matrix_t;
            
    typedef Tpetra::Vector<scalar_t,
                           local_ordinal_t,
                           global_ordinal_t,
                           node_t
            > vector_t;
    
    typedef Tpetra::Operator<scalar_t,
                             local_ordinal_t,
                             global_ordinal_t,
                             node_t
            > operator_t;
    
    typedef Tpetra::MultiVector<scalar_t,
                                local_ordinal_t,
                                global_ordinal_t,
                                node_t
            > multivector_t;
    
    
    typedef Tpetra::Map<local_ordinal_t,
                        global_ordinal_t,
                        node_t
            > dmap_t;
    
    typedef Teuchos::MpiComm<int>    comm_t;
}

#endif
