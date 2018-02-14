#ifndef AMR_MULTI_GRID_DEFS_H
#define AMR_MULTI_GRID_DEFS_H

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
    
#ifdef AMR_MG_SERIAL_NODE
    typedef ::Kokkos::Compat::KokkosSerialWrapperNode node_t;
#elif AMR_MG_PTHREAD_NODE
    typedef ::Kokkos::Compat::KokkosThreadsWrapperNode node_t;
#elif AMR_MG_OPENMP_NODE
    typedef ::Kokkos::Compat::KokkosOpenMPWrapperNode node_t;
#elif AMR_MG_CUDA_NODE
    typedef ::Kokkos::Compat::KokkosCudaWrapperNode node_t;
#else
    typedef KokkosClassic::DefaultNode::DefaultNodeType node_t;
#endif
    
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
