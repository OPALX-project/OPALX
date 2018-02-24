#ifndef MUELU_BOTTOM_SOLVER_H
#define MUELU_BOTTOM_SOLVER_H

#include "BottomSolver.h"

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_Utilities.hpp>

#include "MueLu_Factory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_RepartitionHeuristicFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_Zoltan2Interface.hpp"
#include "MueLu_RebalanceAcFactory.hpp"

class MueLuBottomSolver : public BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                              Teuchos::RCP<amr::multivector_t> >
{
public:
    typedef amr::matrix_t matrix_t;
    typedef amr::vector_t vector_t;
    typedef amr::scalar_t scalar_t;
    typedef amr::multivector_t mv_t;
    typedef amr::operator_t op_t;
    typedef amr::local_ordinal_t lo_t;
    typedef amr::global_ordinal_t go_t;
    typedef amr::node_t node_t;
    
    typedef MueLu::Hierarchy<scalar_t, lo_t, go_t, node_t> hierarchy_t;
    typedef MueLu::Level level_t;
    typedef Xpetra::Matrix<scalar_t, lo_t, go_t, node_t> xmatrix_t;
    typedef Xpetra::MultiVector<scalar_t, lo_t, go_t, node_t> xmv_t;
    typedef MueLu::Utilities<scalar_t, lo_t, go_t, node_t> util_t;
    
    typedef MueLu::TentativePFactory<scalar_t, lo_t, go_t, node_t> tentPFactory_t;
    typedef MueLu::SaPFactory<scalar_t, lo_t, go_t, node_t> saPFactory_t;
    typedef MueLu::TransPFactory<scalar_t, lo_t, go_t, node_t> transPFactory_t;
    typedef MueLu::RAPFactory<scalar_t, lo_t, go_t, node_t> rAPFactory_t;
    typedef MueLu::FactoryManager<scalar_t, lo_t, go_t, node_t> fManager_t;
    
    typedef MueLu::CoordinatesTransferFactory<scalar_t, lo_t, go_t, node_t> coordTransferFactory_t;
    typedef MueLu::RepartitionHeuristicFactory<scalar_t, lo_t, go_t, node_t> repartFactory_t;
    typedef MueLu::Zoltan2Interface<scalar_t, lo_t, go_t, node_t> zoltan_t;
    typedef MueLu::RebalanceTransferFactory<scalar_t, lo_t, go_t, node_t> rebalTransferFactory_t;
    typedef MueLu::RebalanceAcFactory<scalar_t, lo_t, go_t, node_t> rebalAcFactory_t;
    typedef MueLu::SmootherPrototype<scalar_t, lo_t, go_t, node_t> smootherPrototype_t;
    typedef MueLu::SmootherFactory<scalar_t, lo_t, go_t, node_t> smootherFactory_t;
    typedef MueLu::TrilinosSmoother<scalar_t, lo_t, go_t, node_t> smoother_t;
    typedef MueLu::DirectSolver<scalar_t, lo_t, go_t, node_t> solver_t;
    
public:
    
    MueLuBottomSolver();
    
    void solve(const Teuchos::RCP<mv_t>& x,
               const Teuchos::RCP<mv_t>& b);
    
    void setOperator(const Teuchos::RCP<matrix_t>& A);
    
    std::size_t getNumIters();

private:
    void setupAggregation_m();
    void setupSmoother_m();
    void setupCoarser_m();    
    
private:
    Teuchos::RCP<hierarchy_t> hierarchy_mp;     ///< manages the multigrid hierarchy
    
    Teuchos::RCP<level_t> finest_mp;            ///< finest level of hierarchy

    Teuchos::RCP<xmatrix_t> A_mp;               ///< MueLu requires Xpetra

    scalar_t tolerance_m;                       ///< stopping criteria of multigrid iteration
    
};

#endif
