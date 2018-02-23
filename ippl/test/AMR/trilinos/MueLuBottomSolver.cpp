#include "MueLuBottomSolver.h"

#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_RAPFactory.hpp"

MueLuBottomSolver::MueLuBottomSolver()
    : hierarchy_mp(Teuchos::null),
      finest_mp(Teuchos::null),
      A_mp(Teuchos::null),
      tolerance_m(1.0e-4)
{
/*
    hierarchy_mp = Teuchos::rcp(new hierarchy_t());
    
    // WCYCLE also supported
    hierarchy_mp->SetCycle(MueLu::CycleType::VCYCLE);
    hierarchy_mp->SetMaxCoarseSize(200);
    hierarchy_mp->SetPRrebalance(false);
    hierarchy_mp->SetImplicitTranspose(true);
    hierarchy_mp->IsPreconditioner(false);
    
    hierarchy_mp->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    
    finest_mp = hierarchy_mp->GetLevel();
    
    finest_mp->setDefaultVerbLevel(Teuchos::VERB_HIGH);
*/
}


void MueLuBottomSolver::solve(const Teuchos::RCP<mv_t>& x,
                              const Teuchos::RCP<mv_t>& b)
{
    // MueLu requires Xpetra multivectors (wrap them)
    Teuchos::RCP<xmv_t> xx = MueLu::TpetraMultiVector_To_XpetraMultiVector(x);
    Teuchos::RCP<xmv_t> xb = MueLu::TpetraMultiVector_To_XpetraMultiVector(b);

    // InitialGuessIsZero = false
    // startLevel = 0
    hierarchy_mp->Iterate(*xb, *xx, 20/*tolerance_m*/, false, 0);
    
    // put multivector back
    x->assign(*util_t::MV2NonConstTpetraMV2(*xx));
}


void MueLuBottomSolver::setOperator(const Teuchos::RCP<matrix_t>& A) {
    Teuchos::RCP<MueLu::TentativePFactory<scalar_t, lo_t, go_t, node_t> > TentativePFact = Teuchos::rcp( new MueLu::TentativePFactory<scalar_t, lo_t, go_t, node_t>() );
    Teuchos::RCP<MueLu::SaPFactory<scalar_t, lo_t, go_t, node_t> > SaPFact = Teuchos::rcp( new MueLu::SaPFactory<scalar_t, lo_t, go_t, node_t>() );
    Teuchos::RCP<MueLu::TransPFactory<scalar_t, lo_t, go_t, node_t> >     RFact          = Teuchos::rcp( new MueLu::TransPFactory<scalar_t, lo_t, go_t, node_t>());

    MueLu::FactoryManager<scalar_t, lo_t, go_t, node_t> M;
    M.SetFactory("Ptent", TentativePFact);
    M.SetFactory("P",     SaPFact);
    M.SetFactory("R",     RFact);

    M.SetFactory("Smoother", Teuchos::null);      //skips smoother setup
    M.SetFactory("CoarseSolver", Teuchos::null);  //skips coarsest solve setup

    int startLevel = 0;
    int maxLevels = 10;


    // Indicate which Hierarchy operators we want to keep
    hierarchy_mp = Teuchos::rcp(new hierarchy_t());  
    hierarchy_mp->Keep("P", SaPFact.get());  //SaPFact is the generating factory for P.
    hierarchy_mp->Keep("R", RFact.get());    //RFact is the generating factory for R.
    hierarchy_mp->Keep("Ptent", TentativePFact.get());  //SaPFact is the generating factory for P.

    hierarchy_mp->SetCycle(MueLu::CycleType::VCYCLE);
    hierarchy_mp->SetMaxCoarseSize(200);
    hierarchy_mp->SetPRrebalance(false);
    hierarchy_mp->SetImplicitTranspose(true);
    hierarchy_mp->IsPreconditioner(false);
    hierarchy_mp->setDefaultVerbLevel(Teuchos::VERB_HIGH);

    // Create Gauss-Seidel smoother.                                                                                  
    std::string ifpackType = "RELAXATION";
    Teuchos::ParameterList ifpackList;
    ifpackList.set("relaxation: sweeps", (lo_t) 4);
    ifpackList.set("relaxation: damping factor", (scalar_t) 1.0);
    Teuchos::RCP<MueLu::SmootherPrototype<scalar_t, lo_t, go_t, node_t> > smootherPrototype = Teuchos::rcp(new MueLu::TrilinosSmoother<scalar_t, lo_t, go_t, node_t>(ifpackType, ifpackList));

    M.SetFactory("Smoother", Teuchos::rcp(new MueLu::SmootherFactory<scalar_t, lo_t, go_t, node_t>(smootherPrototype)));

    // Create coarsest solver.
    Teuchos::RCP<MueLu::SmootherPrototype<scalar_t, lo_t, go_t, node_t> > coarseSolverPrototype = Teuchos::rcp( new MueLu::DirectSolver<scalar_t, lo_t, go_t, node_t>() );
    Teuchos::RCP<MueLu::SmootherFactory<scalar_t, lo_t, go_t, node_t> >   coarseSolverFact = Teuchos::rcp( new MueLu::SmootherFactory<scalar_t, lo_t, go_t, node_t>(coarseSolverPrototype, Teuchos::null) );
    M.SetFactory("CoarseSolver", coarseSolverFact);

    A_mp = MueLu::TpetraCrs_To_XpetraMatrix<scalar_t, lo_t, go_t, node_t>(A);
    A_mp->SetFixedBlockSize(1); // only 1 DOF per node (pure Laplace problem)

    finest_mp = hierarchy_mp->GetLevel();
    finest_mp->Set("A", A_mp);
    
    Teuchos::RCP<mv_t> nullspace = Teuchos::rcp(new mv_t(A->getRowMap(), 1));
    Teuchos::RCP<xmv_t> xnullspace = MueLu::TpetraMultiVector_To_XpetraMultiVector(nullspace);
    nullspace->putScalar(1.0);
    finest_mp->Set("Nullspace", xnullspace);
    //finest_mp->Set("Coordinates", coordinates)

    hierarchy_mp->Setup(M,startLevel,maxLevels);
}


std::size_t MueLuBottomSolver::getNumIters() {
    return 1;
}


void MueLuBottomSolver::setupAggregation_m() {

}


void MueLuBottomSolver::setupSmoother_m() {


}


void MueLuBottomSolver::setupCoarser_m() {

}
