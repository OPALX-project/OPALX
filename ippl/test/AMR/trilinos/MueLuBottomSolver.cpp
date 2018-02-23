#include "MueLuBottomSolver.h"

MueLuBottomSolver::MueLuBottomSolver()
    : A_mp(Teuchos::null),
      tolerance_m(1.0e-4)
{
    hierarchy_mp = rcp(new hierarchy_t());
    
    // WCYCLE also supported
    hierarchy_mp->SetCycle(MueLu::CycleType::VCYCLE);
    hierarchy_mp->SetMaxCoarseSize(200);
    hierarchy_mp->SetPRrebalance(false);
    hierarchy_mp->SetImplicitTranspose(true);
    hierarchy_mp->IsPreconditioner(false);
    
    hierarchy_mp->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    
    finest_mp = hierarchy_mp->GetLevel();
    
    finest_mp->setDefaultVerbLevel(Teuchos::VERB_HIGH);
}


void MueLuBottomSolver::solve(const Teuchos::RCP<mv_t>& x,
                              const Teuchos::RCP<mv_t>& b)
{
    // MueLu requires Xpetra multivectors (wrap them)
    Teuchos::RCP<xmv_t> xx = MueLu::TpetraMultiVector_To_XpetraMultiVector(x);
    Teuchos::RCP<xmv_t> xb = MueLu::TpetraMultiVector_To_XpetraMultiVector(b);

    // InitialGuessIsZero = false
    // startLevel = 0
    MueLu::ReturnType rtype = hierarchy_mp->Iterate(*xb, *xx, tolerance_m, false, 0);
    if ( rtype != MueLu::ReturnType::Converged )
	std::cout << "not converged" << std::endl;

    // put multivector back
    x->assign(*util_t::MV2NonConstTpetraMV2(*xx));
}


void MueLuBottomSolver::setOperator(const Teuchos::RCP<matrix_t>& A) {
    
    A_mp = MueLu::TpetraCrs_To_XpetraMatrix<scalar_t, lo_t, go_t, node_t>(A);
    A_mp->SetFixedBlockSize(1); // only 1 DOF per node (pure Laplace problem)
    finest_mp->Set("A", A_mp);
    
    Teuchos::RCP<mv_t> nullspace = Teuchos::rcp(new mv_t(A->getRowMap(), 1));
    Teuchos::RCP<xmv_t> xnullspace = MueLu::TpetraMultiVector_To_XpetraMultiVector(nullspace);
    nullspace->putScalar(1.0);
    finest_mp->Set("Nullspace", xnullspace);
    //finest_mp->Set("Coordinates", coordinates)
}


std::size_t MueLuBottomSolver::getNumIters() {
    return 1;
}
