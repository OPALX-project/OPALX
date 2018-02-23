#include "MueLuBottomSolver.h"


MueLuBottomSolver::MueLuBottomSolver()
    : conv_m(1.0e-4)
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
    // InitialGuessIsZero = false
    // startLevel = 0
    hierarchy_mp->Iterate(*b, *x, conv_m, false, 0);
}


void MueLuBottomSolver::setOperator(const Teuchos::RCP<matrix_t>& A) {
    
    finest_mp->Set("A", A);
    Teuchos::RCP<mv_t> nullspace = MultiVectorFactory::Build(map, 1);
    nullspace->putScalar(one);
    finest_mp->Set("Nullspace", nullspace);
//     finest_mp->Set("Coordinates", coordinates)
}


std::size_t MueLuBottomSolver::getNumIters() {
    return 1;
}
