#include "BelosBottomSolver.h"

#include "Utilities/OpalException.h"

BelosBottomSolver::BelosBottomSolver(std::string solvertype,
                                     const std::shared_ptr<prec_t>& prec_p)
    : problem_mp( Teuchos::rcp( new problem_t() ) ),
      prec_mp(prec_p),
      reltol_m(1.0e-9),
      maxiter_m(1000)
{
    this->initSolver_m(solvertype);
}


BelosBottomSolver::~BelosBottomSolver() {
    problem_mp = Teuchos::null;
    params_mp = Teuchos::null;
    solver_mp = Teuchos::null;
}


void BelosBottomSolver::solve(const Teuchos::RCP<mv_t>& x,
                              const Teuchos::RCP<mv_t>& b)
{
    /*
     * solve linear system Ax = b
     */
    
    // change sign of rhs due to change of A in BelosBottomSolver::setOperator
    b->scale(-1.0);
    
    problem_mp->setProblem(x, b);
    
    solver_mp->setProblem(problem_mp);
    
    Belos::ReturnType ret = solver_mp->solve();

/*
    Teuchos::Array<typename Teuchos::ScalarTraits<scalar_t>::magnitudeType> normVec(1);
    mv_t Ax(b->getMap(),1);
    mv_t residual(b->getMap(),1);
    problem_mp->getOperator()->apply(*x, residual);
    residual.update(1.0, *b, -1.0);
    residual.norm2(normVec);

    std::cout << normVec[0] << std::endl;
*/  
    if ( ret != Belos::Converged ) {
        //TODO When put into OPAL: Replace with gmsg
        std::cerr << "Warning: Bottom solver not converged. Achieved tolerance"
                  << " after " << solver_mp->getNumIters() << " iterations is "
                  << solver_mp->achievedTol() << "." << std::endl;
    }
    
    // undo sign change
    b->scale(-1.0);
}


void BelosBottomSolver::setOperator(const Teuchos::RCP<matrix_t>& A) {
    
    // make positive definite --> rhs has to change sign as well
    A->resumeFill();
    A->scale(-1.0);
    A->fillComplete();
    
    if ( problem_mp == Teuchos::null )
        throw OpalException("BelosBottomSolver::setOperator()",
                            "No problem defined.");
    
    problem_mp->setOperator(A);
    
    if ( prec_mp != nullptr ) {
	prec_mp->create(A);
        problem_mp->setLeftPrec(prec_mp->get());
    }
}


void BelosBottomSolver::initSolver_m(std::string solvertype) {
    
    Belos::SolverFactory<scalar_t, mv_t, op_t> factory;
    
    params_mp = Teuchos::parameterList();
    
    params_mp->set("Block Size", 1);
    params_mp->set("Convergence Tolerance", reltol_m);
    params_mp->set("Adaptive Block Size", false);
    params_mp->set("Use Single Reduction", true);
    params_mp->set("Explicit Residual Scaling", "Norm of RHS");
    params_mp->set("Maximum Iterations", maxiter_m);
    params_mp->set("Verbosity", Belos::Errors + Belos::Warnings +
                                Belos::TimingDetails + Belos::FinalSummary +
                                Belos::StatusTestDetails);
    params_mp->set("Output Frequency", 10);
    
    solver_mp = factory.create(solvertype, params_mp);
}
