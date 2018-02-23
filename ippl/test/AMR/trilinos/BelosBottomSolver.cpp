#include "BelosBottomSolver.h"

#include "Ippl.h"
#include "Utilities/OpalException.h"

BelosBottomSolver::BelosBottomSolver(std::string solvertype,
                                     const std::shared_ptr<prec_t>& prec_p)
    : problem_mp( Teuchos::rcp( new problem_t() ) ),
      prec_mp(prec_p),
      reltol_m(1.0e-9),
      maxiter_m(100)
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
    
    static IpplTimings::TimerRef precTimer = IpplTimings::getTimer("AMR MG prec setup");

    if ( prec_mp != nullptr ) {
        IpplTimings::startTimer(precTimer);
        prec_mp->create(A);
        IpplTimings::stopTimer(precTimer);
        problem_mp->setLeftPrec(prec_mp->get());
    }
}


std::size_t BelosBottomSolver::getNumIters() {
    if ( solver_mp == Teuchos::null )
        throw OpalException("BelosBottomSolver::getNumIters()",
                            "No solver initialized.");
    
    return solver_mp->getNumIters();
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
    params_mp->set("Verbosity", Belos::Errors + Belos::Warnings);
    params_mp->set("Output Frequency", 1);
    
    solver_mp = factory.create(solvertype, params_mp);
}
