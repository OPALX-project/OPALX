#include "BelosBottomSolver.h"

#include "Utilities/OpalException.h"

BelosBottomSolver::BelosBottomSolver(std::string solvertype, Preconditioner precond)
    : problem_mp( Teuchos::rcp( new problem_t() ) ),
      precond_m(precond),
      prectype_m(""),
      P_mp(Teuchos::null),
      reltol_m(1.0e-9),
      maxiter_m(1000)
{
    this->initSolver_m(solvertype);
    this->initPreconditioner_m();
}


BelosBottomSolver::~BelosBottomSolver() {
    problem_mp = Teuchos::null;
    params_mp = Teuchos::null;
    solver_mp = Teuchos::null;
    prec_mp = Teuchos::null;
    P_mp = Teuchos::null;
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
    
    if ( precond_m != Preconditioner::NONE ) {
        this->computePrecond_m(A);
        problem_mp->setLeftPrec(P_mp);
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
    
    params_mp->set("Convergence Tolerance", reltol_m);
    params_mp->set("Maximum Iterations", maxiter_m);
    
    solver_mp = factory.create(solvertype, params_mp);
}


void BelosBottomSolver::initPreconditioner_m()
{
    prec_mp = Teuchos::parameterList();
    
    switch ( precond_m ) {
        case Preconditioner::ILUT:
            // inclomplete LU
            prectype_m = "ILUT";
            
            prec_mp->set("fact: ilut level-of-fill", 5.0);
            prec_mp->set("fact: drop tolerance", 0.0);
            
            break;
        case Preconditioner::CHEBYSHEV:
            prectype_m = "CHEBYSHEV";
            prec_mp->set("chebyshev: degree", 1);
            break;
        case Preconditioner::RILUK:
            prectype_m = "RILUK";
            prec_mp->set("fact: iluk level-of-fill", 0);
            prec_mp->set("fact: relax value", 0.0);
            prec_mp->set("fact: absolute threshold", 0.0);
            prec_mp->set("fact: relative threshold", 1.0);
            break;
        case Preconditioner::NONE:
            prectype_m = "";
            break;
        default:
            throw OpalException("BelosBottomSolver::initPreconditioner_m()",
                                "This type of Ifpack2 preconditioner not supported.");
    }
    
}


void BelosBottomSolver::computePrecond_m(const Teuchos::RCP<const matrix_t>& A)
{
    Ifpack2::Factory factory;
    
    P_mp = factory.create(prectype_m, A);
    P_mp->setParameters(*prec_mp);
    P_mp->initialize();
    P_mp->compute();
}
