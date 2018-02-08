#include "AmesosBottomSolver.h"

AmesosBottomSolver::AmesosBottomSolver(std::string solvertype)
    : solvertype_m(solvertype)
{ }


AmesosBottomSolver::~AmesosBottomSolver() {
    solver_mp = Teuchos::null;
}


void AmesosBottomSolver::solve(const Teuchos::RCP<mv_t>& x,
                               const Teuchos::RCP<mv_t>& b)
{
    /*
     * solve linear system Ax = b
     */
    solver_mp->solve(x.get(), b.get());
}


void AmesosBottomSolver::setOperator(const Teuchos::RCP<matrix_t>& A) {
    try {
        solver_mp = Amesos2::create<matrix_t, mv_t>(solvertype_m, A);
    } catch(const std::invalid_argument& ex) {
        //TODO change to gmsg IPPL when built-in in OPAL
        std::cerr << ex.what() << std::endl;
    }
    
    solver_mp->symbolicFactorization();
    solver_mp->numericFactorization();
}


std::size_t AmesosBottomSolver::getNumIters() {
    return 1;   // direct solvers do only one step
}
