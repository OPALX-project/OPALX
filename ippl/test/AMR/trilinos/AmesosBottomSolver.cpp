#include "AmesosBottomSolver.h"

AmesosBottomSolver::AmesosBottomSolver(std::string solvertype)
    : solvertype_m(solvertype)
{ }


AmesosBottomSolver::~AmesosBottomSolver() {
    solver_mp = Teuchos::null;
}


void AmesosBottomSolver::solve(const Teuchos::RCP<vector_t>& x,
                               const Teuchos::RCP<vector_t>& b)
{
    /*
     * solve linear system Ax = b
     */
    solver_mp->solve(x.get(), b.get());
}


void AmesosBottomSolver::setOperator(const Teuchos::RCP<matrix_t>& A) {
    solver_mp = Amesos2::create<matrix_t, vector_t>(solvertype_m, A,
                                                    Teuchos::RCP<vector_t>(),
                                                    Teuchos::RCP<vector_t>());
}
