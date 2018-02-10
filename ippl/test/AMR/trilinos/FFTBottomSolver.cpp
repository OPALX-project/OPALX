#include "FFTBottomSolver.h"


FFTBottomSolver::FFTBottomSolver(Mesh_t *mesh,
                                 FieldLayout_t *fl,
                                 std::string greensFunction,
                                 std::string bcz)
    : FFTPoissonSolver(mesh, fl, greensFunction, bcz)
{ }


void FFTBottomSolver::solve(const Teuchos::RCP<mv_t>& x,
                            const Teuchos::RCP<mv_t>& b )
{
    
    
}


void FFTBottomSolver::setOperator(const Teuchos::RCP<matrix_t>& A)
{
    // do nothing here
};
