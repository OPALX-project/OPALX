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
    Field_t lhs;
    Vector_t hr;
    
    this->vector2field_m(lhs, b);
    
    FFTPoissonSolver::computePotential(lhs, hr);
    
    this->field2vector_m(lhs, x);
}


void FFTBottomSolver::setOperator(const Teuchos::RCP<matrix_t>& A)
{
    // do nothing here
};


void FFTBottomSolver::field2vector_m(const Field_t& field,
                                     const Teuchos::RCP<mv_t>& vector)
{
    
}
    

void FFTBottomSolver::vector2field_m(Field_t& field,
                                     const Teuchos::RCP<mv_t>& vector)
{
    
    
}
