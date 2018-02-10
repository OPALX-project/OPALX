#ifndef FFT_SOLVER_H
#define FFT_SOLVER_H

#include "AmrMultiGridCore.h"

#include "Solvers/FFTPoissonSolver.h"

class FFTBottomSolver : public BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                            Teuchos::RCP<amr::multivector_t> >,
                        public FFTPoissonSolver
{
public:
    typedef amr::matrix_t matrix_t;
    typedef amr::vector_t vector_t;
    typedef amr::scalar_t scalar_t;
    typedef amr::multivector_t mv_t;
    typedef amr::operator_t op_t;
    typedef amr::local_ordinal_t lo_t;
    typedef amr::global_ordinal_t go_t;
    
public:
    
    FFTBottomSolver(Mesh_t *mesh,
                    FieldLayout_t *fl,
                    std::string greensFunction,
                    std::string bcz);
    
    void solve(const Teuchos::RCP<mv_t>& x,
               const Teuchos::RCP<mv_t>& b);
    
    void setOperator(const Teuchos::RCP<matrix_t>& A);
    
    
private:
    void field2vector_m(const Field_t& field,
                        const Teuchos::RCP<mv_t>& vector);
    
    void vector2field_m(Field_t& field,
                        const Teuchos::RCP<mv_t>& vector);
};

#endif
