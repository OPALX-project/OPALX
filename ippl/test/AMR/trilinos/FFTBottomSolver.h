#ifndef FFT_SOLVER_H
#define FFT_SOLVER_H

#include "AmrMultiGridCore.h"

#include "Solver/FFTPoissonSolver.h"

class FFTBottomSolver : public BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                            Teuchos::RCP<amr::multivector_t> >,
                        public FFTPoissonSolver
{
public:
    
    FFTBottomSolver() { };
    
    void solver(const Teuchos::RCP<mv_t>& x,
                const Teuchos::RCP<mv_t>& b) { };
    
    void setOperator(const Teuchos::RCP<matrix_t>& A) { };
    
    
private:
    
    
};

#endif
