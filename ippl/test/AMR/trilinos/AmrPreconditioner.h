#ifndef AMR_PRECONDITIONER
#define AMR_PRECONDITIONER

#include "AmrMultiGridDefs.h"

namespace amr {
    enum Preconditioner {
        NONE,
        ILUT,       //incomplete LU
        CHEBYSHEV,
        RILUK,      // ILU(k)
        SA          // smoothed aggregation multigrid
    };
}

/// Bottom solver preconditioners
template <class MatrixType>
class AmrPreconditioner
{
public:
    typedef amr::operator_t operator_t;
    
public:
    
    /*!
     * Instantiate the preconditioner matrix
     * @param A matrix for which to create preconditioner
     */
    virtual void create(Teuchos::RCP<MatrixType>& A) = 0;
    
    virtual Teuchos::RCP<operator_t> get() = 0;
};


#endif
