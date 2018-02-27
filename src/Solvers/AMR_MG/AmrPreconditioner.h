#ifndef AMR_PRECONDITIONER_H
#define AMR_PRECONDITIONER_H

#include "AmrMultiGridDefs.h"

namespace amr {
    enum Preconditioner {
        NONE,
        ILUT,           ///< incomplete LU
        CHEBYSHEV,
        RILUK,          ///< ILU(k)
        SA,             ///< smoothed aggregation multigrid
        JACOBI,         ///< Jacobi point relaxation
        BLOCK_JACOBI,   ///< Jacobi block relaxation
        GS,             ///< Gauss-Seidel point relaxation
        BLOCK_GS        ///< Gauss-Seidel block relaxation
    };
}

/// Bottom solver preconditioners
template <class Matrix, class Level>
class AmrPreconditioner
{
public:
    typedef amr::operator_t operator_t;
    
public:
    
    /*!
     * Instantiate the preconditioner matrix
     * @param A matrix for which to create preconditioner
     * @param level_p bottom level if necessary to build preconditioner
     */
    virtual void create(const Teuchos::RCP<Matrix>& A, Level* level_p = nullptr) = 0;
    
    /*!
     * @returns the preconditioner
     */
    virtual Teuchos::RCP<operator_t> get() = 0;
};


#endif
