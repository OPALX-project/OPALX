#ifndef BOTTOM_SOLVER_H
#define BOTTOM_SOLVER_H

#include "AmrMultiGridDefs.h"

/// Abstract base class for all base level solvers
template <class MatrixType, class VectorType>
class BottomSolver {
    
public:
    /*!
     * Solves
     * \f[
     *      Ax = b
     * \f]
     * @param x left-hand side
     * @param b right-hand side
     */
    virtual void solve(const VectorType& x,
                       const VectorType& b) = 0;
    
    /*!
     * Set the system matrix
     * @param A system matrix
     */
    virtual void setOperator(const MatrixType& A) = 0;
    
    
    /*!
     * @returns the number of required iterations
     */
    virtual std::size_t getNumIters() = 0;
};

#endif
