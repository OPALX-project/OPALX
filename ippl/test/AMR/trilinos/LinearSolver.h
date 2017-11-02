#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

template <class MatrixType, class VectorType>
class LinearSolver {
    
public:
    virtual void solve(const MatrixType& A,
                       VectorType& x,
                       const VectorType& b, double tol) = 0;
    
    virtual void residual(VectorType& r, const VectorType& x, const VectorType& b) = 0;
    
};

#endif
