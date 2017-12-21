#ifndef BOTTOM_SOLVER_H
#define BOTTOM_SOLVER_H

namespace amr {
    enum Preconditioner {
        ILUT,       // incomplete LU
        CHEBYSHEV,
        NONE
    };
}

template <class MatrixType, class VectorType>
class BottomSolver {
    
public:
    virtual void solve(const VectorType& x,
                       const VectorType& b) = 0;
    
    virtual void setOperator(const MatrixType& A) = 0;
};

#endif
