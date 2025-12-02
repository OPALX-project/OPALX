// Implements 3x3 Matrix with the Matrix - Vector Product 
// as a KOKKOS_INLINE_FUNCTION

#ifndef OPAL_MATRIX_HH
#define OPAL_MATRIX_HH

#include <Kokkos_Core.hpp>

struct Matrix_t {
    double m[3][3];

    // Default constructor  
    KOKKOS_INLINE_FUNCTION
    Matrix_t() {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                m[i][j] = 0.0;
            }
        }
        m[0][0] = 1.0;
        m[1][1] = 1.0;
        m[2][2] = 1.0;
    }

    // Operators 
    KOKKOS_INLINE_FUNCTION
    double& operator()(int r, int c) {
        return m[r][c];
    }

    KOKKOS_INLINE_FUNCTION
    double operator()(int r, int c) const {
        return m[r][c];
    }
};

// Takes the Matrix-Vector Product
template <class T>
KOKKOS_INLINE_FUNCTION
T prod_vector(const Matrix_t& rotation, const T& vect) {
    
    T prodVector(0.0);

    for (int i = 0; i < 3; ++i) {
        double val = 0.0;
        for (int j = 0; j < 3; ++j) {
            val += rotation(i,j) * vect[j];
        }
        prodVector[i] = val;
    }

    return prodVector;
}

// Takes the Matrix^T - Vector Product
T prod_vector_transpose(const Matrix_t& rotation, const T& vect) {
    
    T prodVector(0.0);

    for (int i = 0; i < 3; ++i) {
        double val = 0.0;
        for (int j = 0; j < 3; ++j) {
            val += rotation(j,i) * vect[j];
        }
        prodVector[i] = val;
    }

    return prodVector;
}

/*
#include <boost/numeric/ublas/matrix.hpp>

typedef boost::numeric::ublas::matrix<double> matrix_t;

template <class T>
T prod_boost_vector(boost::numeric::ublas::matrix<double> rotation, const T& vect) {
    boost::numeric::ublas::vector<double> boostVector(3);

    boostVector(0) = vect[0];
    boostVector(1) = vect[1];
    boostVector(2) = vect[2];

    boostVector = boost::numeric::ublas::prod(rotation, boostVector);

    T prodVector(0.0);

    prodVector[0] = boostVector(0);
    prodVector[1] = boostVector(1);
    prodVector[2] = boostVector(2);

    return prodVector;
}
*/
#endif
