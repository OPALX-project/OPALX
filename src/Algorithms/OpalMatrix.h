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

    // copy constructor
    KOKKOS_INLINE_FUNCTION
    Matrix_t(const Matrix_t& right){
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                m[i][j] = right(i,j); 
            }
        }
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
template <class T>
KOKKOS_INLINE_FUNCTION
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

KOKKOS_INLINE_FUNCTION
Matrix_t get_transpose(const Matrix_t& rotation){
    Matrix_t transpose;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            transpose(j,i) = rotation(i,j);
        }
    }
    return transpose;
}

#endif
