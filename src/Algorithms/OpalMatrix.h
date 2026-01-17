#ifndef OPAL_MATRIX_HH
#define OPAL_MATRIX_HH

#include <Kokkos_Core.hpp>

/**
 * @brief Simple 3x3 Matrix class for Opalx
 * 
 * @note Implements the necessary oparations for coordinate system transformations.
 * @note For now we do not use the Kokkos BLAS API.
 */
struct OpalMatrix_t {
    double m[3][3];

    /**
     * @brief Default constructor - initializes to Identity Matrix
     */
    KOKKOS_INLINE_FUNCTION
    OpalMatrix_t() {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                m[i][j] = 0.0;
            }
        }
        m[0][0] = 1.0;
        m[1][1] = 1.0;
        m[2][2] = 1.0;
    }

    /**
     * @brief Copy constructor
     * @param right Matrix to copy
     */
    KOKKOS_INLINE_FUNCTION
    OpalMatrix_t(const OpalMatrix_t& right){
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                m[i][j] = right(i,j); 
            }
        }
    }

    /**
     * @brief Access operator
     * @param r Row index
     * @param c Column index
     * @returns Reference to element at (r,c)
     */
    KOKKOS_INLINE_FUNCTION
    double& operator()(int r, int c) {
        return m[r][c];
    }

    /**
     * @brief Const access operator
     * @param r Row index
     * @param c Column index
     * @returns Value of element at (r,c)
     */
    KOKKOS_INLINE_FUNCTION
    double operator()(int r, int c) const {
        return m[r][c];
    }
};

/**
 * @brief Takes the Matrix - Vector Product
 * @param rotation 3x3 Matrix
 * @param vect 3D Vector
 * @returns Product Vector
 */
template <class T>
KOKKOS_INLINE_FUNCTION
T prod_vector(const OpalMatrix_t& rotation, const T& vect) {
    
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

/**
 * @brief Takes the Matrix^T - Vector Product
 * @param rotation 3x3 Matrix
 * @param vect 3D Vector
 * @returns Product Vector 
 */
template <class T>
KOKKOS_INLINE_FUNCTION
T prod_vector_transpose(const OpalMatrix_t& rotation, const T& vect) {
    
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

/**
 * @brief Returns the transpose of a given OpalMatrix_t
 * @param rotation 3x3 Matrix
 * @returns Transpose Matrix
 */
KOKKOS_INLINE_FUNCTION
OpalMatrix_t get_transpose(const OpalMatrix_t& rotation){
    OpalMatrix_t transpose;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            transpose(j,i) = rotation(i,j);
        }
    }
    return transpose;
}

#endif
