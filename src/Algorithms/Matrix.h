//
// This header file provides matrix and vector operations to replace boost::numeric::ublas.
// Kokkos-friendly implementation with fixed-size matrices for GPU/CPU portability.
//
// Copyright (c) 2023, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef OPAL_MATRIX_HH
#define OPAL_MATRIX_HH

#include <Kokkos_Core.hpp>
#include <cassert>
#include <vector>
#include "Ippl.h"

// Kokkos-friendly fixed-size matrix structure
template <int Rows, int Cols>
struct matrix_t {
    double m[Rows][Cols];

    // Default constructor - initializes to identity if square, zero otherwise
    KOKKOS_INLINE_FUNCTION
    matrix_t() {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                m[i][j] = 0.0;
            }
        }
        // Set identity for square matrices
        if (Rows == Cols) {
            for (int i = 0; i < Rows; ++i) {
                m[i][i] = 1.0;
            }
        }
    }

    // Constructor with initial value
    KOKKOS_INLINE_FUNCTION
    matrix_t(double init_value) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                m[i][j] = init_value;
            }
        }
    }

    // Copy constructor
    KOKKOS_INLINE_FUNCTION
    matrix_t(const matrix_t& right) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                m[i][j] = right(i, j);
            }
        }
    }

    // Assignment operator
    KOKKOS_INLINE_FUNCTION
    matrix_t& operator=(const matrix_t& right) {
        if (this != &right) {
            for (int i = 0; i < Rows; ++i) {
                for (int j = 0; j < Cols; ++j) {
                    m[i][j] = right(i, j);
                }
            }
        }
        return *this;
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

    KOKKOS_INLINE_FUNCTION
    int size1() const { return Rows; }

    KOKKOS_INLINE_FUNCTION
    int size2() const { return Cols; }
    
};

// Type aliases for common matrix sizes
using matrix3x3_t = matrix_t<3, 3>;
using matrix6x6_t = matrix_t<6, 6>;

// Kokkos-friendly vector_t using ippl::Vector
// Note: ippl::Vector uses unsigned int for dimension, so we use unsigned here too
template <unsigned Size>
using vector_t = ippl::Vector<double, Size>;

// Kokkos-friendly matrix-vector product for fixed-size matrices
template <int Rows, int Cols, class T>
KOKKOS_INLINE_FUNCTION
T prod_vector(const matrix_t<Rows, Cols>& rotation, const T& vect) {
    T prodVector(0.0);
    
    for (int i = 0; i < Rows; ++i) {
        double val = 0.0;
        for (int j = 0; j < Cols; ++j) {
            val += rotation(i, j) * vect[j];
        }
        prodVector[i] = val;
    }
    
    return prodVector;
}

// Kokkos-friendly matrix^T-vector product for fixed-size matrices
template <int Rows, int Cols, class T>
KOKKOS_INLINE_FUNCTION
T prod_vector_transpose(const matrix_t<Rows, Cols>& rotation, const T& vect) {
    T prodVector(0.0);
    
    for (int i = 0; i < Cols; ++i) {
        double val = 0.0;
        for (int j = 0; j < Rows; ++j) {
            val += rotation(j, i) * vect[j];
        }
        prodVector[i] = val;
    }
    
    return prodVector;
}

// Kokkos-friendly matrix transpose for fixed-size matrices
template <int Rows, int Cols>
KOKKOS_INLINE_FUNCTION
matrix_t<Cols, Rows> get_transpose(const matrix_t<Rows, Cols>& rotation) {
    matrix_t<Cols, Rows> transpose;
    for (int i = 0; i < Rows; ++i) {
        for (int j = 0; j < Cols; ++j) {
            transpose(j, i) = rotation(i, j);
        }
    }
    return transpose;
}

// Kokkos-friendly matrix-matrix multiplication
template <int Rows1, int Cols1, int Rows2, int Cols2>
KOKKOS_INLINE_FUNCTION
matrix_t<Rows1, Cols2> prod(const matrix_t<Rows1, Cols1>& a, const matrix_t<Rows2, Cols2>& b) {
    matrix_t<Rows1, Cols2> result(0.0);
    for (int i = 0; i < Rows1; ++i) {
        for (int j = 0; j < Cols2; ++j) {
            double sum = 0.0;
            for (int k = 0; k < Cols1; ++k) {
                sum += a(i, k) * b(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

// Kokkos-friendly matrix-vector multiplication (matrix * vector)
// Note: Casts int template parameters to unsigned for ippl::Vector compatibility
template <int Rows, int Cols>
KOKKOS_INLINE_FUNCTION
ippl::Vector<double, static_cast<unsigned>(Rows)> prod_matrix_vector(const matrix_t<Rows, Cols>& m, const ippl::Vector<double, static_cast<unsigned>(Cols)>& v) {
    ippl::Vector<double, static_cast<unsigned>(Rows)> result(0.0);
    for (int i = 0; i < Rows; ++i) {
        double sum = 0.0;
        for (int j = 0; j < Cols; ++j) {
            sum += m(i, j) * v[j];
        }
        result[i] = sum;
    }
    return result;
}

#endif
