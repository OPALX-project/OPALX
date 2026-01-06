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

// Kokkos-friendly fixed-size matrix structure
template <int Rows, int Cols>
struct OpalMatrix_t {
    double m[Rows][Cols];

    // Default constructor - initializes to identity if square, zero otherwise
    KOKKOS_INLINE_FUNCTION
    OpalMatrix_t() {
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
    OpalMatrix_t(double init_value) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                m[i][j] = init_value;
            }
        }
    }

    // Copy constructor
    KOKKOS_INLINE_FUNCTION
    OpalMatrix_t(const OpalMatrix_t& right) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                m[i][j] = right(i, j);
            }
        }
    }

    // Assignment operator
    KOKKOS_INLINE_FUNCTION
    OpalMatrix_t& operator=(const OpalMatrix_t& right) {
        if (this != &right) {
            for (int i = 0; i < Rows; ++i) {
                for (int j = 0; j < Cols; ++j) {
                    m[i][j] = right(i, j);
                }
            }
        }
        return *this;
    }

    // Conversion constructor from matrix_t (non-Kokkos, host-only)
    // This will be defined after matrix_t is fully defined
    // Forward declaration of matrix_t is needed
    OpalMatrix_t(const class matrix_t& m);
    
    // Assignment from matrix_t (non-Kokkos, host-only)
    OpalMatrix_t& operator=(const class matrix_t& m);

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
using OpalMatrix3x3_t = OpalMatrix_t<3, 3>;
using OpalMatrix6x6_t = OpalMatrix_t<6, 6>;

// Forward declaration
class matrix_t;

// Legacy matrix_t class for backward compatibility
// Uses dynamic allocation for variable sizes
class matrix_t {
public:
    matrix_t() : rows_(0), cols_(0) {}
    
    matrix_t(size_t rows, size_t cols, double init_value = 0.0)
        : rows_(rows), cols_(cols), data_(rows * cols, init_value) {}
    
    double& operator()(size_t i, size_t j) {
        assert(i < rows_ && j < cols_);
        return data_[i * cols_ + j];
    }
    
    const double& operator()(size_t i, size_t j) const {
        assert(i < rows_ && j < cols_);
        return data_[i * cols_ + j];
    }
    
    size_t size1() const { return rows_; }
    size_t size2() const { return cols_; }
    
    void resize(size_t rows, size_t cols, double init_value = 0.0) {
        rows_ = rows;
        cols_ = cols;
        data_.resize(rows * cols, init_value);
    }

    // Conversion from OpalMatrix3x3_t
    matrix_t(const OpalMatrix3x3_t& m) : rows_(3), cols_(3), data_(9) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                data_[i * 3 + j] = m(i, j);
            }
        }
    }

    // Conversion from OpalMatrix6x6_t
    matrix_t(const OpalMatrix6x6_t& m) : rows_(6), cols_(6), data_(36) {
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                data_[i * 6 + j] = m(i, j);
            }
        }
    }
    
private:
    size_t rows_;
    size_t cols_;
    std::vector<double> data_;
};

// Define the conversion constructor and assignment operator for OpalMatrix_t from matrix_t
// These are defined after matrix_t is fully defined
template <int Rows, int Cols>
OpalMatrix_t<Rows, Cols>::OpalMatrix_t(const matrix_t& m) {
    // Only allow conversion if sizes match
    if (m.size1() == Rows && m.size2() == Cols) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                this->m[i][j] = m(i, j);
            }
        }
    } else {
        // Initialize to zero if sizes don't match
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                this->m[i][j] = 0.0;
            }
        }
    }
}

template <int Rows, int Cols>
OpalMatrix_t<Rows, Cols>& OpalMatrix_t<Rows, Cols>::operator=(const matrix_t& m) {
    // Only allow assignment if sizes match
    if (m.size1() == Rows && m.size2() == Cols) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                this->m[i][j] = m(i, j);
            }
        }
    }
    return *this;
}

// Simple vector class to replace boost::numeric::ublas::vector
class vector_t {
public:
    vector_t() : size_(0) {}
    
    explicit vector_t(size_t size, double init_value = 0.0)
        : size_(size), data_(size, init_value) {}
    
    double& operator()(size_t i) {
        assert(i < size_);
        return data_[i];
    }
    
    const double& operator()(size_t i) const {
        assert(i < size_);
        return data_[i];
    }
    
    size_t size() const { return size_; }
    
    void resize(size_t size, double init_value = 0.0) {
        size_ = size;
        data_.resize(size, init_value);
    }
    
private:
    size_t size_;
    std::vector<double> data_;
};

// Kokkos-friendly matrix-vector product for fixed-size matrices
template <int Rows, int Cols, class T>
KOKKOS_INLINE_FUNCTION
T prod_vector(const OpalMatrix_t<Rows, Cols>& rotation, const T& vect) {
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
T prod_vector_transpose(const OpalMatrix_t<Rows, Cols>& rotation, const T& vect) {
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
OpalMatrix_t<Cols, Rows> get_transpose(const OpalMatrix_t<Rows, Cols>& rotation) {
    OpalMatrix_t<Cols, Rows> transpose;
    for (int i = 0; i < Rows; ++i) {
        for (int j = 0; j < Cols; ++j) {
            transpose(j, i) = rotation(i, j);
        }
    }
    return transpose;
}

// Matrix-vector multiplication (legacy)
inline vector_t prod(const matrix_t& m, const vector_t& v) {
    assert(m.size2() == v.size());
    vector_t result(m.size1());
    for (size_t i = 0; i < m.size1(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < m.size2(); ++j) {
            sum += m(i, j) * v(j);
        }
        result(i) = sum;
    }
    return result;
}

// Overload for matrix-vector when second argument is matrix_t (treat as column vector)
inline matrix_t prod_matrix_vector(const matrix_t& m, const matrix_t& v) {
    assert(m.size2() == v.size1() && v.size2() == 1);
    matrix_t result(m.size1(), 1);
    for (size_t i = 0; i < m.size1(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < m.size2(); ++j) {
            sum += m(i, j) * v(j, 0);
        }
        result(i, 0) = sum;
    }
    return result;
}

// Matrix transpose (legacy)
inline matrix_t trans(const matrix_t& m) {
    matrix_t result(m.size2(), m.size1());
    for (size_t i = 0; i < m.size1(); ++i) {
        for (size_t j = 0; j < m.size2(); ++j) {
            result(j, i) = m(i, j);
        }
    }
    return result;
}

// Matrix-matrix multiplication (legacy)
inline matrix_t prod(const matrix_t& a, const matrix_t& b) {
    assert(a.size2() == b.size1());
    matrix_t result(a.size1(), b.size2());
    for (size_t i = 0; i < a.size1(); ++i) {
        for (size_t j = 0; j < b.size2(); ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < a.size2(); ++k) {
                sum += a(i, k) * b(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

// Kokkos-friendly version for 3D vectors (commonly used)
template <class T>
KOKKOS_INLINE_FUNCTION
T prod_boost_vector(const OpalMatrix3x3_t& rotation, const T& vect) {
    return prod_vector(rotation, vect);
}

// Legacy version for backward compatibility
template <class T>
T prod_boost_vector(const matrix_t& rotation, const T& vect) {
    vector_t boostVector(3);
    
    boostVector(0) = vect[0];
    boostVector(1) = vect[1];
    boostVector(2) = vect[2];
    
    boostVector = prod(rotation, boostVector);
    
    T prodVector(0.0);
    
    prodVector[0] = boostVector(0);
    prodVector[1] = boostVector(1);
    prodVector[2] = boostVector(2);
    
    return prodVector;
}

#endif
