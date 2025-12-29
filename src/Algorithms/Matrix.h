//
// This header file provides matrix and vector operations to replace boost::numeric::ublas.
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

#include <vector>
#include <cassert>

// Simple matrix class to replace boost::numeric::ublas::matrix
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
    
private:
    size_t rows_;
    size_t cols_;
    std::vector<double> data_;
};

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

// Matrix-vector multiplication
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

// Matrix transpose
inline matrix_t trans(const matrix_t& m) {
    matrix_t result(m.size2(), m.size1());
    for (size_t i = 0; i < m.size1(); ++i) {
        for (size_t j = 0; j < m.size2(); ++j) {
            result(j, i) = m(i, j);
        }
    }
    return result;
}

// Matrix-matrix multiplication
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
