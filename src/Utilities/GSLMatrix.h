//
// GSL Matrix/Vector compatibility to replace gsl_matrix, gsl_vector
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
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef OPAL_GSL_MATRIX_HH
#define OPAL_GSL_MATRIX_HH

#include <vector>
#include <cstring>
#include <stdexcept>
#include "Utilities/GSLComplex.h"

// GSL matrix structure
struct gsl_matrix {
    size_t size1;
    size_t size2;
    size_t tda;  // physical row dimension
    double* data;
    size_t owner;  // ownership flag
    
    gsl_matrix() : size1(0), size2(0), tda(0), data(nullptr), owner(0) {}
    
    double* ptr(size_t i, size_t j) {
        return data + i * tda + j;
    }
    
    const double* ptr(size_t i, size_t j) const {
        return data + i * tda + j;
    }
};

// GSL complex matrix structure
struct gsl_matrix_complex {
    size_t size1;
    size_t size2;
    size_t tda;
    gsl_complex* data;
    size_t owner;
    
    gsl_matrix_complex() : size1(0), size2(0), tda(0), data(nullptr), owner(0) {}
    
    gsl_complex* ptr(size_t i, size_t j) {
        return data + i * tda + j;
    }
    
    const gsl_complex* ptr(size_t i, size_t j) const {
        return data + i * tda + j;
    }
};

// GSL vector structure
struct gsl_vector {
    size_t size;
    size_t stride;
    double* data;
    size_t owner;
    
    gsl_vector() : size(0), stride(0), data(nullptr), owner(0) {}
    
    double* ptr(size_t i) {
        return data + i * stride;
    }
    
    const double* ptr(size_t i) const {
        return data + i * stride;
    }
};

// GSL complex vector structure
struct gsl_vector_complex {
    size_t size;
    size_t stride;
    gsl_complex* data;
    size_t owner;
    
    gsl_vector_complex() : size(0), stride(0), data(nullptr), owner(0) {}
    
    gsl_complex* ptr(size_t i) {
        return data + i * stride;
    }
    
    const gsl_complex* ptr(size_t i) const {
        return data + i * stride;
    }
};

// Matrix allocation
inline gsl_matrix* gsl_matrix_alloc(size_t n1, size_t n2) {
    gsl_matrix* m = new gsl_matrix();
    m->size1 = n1;
    m->size2 = n2;
    m->tda = n2;
    m->data = new double[n1 * n2];
    m->owner = 1;
    std::fill(m->data, m->data + n1 * n2, 0.0);
    return m;
}

inline void gsl_matrix_free(gsl_matrix* m) {
    if (m && m->owner) {
        delete[] m->data;
    }
    delete m;
}

inline gsl_matrix_complex* gsl_matrix_complex_alloc(size_t n1, size_t n2) {
    gsl_matrix_complex* m = new gsl_matrix_complex();
    m->size1 = n1;
    m->size2 = n2;
    m->tda = n2;
    m->data = new gsl_complex[n1 * n2];
    m->owner = 1;
    for (size_t i = 0; i < n1 * n2; ++i) {
        m->data[i] = gsl_complex(0.0, 0.0);
    }
    return m;
}

inline void gsl_matrix_complex_free(gsl_matrix_complex* m) {
    if (m && m->owner) {
        delete[] m->data;
    }
    delete m;
}

// Vector allocation
inline gsl_vector* gsl_vector_alloc(size_t n) {
    gsl_vector* v = new gsl_vector();
    v->size = n;
    v->stride = 1;
    v->data = new double[n];
    v->owner = 1;
    std::fill(v->data, v->data + n, 0.0);
    return v;
}

inline void gsl_vector_free(gsl_vector* v) {
    if (v && v->owner) {
        delete[] v->data;
    }
    delete v;
}

inline gsl_vector_complex* gsl_vector_complex_alloc(size_t n) {
    gsl_vector_complex* v = new gsl_vector_complex();
    v->size = n;
    v->stride = 1;
    v->data = new gsl_complex[n];
    v->owner = 1;
    for (size_t i = 0; i < n; ++i) {
        v->data[i] = gsl_complex(0.0, 0.0);
    }
    return v;
}

inline void gsl_vector_complex_free(gsl_vector_complex* v) {
    if (v && v->owner) {
        delete[] v->data;
    }
    delete v;
}

// Matrix access
inline double* gsl_matrix_ptr(gsl_matrix* m, size_t i, size_t j) {
    return m->ptr(i, j);
}

inline const double* gsl_matrix_ptr(const gsl_matrix* m, size_t i, size_t j) {
    return const_cast<gsl_matrix*>(m)->ptr(i, j);
}

inline gsl_complex* gsl_matrix_complex_ptr(gsl_matrix_complex* m, size_t i, size_t j) {
    return m->ptr(i, j);
}

inline const gsl_complex* gsl_matrix_complex_ptr(const gsl_matrix_complex* m, size_t i, size_t j) {
    return const_cast<gsl_matrix_complex*>(m)->ptr(i, j);
}

// Vector access
inline double* gsl_vector_ptr(gsl_vector* v, size_t i) {
    return v->ptr(i);
}

inline const double* gsl_vector_ptr(const gsl_vector* v, size_t i) {
    return const_cast<gsl_vector*>(v)->ptr(i);
}

inline gsl_complex* gsl_vector_complex_ptr(gsl_vector_complex* v, size_t i) {
    return v->ptr(i);
}

inline const gsl_complex* gsl_vector_complex_ptr(const gsl_vector_complex* v, size_t i) {
    return const_cast<gsl_vector_complex*>(v)->ptr(i);
}

// Matrix operations
inline void gsl_matrix_memcpy(gsl_matrix* dest, const gsl_matrix* src) {
    if (dest->size1 != src->size1 || dest->size2 != src->size2) {
        throw std::runtime_error("gsl_matrix_memcpy: size mismatch");
    }
    std::memcpy(dest->data, src->data, src->size1 * src->size2 * sizeof(double));
}

inline void gsl_matrix_complex_memcpy(gsl_matrix_complex* dest, const gsl_matrix_complex* src) {
    if (dest->size1 != src->size1 || dest->size2 != src->size2) {
        throw std::runtime_error("gsl_matrix_complex_memcpy: size mismatch");
    }
    std::memcpy(dest->data, src->data, src->size1 * src->size2 * sizeof(gsl_complex));
}

inline void gsl_matrix_transpose_memcpy(gsl_matrix* dest, const gsl_matrix* src) {
    if (dest->size1 != src->size2 || dest->size2 != src->size1) {
        throw std::runtime_error("gsl_matrix_transpose_memcpy: size mismatch");
    }
    for (size_t i = 0; i < src->size1; ++i) {
        for (size_t j = 0; j < src->size2; ++j) {
            *gsl_matrix_ptr(dest, j, i) = *gsl_matrix_ptr(src, i, j);
        }
    }
}

inline void gsl_matrix_complex_transpose_memcpy(gsl_matrix_complex* dest, const gsl_matrix_complex* src) {
    if (dest->size1 != src->size2 || dest->size2 != src->size1) {
        throw std::runtime_error("gsl_matrix_complex_transpose_memcpy: size mismatch");
    }
    for (size_t i = 0; i < src->size1; ++i) {
        for (size_t j = 0; j < src->size2; ++j) {
            *gsl_matrix_complex_ptr(dest, j, i) = *gsl_matrix_complex_ptr(src, i, j);
        }
    }
}

inline void gsl_matrix_scale(gsl_matrix* m, double x) {
    for (size_t i = 0; i < m->size1 * m->size2; ++i) {
        m->data[i] *= x;
    }
}

inline void gsl_matrix_complex_scale(gsl_matrix_complex* m, gsl_complex x) {
    for (size_t i = 0; i < m->size1 * m->size2; ++i) {
        m->data[i] = gsl_complex_mul(m->data[i], x);
    }
}

inline void gsl_matrix_add(gsl_matrix* a, const gsl_matrix* b) {
    if (a->size1 != b->size1 || a->size2 != b->size2) {
        throw std::runtime_error("gsl_matrix_add: size mismatch");
    }
    for (size_t i = 0; i < a->size1 * a->size2; ++i) {
        a->data[i] += b->data[i];
    }
}

inline void gsl_matrix_complex_add(gsl_matrix_complex* a, const gsl_matrix_complex* b) {
    if (a->size1 != b->size1 || a->size2 != b->size2) {
        throw std::runtime_error("gsl_matrix_complex_add: size mismatch");
    }
    for (size_t i = 0; i < a->size1 * a->size2; ++i) {
        a->data[i] = gsl_complex_add(a->data[i], b->data[i]);
    }
}

// Vector operations
inline void gsl_vector_memcpy(gsl_vector* dest, const gsl_vector* src) {
    if (dest->size != src->size) {
        throw std::runtime_error("gsl_vector_memcpy: size mismatch");
    }
    for (size_t i = 0; i < src->size; ++i) {
        *gsl_vector_ptr(dest, i) = *gsl_vector_ptr(src, i);
    }
}

inline void gsl_vector_complex_memcpy(gsl_vector_complex* dest, const gsl_vector_complex* src) {
    if (dest->size != src->size) {
        throw std::runtime_error("gsl_vector_complex_memcpy: size mismatch");
    }
    for (size_t i = 0; i < src->size; ++i) {
        *gsl_vector_complex_ptr(dest, i) = *gsl_vector_complex_ptr(src, i);
    }
}

// Vector scaling
inline void gsl_vector_scale(gsl_vector* v, double x) {
    for (size_t i = 0; i < v->size; ++i) {
        *gsl_vector_ptr(v, i) *= x;
    }
}

inline void gsl_vector_complex_scale(gsl_vector_complex* v, gsl_complex x) {
    for (size_t i = 0; i < v->size; ++i) {
        *gsl_vector_complex_ptr(v, i) = gsl_complex_mul(*gsl_vector_complex_ptr(v, i), x);
    }
}

// Vector addition
inline void gsl_vector_add(gsl_vector* a, const gsl_vector* b) {
    if (a->size != b->size) {
        throw std::runtime_error("gsl_vector_add: size mismatch");
    }
    for (size_t i = 0; i < a->size; ++i) {
        *gsl_vector_ptr(a, i) += *gsl_vector_ptr(b, i);
    }
}

inline void gsl_vector_complex_add(gsl_vector_complex* a, const gsl_vector_complex* b) {
    if (a->size != b->size) {
        throw std::runtime_error("gsl_vector_complex_add: size mismatch");
    }
    for (size_t i = 0; i < a->size; ++i) {
        *gsl_vector_complex_ptr(a, i) = gsl_complex_add(*gsl_vector_complex_ptr(a, i), 
                                                         *gsl_vector_complex_ptr(b, i));
    }
}

// Matrix set/get functions
inline void gsl_matrix_set(gsl_matrix* m, size_t i, size_t j, double x) {
    *gsl_matrix_ptr(m, i, j) = x;
}

inline double gsl_matrix_get(const gsl_matrix* m, size_t i, size_t j) {
    return *gsl_matrix_ptr(m, i, j);
}

inline void gsl_matrix_set_all(gsl_matrix* m, double x) {
    for (size_t i = 0; i < m->size1 * m->size2; ++i) {
        m->data[i] = x;
    }
}

inline void gsl_matrix_set_zero(gsl_matrix* m) {
    gsl_matrix_set_all(m, 0.0);
}

inline void gsl_matrix_set_identity(gsl_matrix* m) {
    gsl_matrix_set_zero(m);
    size_t n = (m->size1 < m->size2) ? m->size1 : m->size2;
    for (size_t i = 0; i < n; ++i) {
        gsl_matrix_set(m, i, i, 1.0);
    }
}

// Complex matrix set/get functions
inline void gsl_matrix_complex_set(gsl_matrix_complex* m, size_t i, size_t j, gsl_complex x) {
    *gsl_matrix_complex_ptr(m, i, j) = x;
}

inline gsl_complex gsl_matrix_complex_get(const gsl_matrix_complex* m, size_t i, size_t j) {
    return *gsl_matrix_complex_ptr(m, i, j);
}

inline void gsl_matrix_complex_set_zero(gsl_matrix_complex* m) {
    for (size_t i = 0; i < m->size1 * m->size2; ++i) {
        m->data[i] = gsl_complex(0.0, 0.0);
    }
}

// Vector set/get functions
inline void gsl_vector_set(gsl_vector* v, size_t i, double x) {
    *gsl_vector_ptr(v, i) = x;
}

inline double gsl_vector_get(const gsl_vector* v, size_t i) {
    return *gsl_vector_ptr(v, i);
}

inline void gsl_vector_set_all(gsl_vector* v, double x) {
    for (size_t i = 0; i < v->size; ++i) {
        *gsl_vector_ptr(v, i) = x;
    }
}

inline void gsl_vector_set_zero(gsl_vector* v) {
    gsl_vector_set_all(v, 0.0);
}

// Complex vector set/get functions
inline void gsl_vector_complex_set(gsl_vector_complex* v, size_t i, gsl_complex x) {
    *gsl_vector_complex_ptr(v, i) = x;
}

inline gsl_complex gsl_vector_complex_get(const gsl_vector_complex* v, size_t i) {
    return *gsl_vector_complex_ptr(v, i);
}

inline void gsl_vector_complex_set_zero(gsl_vector_complex* v) {
    for (size_t i = 0; i < v->size; ++i) {
        *gsl_vector_complex_ptr(v, i) = gsl_complex(0.0, 0.0);
    }
}

#endif // OPAL_GSL_MATRIX_HH

