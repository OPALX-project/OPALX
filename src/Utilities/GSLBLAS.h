//
// GSL BLAS compatibility to replace gsl_blas
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

#ifndef OPAL_GSL_BLAS_HH
#define OPAL_GSL_BLAS_HH

#include <cstring>
#include "Utilities/GSLComplex.h"
#include "Utilities/GSLMatrix.h"

/// \see Documentation on https://www.gnu.org/software/gsl/doc/html/blas.html
/// \see Implementation on https://www.gnu.org/software/gsl/
/// \brief Transpose operation selector for BLAS routines.
/// \details
/// - \c CblasNoTrans: \f$\mathrm{op}(A) = A\f$
/// - \c CblasTrans: \f$\mathrm{op}(A) = A^T\f$
/// - \c CblasConjTrans: \f$\mathrm{op}(A) = A^H\f$ (conjugate transpose)
enum CBLAS_TRANSPOSE { CblasNoTrans = 111, CblasTrans = 112, CblasConjTrans = 113 };

/// \brief Real matrix-matrix multiply and accumulate.
/// \details Computes \f$C \leftarrow \alpha \,\mathrm{op}(A)\,\mathrm{op}(B) + \beta\, C\f$
/// where \f$\mathrm{op}(X)\f$ is chosen by \p TransA and \p TransB.
/// \param TransA Input: operation applied to \p A.
/// \param TransB Input: operation applied to \p B.
/// \param alpha Input: scalar multiplier for \f$\mathrm{op}(A)\mathrm{op}(B)\f$.
/// \param A Input: left matrix operand.
/// \param B Input: right matrix operand.
/// \param beta Input: scalar multiplier for existing \p C.
/// \param C Input/Output: accumulation target updated in-place.
inline void gsl_blas_dgemm(
    CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, double alpha, const gsl_matrix* A,
    const gsl_matrix* B, double beta, gsl_matrix* C
) {
    size_t m = (TransA == CblasNoTrans) ? A->size1 : A->size2;
    size_t n = (TransB == CblasNoTrans) ? B->size2 : B->size1;
    size_t k = (TransA == CblasNoTrans) ? A->size2 : A->size1;

    if (C->size1 != m || C->size2 != n) {
        throw std::runtime_error("gsl_blas_dgemm: size mismatch");
    }

    // Initialize C with beta*C
    if (beta != 1.0) {
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j) {
                *gsl_matrix_ptr(C, i, j) *= beta;
            }
        }
    }

    // Perform multiplication
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (size_t l = 0; l < k; ++l) {
                double a_val =
                    (TransA == CblasNoTrans) ? *gsl_matrix_ptr(A, i, l) : *gsl_matrix_ptr(A, l, i);
                double b_val =
                    (TransB == CblasNoTrans) ? *gsl_matrix_ptr(B, l, j) : *gsl_matrix_ptr(B, j, l);
                sum += a_val * b_val;
            }
            *gsl_matrix_ptr(C, i, j) += alpha * sum;
        }
    }
}

/// \brief Complex matrix-matrix multiply and accumulate.
/// \details Computes \f$C \leftarrow \alpha \,\mathrm{op}(A)\,\mathrm{op}(B) + \beta\, C\f$
/// with \f$\mathrm{op}(X)\f$ defined by \p TransA/\p TransB (including conjugate transpose).
/// \param TransA Input: operation applied to \p A.
/// \param TransB Input: operation applied to \p B.
/// \param alpha Input: scalar multiplier for \f$\mathrm{op}(A)\mathrm{op}(B)\f$.
/// \param A Input: left matrix operand.
/// \param B Input: right matrix operand.
/// \param beta Input: scalar multiplier for existing \p C.
/// \param C Input/Output: accumulation target updated in-place.
inline void gsl_blas_zgemm(
    CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, gsl_complex alpha, const gsl_matrix_complex* A,
    const gsl_matrix_complex* B, gsl_complex beta, gsl_matrix_complex* C
) {
    size_t m = (TransA == CblasNoTrans) ? A->size1 : A->size2;
    size_t n = (TransB == CblasNoTrans) ? B->size2 : B->size1;
    size_t k = (TransA == CblasNoTrans) ? A->size2 : A->size1;

    if (C->size1 != m || C->size2 != n) {
        throw std::runtime_error("gsl_blas_zgemm: size mismatch");
    }

    // Initialize C with beta*C
    if (beta.dat[0] != 1.0 || beta.dat[1] != 0.0) {
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j) {
                *gsl_matrix_complex_ptr(C, i, j) =
                    gsl_complex_mul(*gsl_matrix_complex_ptr(C, i, j), beta);
            }
        }
    }

    // Perform multiplication
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            gsl_complex sum = gsl_complex(0.0, 0.0);
            for (size_t l = 0; l < k; ++l) {
                gsl_complex a_val = (TransA == CblasNoTrans) ? *gsl_matrix_complex_ptr(A, i, l)
                                                             : *gsl_matrix_complex_ptr(A, l, i);
                if (TransA == CblasConjTrans) {
                    a_val = gsl_complex_conjugate(a_val);
                }
                gsl_complex b_val = (TransB == CblasNoTrans) ? *gsl_matrix_complex_ptr(B, l, j)
                                                             : *gsl_matrix_complex_ptr(B, j, l);
                if (TransB == CblasConjTrans) {
                    b_val = gsl_complex_conjugate(b_val);
                }
                sum = gsl_complex_add(sum, gsl_complex_mul(a_val, b_val));
            }
            *gsl_matrix_complex_ptr(C, i, j) =
                gsl_complex_add(*gsl_matrix_complex_ptr(C, i, j), gsl_complex_mul(alpha, sum));
        }
    }
}

/// \brief Real matrix-vector multiply and accumulate.
/// \details Computes \f$y \leftarrow \alpha\,\mathrm{op}(A)\,x + \beta\,y\f$
/// with \f$\mathrm{op}(A) = A\f$ or \f$A^T\f$.
/// \param TransA Input: operation applied to \p A.
/// \param alpha Input: scalar multiplier for \f$\mathrm{op}(A)x\f$.
/// \param A Input: matrix operand.
/// \param x Input: vector operand.
/// \param beta Input: scalar multiplier for existing \p y.
/// \param y Input/Output: result vector updated in-place.
inline void gsl_blas_dgemv(
    CBLAS_TRANSPOSE TransA, double alpha, const gsl_matrix* A, const gsl_vector* x, double beta,
    gsl_vector* y
) {
    size_t m = A->size1;
    size_t n = A->size2;

    if (TransA == CblasNoTrans) {
        if (x->size != n || y->size != m) {
            throw std::runtime_error("gsl_blas_dgemv: size mismatch");
        }
        // y = beta*y
        for (size_t i = 0; i < m; ++i) {
            *gsl_vector_ptr(y, i) *= beta;
        }
        // y += alpha*A*x
        for (size_t i = 0; i < m; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < n; ++j) {
                sum += *gsl_matrix_ptr(A, i, j) * *gsl_vector_ptr(x, j);
            }
            *gsl_vector_ptr(y, i) += alpha * sum;
        }
    } else {
        if (x->size != m || y->size != n) {
            throw std::runtime_error("gsl_blas_dgemv: size mismatch");
        }
        // y = beta*y
        for (size_t i = 0; i < n; ++i) {
            *gsl_vector_ptr(y, i) *= beta;
        }
        // y += alpha*A^T*x
        for (size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < m; ++j) {
                sum += *gsl_matrix_ptr(A, j, i) * *gsl_vector_ptr(x, j);
            }
            *gsl_vector_ptr(y, i) += alpha * sum;
        }
    }
}

/// \brief Complex matrix-vector multiply and accumulate.
/// \details Computes \f$y \leftarrow \alpha\,\mathrm{op}(A)\,x + \beta\,y\f$
/// with \f$\mathrm{op}(A) = A, A^T,\f$ or \f$A^H\f$.
/// \param TransA Input: operation applied to \p A.
/// \param alpha Input: scalar multiplier for \f$\mathrm{op}(A)x\f$.
/// \param A Input: matrix operand.
/// \param x Input: vector operand.
/// \param beta Input: scalar multiplier for existing \p y.
/// \param y Input/Output: result vector updated in-place.
inline void gsl_blas_zgemv(
    CBLAS_TRANSPOSE TransA, gsl_complex alpha, const gsl_matrix_complex* A,
    const gsl_vector_complex* x, gsl_complex beta, gsl_vector_complex* y
) {
    size_t m = A->size1;
    size_t n = A->size2;

    if (TransA == CblasNoTrans) {
        if (x->size != n || y->size != m) {
            throw std::runtime_error("gsl_blas_zgemv: size mismatch");
        }
        // y = beta*y
        for (size_t i = 0; i < m; ++i) {
            *gsl_vector_complex_ptr(y, i) = gsl_complex_mul(*gsl_vector_complex_ptr(y, i), beta);
        }
        // y += alpha*A*x
        for (size_t i = 0; i < m; ++i) {
            gsl_complex sum = gsl_complex(0.0, 0.0);
            for (size_t j = 0; j < n; ++j) {
                sum = gsl_complex_add(
                    sum,
                    gsl_complex_mul(*gsl_matrix_complex_ptr(A, i, j), *gsl_vector_complex_ptr(x, j))
                );
            }
            *gsl_vector_complex_ptr(y, i) =
                gsl_complex_add(*gsl_vector_complex_ptr(y, i), gsl_complex_mul(alpha, sum));
        }
    } else {
        if (x->size != m || y->size != n) {
            throw std::runtime_error("gsl_blas_zgemv: size mismatch");
        }
        // y = beta*y
        for (size_t i = 0; i < n; ++i) {
            *gsl_vector_complex_ptr(y, i) = gsl_complex_mul(*gsl_vector_complex_ptr(y, i), beta);
        }
        // y += alpha*A^T*x or A^H*x
        for (size_t i = 0; i < n; ++i) {
            gsl_complex sum = gsl_complex(0.0, 0.0);
            for (size_t j = 0; j < m; ++j) {
                gsl_complex a_val = *gsl_matrix_complex_ptr(A, j, i);
                if (TransA == CblasConjTrans) {
                    a_val = gsl_complex_conjugate(a_val);
                }
                sum = gsl_complex_add(sum, gsl_complex_mul(a_val, *gsl_vector_complex_ptr(x, j)));
            }
            *gsl_vector_complex_ptr(y, i) =
                gsl_complex_add(*gsl_vector_complex_ptr(y, i), gsl_complex_mul(alpha, sum));
        }
    }
}

#endif  // OPAL_GSL_BLAS_HH
