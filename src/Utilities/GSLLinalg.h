//
// GSL Linear Algebra compatibility to replace gsl_linalg
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

#ifndef OPAL_GSL_LINALG_HH
#define OPAL_GSL_LINALG_HH

#include "Utilities/GSLMatrix.h"
#include "Utilities/GSLComplex.h"
#include <vector>
#include <algorithm>
#include <cmath>

/// \see Documentation on https://www.gnu.org/software/gsl/doc/html/linalg.htm
/// \see Implementation on https://www.gnu.org/software/gsl/
/// \brief Permutation vector for LU decomposition.
struct gsl_permutation {
    size_t size;
    size_t* data;
    
    gsl_permutation() : size(0), data(nullptr) {}
};

/// \brief Allocate an identity permutation of size \f$n\f$.
/// \param n Input: permutation size.
/// \return Output: pointer to newly allocated permutation.
inline gsl_permutation* gsl_permutation_alloc(size_t n) {
    gsl_permutation* p = new gsl_permutation();
    p->size = n;
    p->data = new size_t[n];
    for (size_t i = 0; i < n; ++i) {
        p->data[i] = i;
    }
    return p;
}

/// \brief Free a permutation allocated by \c gsl_permutation_alloc.
/// \param p Input: permutation to release (can be null).
inline void gsl_permutation_free(gsl_permutation* p) {
    if (p) {
        delete[] p->data;
        delete p;
    }
}

/// \brief In-place LU decomposition with partial pivoting (real).
/// \details Computes \f$P A = L U\f$, storing \f$L\f$ (unit lower) and \f$U\f$
/// in \p A. The permutation \p p encodes \f$P\f$.
/// \param A Input/Output: matrix overwritten with LU factors.
/// \param p Output: permutation encoding row swaps.
/// \param signum Output: sign of the permutation (\f$\pm 1\f$).
/// \return Output: 0 on success, -1 on failure (e.g. non-square or singular).
inline int gsl_linalg_LU_decomp(gsl_matrix* A, gsl_permutation* p, int* signum) {
    if (A->size1 != A->size2) {
        return -1;  // Error: not square
    }
    
    size_t n = A->size1;
    *signum = 1;
    
    // Initialize permutation
    for (size_t i = 0; i < n; ++i) {
        p->data[i] = i;
    }
    
    // Perform LU decomposition with partial pivoting
    for (size_t k = 0; k < n - 1; ++k) {
        // Find pivot
        size_t max_row = k;
        double max_val = std::abs(*gsl_matrix_ptr(A, k, k));
        for (size_t i = k + 1; i < n; ++i) {
            double val = std::abs(*gsl_matrix_ptr(A, i, k));
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }
        
        // Swap rows if needed
        if (max_row != k) {
            *signum = -*signum;
            std::swap(p->data[k], p->data[max_row]);
            for (size_t j = 0; j < n; ++j) {
                std::swap(*gsl_matrix_ptr(A, k, j), *gsl_matrix_ptr(A, max_row, j));
            }
        }
        
        // Check for singularity
        if (std::abs(*gsl_matrix_ptr(A, k, k)) < 1e-15) {
            return -1;  // Singular matrix
        }
        
        // Eliminate
        for (size_t i = k + 1; i < n; ++i) {
            double factor = *gsl_matrix_ptr(A, i, k) / *gsl_matrix_ptr(A, k, k);
            *gsl_matrix_ptr(A, i, k) = factor;
            for (size_t j = k + 1; j < n; ++j) {
                *gsl_matrix_ptr(A, i, j) -= factor * *gsl_matrix_ptr(A, k, j);
            }
        }
    }
    
    return 0;  // Success
}

/// \brief Determinant from LU decomposition (real).
/// \details Computes \f$\det(A) = \mathrm{signum} \prod_i U_{ii}\f$.
/// \param LU Input: LU-decomposed matrix.
/// \param signum Input: sign from \c gsl_linalg_LU_decomp.
/// \return Output: determinant value.
inline double gsl_linalg_LU_det(const gsl_matrix* LU, int signum) {
    size_t n = LU->size1;
    double det = static_cast<double>(signum);
    for (size_t i = 0; i < n; ++i) {
        det *= *gsl_matrix_ptr(LU, i, i);
    }
    return det;
}

/// \brief Invert a matrix from its LU decomposition (real).
/// \details Computes \f$A^{-1}\f$ using the stored \f$L,U\f$ factors and permutation.
/// \param LU Input: LU-decomposed matrix.
/// \param p Input: permutation from \c gsl_linalg_LU_decomp.
/// \param inverse Output: matrix to receive \f$A^{-1}\f$.
/// \return Output: 0 on success, -1 on size mismatch.
inline int gsl_linalg_LU_invert(const gsl_matrix* LU, const gsl_permutation* p, gsl_matrix* inverse) {
    size_t n = LU->size1;
    if (inverse->size1 != n || inverse->size2 != n) {
        return -1;
    }
    
    // Initialize inverse as permuted identity: P*I
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < n; ++i) {
            *gsl_matrix_ptr(inverse, i, j) = (p->data[i] == j) ? 1.0 : 0.0;
        }
    }
    
    // Forward substitution: L*y = P*I (solve for y)
    // L is unit lower triangular stored in LU below diagonal
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (size_t k = 0; k < i; ++k) {
                sum += *gsl_matrix_ptr(LU, i, k) * *gsl_matrix_ptr(inverse, k, j);
            }
            *gsl_matrix_ptr(inverse, i, j) -= sum;
        }
    }
    
    // Back substitution: U*x = y (solve for x, which is the inverse)
    // U is upper triangular stored in LU on and above diagonal
    for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
        for (size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = i + 1; k < n; ++k) {
                sum += *gsl_matrix_ptr(LU, i, k) * *gsl_matrix_ptr(inverse, k, j);
            }
            *gsl_matrix_ptr(inverse, i, j) = (*gsl_matrix_ptr(inverse, i, j) - sum) / *gsl_matrix_ptr(LU, i, i);
        }
    }
    
    return 0;
}

/// \brief In-place LU decomposition with partial pivoting (complex).
/// \details Computes \f$P A = L U\f$, storing \f$L\f$ and \f$U\f$ in \p A.
/// \param A Input/Output: matrix overwritten with LU factors.
/// \param p Output: permutation encoding row swaps.
/// \param signum Output: sign of the permutation (\f$\pm 1\f$).
/// \return Output: 0 on success, -1 on failure (e.g. non-square or singular).
inline int gsl_linalg_complex_LU_decomp(gsl_matrix_complex* A, gsl_permutation* p, int* signum) {
    if (A->size1 != A->size2) {
        return -1;
    }
    
    size_t n = A->size1;
    *signum = 1;
    
    // Initialize permutation
    for (size_t i = 0; i < n; ++i) {
        p->data[i] = i;
    }
    
    // Perform LU decomposition with partial pivoting
    for (size_t k = 0; k < n - 1; ++k) {
        // Find pivot
        size_t max_row = k;
        double max_val = gsl_complex_abs(*gsl_matrix_complex_ptr(A, k, k));
        for (size_t i = k + 1; i < n; ++i) {
            double val = gsl_complex_abs(*gsl_matrix_complex_ptr(A, i, k));
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }
        
        // Swap rows if needed
        if (max_row != k) {
            *signum = -*signum;
            std::swap(p->data[k], p->data[max_row]);
            for (size_t j = 0; j < n; ++j) {
                std::swap(*gsl_matrix_complex_ptr(A, k, j), *gsl_matrix_complex_ptr(A, max_row, j));
            }
        }
        
        // Check for singularity
        if (gsl_complex_abs(*gsl_matrix_complex_ptr(A, k, k)) < 1e-15) {
            return -1;
        }
        
        // Eliminate
        for (size_t i = k + 1; i < n; ++i) {
            gsl_complex factor = gsl_complex_div(*gsl_matrix_complex_ptr(A, i, k), 
                                                 *gsl_matrix_complex_ptr(A, k, k));
            *gsl_matrix_complex_ptr(A, i, k) = factor;
            for (size_t j = k + 1; j < n; ++j) {
                *gsl_matrix_complex_ptr(A, i, j) = gsl_complex_sub(*gsl_matrix_complex_ptr(A, i, j),
                                                                    gsl_complex_mul(factor, *gsl_matrix_complex_ptr(A, k, j)));
            }
        }
    }
    
    return 0;
}

/// \brief Determinant from LU decomposition (complex).
/// \details Computes \f$\det(A) = \mathrm{signum} \prod_i U_{ii}\f$.
/// \param LU Input: LU-decomposed matrix.
/// \param signum Input: sign from \c gsl_linalg_complex_LU_decomp.
/// \return Output: determinant value.
inline gsl_complex gsl_linalg_complex_LU_det(const gsl_matrix_complex* LU, int signum) {
    size_t n = LU->size1;
    gsl_complex det = gsl_complex(static_cast<double>(signum), 0.0);
    for (size_t i = 0; i < n; ++i) {
        det = gsl_complex_mul(det, *gsl_matrix_complex_ptr(LU, i, i));
    }
    return det;
}

/// \brief Invert a matrix from its LU decomposition (complex).
/// \details Computes \f$A^{-1}\f$ using the stored \f$L,U\f$ factors and permutation.
/// \param LU Input: LU-decomposed matrix.
/// \param p Input: permutation from \c gsl_linalg_complex_LU_decomp.
/// \param inverse Output: matrix to receive \f$A^{-1}\f$.
/// \return Output: 0 on success, -1 on size mismatch.
inline int gsl_linalg_complex_LU_invert(const gsl_matrix_complex* LU, const gsl_permutation* p, gsl_matrix_complex* inverse) {
    size_t n = LU->size1;
    if (inverse->size1 != n || inverse->size2 != n) {
        return -1;
    }
    
    // Initialize inverse as permuted identity: P*I
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < n; ++i) {
            *gsl_matrix_complex_ptr(inverse, i, j) = (p->data[i] == j) ? gsl_complex(1.0, 0.0) : gsl_complex(0.0, 0.0);
        }
    }
    
    // Forward substitution: L*y = P*I (solve for y)
    // L is unit lower triangular stored in LU below diagonal
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < n; ++i) {
            gsl_complex sum = gsl_complex(0.0, 0.0);
            for (size_t k = 0; k < i; ++k) {
                sum = gsl_complex_add(sum, gsl_complex_mul(*gsl_matrix_complex_ptr(LU, i, k), 
                                                           *gsl_matrix_complex_ptr(inverse, k, j)));
            }
            *gsl_matrix_complex_ptr(inverse, i, j) = gsl_complex_sub(*gsl_matrix_complex_ptr(inverse, i, j), sum);
        }
    }
    
    // Back substitution: U*x = y (solve for x, which is the inverse)
    // U is upper triangular stored in LU on and above diagonal
    for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
        for (size_t j = 0; j < n; ++j) {
            gsl_complex sum = gsl_complex(0.0, 0.0);
            for (size_t k = i + 1; k < n; ++k) {
                sum = gsl_complex_add(sum, gsl_complex_mul(*gsl_matrix_complex_ptr(LU, i, k), 
                                                           *gsl_matrix_complex_ptr(inverse, k, j)));
            }
            *gsl_matrix_complex_ptr(inverse, i, j) = gsl_complex_div(gsl_complex_sub(*gsl_matrix_complex_ptr(inverse, i, j), sum),
                                                                      *gsl_matrix_complex_ptr(LU, i, i));
        }
    }
    
    return 0;
}

/// \brief Alias for \c gsl_linalg_complex_LU_decomp.
/// \param A Input/Output: matrix overwritten with LU factors.
/// \param p Output: permutation encoding row swaps.
/// \param signum Output: sign of the permutation (\f$\pm 1\f$).
/// \return Output: 0 on success, -1 on failure.
inline int gsl_linalg_LU_decomp_complex(gsl_matrix_complex* A, gsl_permutation* p, int* signum) {
    return gsl_linalg_complex_LU_decomp(A, p, signum);
}

/// \brief Alias for \c gsl_linalg_complex_LU_det.
/// \param LU Input: LU-decomposed matrix.
/// \param signum Input: sign from decomposition.
/// \return Output: determinant value.
inline gsl_complex gsl_linalg_LU_det_complex(const gsl_matrix_complex* LU, int signum) {
    return gsl_linalg_complex_LU_det(LU, signum);
}

/// \brief Alias for \c gsl_linalg_complex_LU_invert.
/// \param LU Input: LU-decomposed matrix.
/// \param p Input: permutation from decomposition.
/// \param inverse Output: matrix to receive \f$A^{-1}\f$.
/// \return Output: 0 on success, -1 on size mismatch.
inline int gsl_linalg_LU_invert_complex(const gsl_matrix_complex* LU, const gsl_permutation* p, gsl_matrix_complex* inverse) {
    return gsl_linalg_complex_LU_invert(LU, p, inverse);
}

#endif // OPAL_GSL_LINALG_HH

