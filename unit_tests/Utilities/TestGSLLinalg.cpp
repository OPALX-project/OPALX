/**
 * \file TestGSLLinalg.cpp
 * \brief Unit tests for GSL linear algebra routines
 *
 * This file contains unit tests for GSL-compatible linear algebra operations
 * including LU decomposition, determinant computation, and matrix inversion
 * for both real and complex matrices.
 *
 * \test GSLLinalgTest::LUDecomp
 * Tests LU decomposition with partial pivoting: \f$\mathbf{P} \mathbf{A} = \mathbf{L}
 * \mathbf{U}\f$, where \f$\mathbf{P}\f$ is a permutation matrix, \f$\mathbf{L}\f$ is lower
 * triangular, and \f$\mathbf{U}\f$ is upper triangular. Verifies that the
 * decomposition completes successfully for a tridiagonal matrix.
 *
 * \test GSLLinalgTest::LUDet
 * Tests determinant computation using LU decomposition:
 * \f$\det(\mathbf{A}) = \det(\mathbf{P}) \det(\mathbf{L}) \det(\mathbf{U})\f$.
 * Verifies correct computation of determinant for a 2×2 matrix with known
 * determinant value.
 *
 * \test GSLLinalgTest::LUInvert
 * Tests matrix inversion using LU decomposition: \f$\mathbf{A}^{-1}\f$.
 * Verifies that the computed inverse satisfies \f$\mathbf{A} \mathbf{A}^{-1} = \mathbf{I}\f$
 * by checking individual elements of the inverse matrix.
 *
 * \test GSLLinalgTest::LUInvertIdentity
 * Tests inversion of the identity matrix. Verifies that \f$\mathbf{I}^{-1} = \mathbf{I}\f$,
 * confirming that the inversion algorithm handles the special case correctly.
 *
 * \test GSLLinalgTest::ComplexLUDecomp
 * Tests LU decomposition for complex matrices. Verifies that the decomposition
 * algorithm works correctly with complex-valued matrices, including the identity
 * matrix case.
 *
 * \test GSLLinalgTest::ComplexLUDet
 * Tests determinant computation for complex matrices. Verifies that the
 * determinant of the identity matrix is correctly computed as 1 (real part)
 * with zero imaginary part.
 *
 * \test GSLLinalgTest::ComplexLUInvert
 * Tests matrix inversion for complex matrices. Verifies that complex matrix
 * inversion works correctly, including the identity matrix case.
 */

#include <gtest/gtest.h>
#include <cmath>
#include "Utilities/GSLComplex.h"
#include "Utilities/GSLLinalg.h"
#include "Utilities/GSLMatrix.h"

class GSLLinalgTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Test setup
    }
};

TEST_F(GSLLinalgTest, LUDecomp) {
    gsl_matrix* A      = gsl_matrix_alloc(3, 3);
    gsl_permutation* p = gsl_permutation_alloc(3);
    int signum;

    // A = [2 1 0]
    //     [1 2 1]
    //     [0 1 2]
    gsl_matrix_set(A, 0, 0, 2.0);
    gsl_matrix_set(A, 0, 1, 1.0);
    gsl_matrix_set(A, 0, 2, 0.0);
    gsl_matrix_set(A, 1, 0, 1.0);
    gsl_matrix_set(A, 1, 1, 2.0);
    gsl_matrix_set(A, 1, 2, 1.0);
    gsl_matrix_set(A, 2, 0, 0.0);
    gsl_matrix_set(A, 2, 1, 1.0);
    gsl_matrix_set(A, 2, 2, 2.0);

    int result = gsl_linalg_LU_decomp(A, p, &signum);
    EXPECT_EQ(result, 0);

    gsl_matrix_free(A);
    gsl_permutation_free(p);
}

TEST_F(GSLLinalgTest, LUDet) {
    gsl_matrix* A      = gsl_matrix_alloc(2, 2);
    gsl_permutation* p = gsl_permutation_alloc(2);
    int signum;

    // A = [1 2]
    //     [3 4]
    // det(A) = 1*4 - 2*3 = -2
    gsl_matrix_set(A, 0, 0, 1.0);
    gsl_matrix_set(A, 0, 1, 2.0);
    gsl_matrix_set(A, 1, 0, 3.0);
    gsl_matrix_set(A, 1, 1, 4.0);

    gsl_linalg_LU_decomp(A, p, &signum);
    double det = gsl_linalg_LU_det(A, signum);

    EXPECT_NEAR(det, -2.0, 1e-10);

    gsl_matrix_free(A);
    gsl_permutation_free(p);
}

TEST_F(GSLLinalgTest, LUInvert) {
    gsl_matrix* A      = gsl_matrix_alloc(2, 2);
    gsl_matrix* invA   = gsl_matrix_alloc(2, 2);
    gsl_permutation* p = gsl_permutation_alloc(2);
    int signum;

    // A = [1 2]
    //     [3 4]
    // A^-1 = [-2  1]
    //         [1.5 -0.5]
    gsl_matrix_set(A, 0, 0, 1.0);
    gsl_matrix_set(A, 0, 1, 2.0);
    gsl_matrix_set(A, 1, 0, 3.0);
    gsl_matrix_set(A, 1, 1, 4.0);

    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_invert(A, p, invA);

    EXPECT_NEAR(gsl_matrix_get(invA, 0, 0), -2.0, 1e-10);
    EXPECT_NEAR(gsl_matrix_get(invA, 0, 1), 1.0, 1e-10);
    EXPECT_NEAR(gsl_matrix_get(invA, 1, 0), 1.5, 1e-10);
    EXPECT_NEAR(gsl_matrix_get(invA, 1, 1), -0.5, 1e-10);

    gsl_matrix_free(A);
    gsl_matrix_free(invA);
    gsl_permutation_free(p);
}

TEST_F(GSLLinalgTest, LUInvertIdentity) {
    gsl_matrix* A      = gsl_matrix_alloc(3, 3);
    gsl_matrix* invA   = gsl_matrix_alloc(3, 3);
    gsl_permutation* p = gsl_permutation_alloc(3);
    int signum;

    gsl_matrix_set_identity(A);

    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_invert(A, p, invA);

    // I^-1 = I
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(gsl_matrix_get(invA, i, j), expected, 1e-10);
        }
    }

    gsl_matrix_free(A);
    gsl_matrix_free(invA);
    gsl_permutation_free(p);
}

TEST_F(GSLLinalgTest, ComplexLUDecomp) {
    gsl_matrix_complex* A = gsl_matrix_complex_alloc(2, 2);
    gsl_permutation* p    = gsl_permutation_alloc(2);
    int signum;

    gsl_complex z1 = gsl_complex_rect(1.0, 0.0);
    gsl_complex z2 = gsl_complex_rect(0.0, 0.0);

    gsl_matrix_complex_set(A, 0, 0, z1);
    gsl_matrix_complex_set(A, 0, 1, z2);
    gsl_matrix_complex_set(A, 1, 0, z2);
    gsl_matrix_complex_set(A, 1, 1, z1);

    int result = gsl_linalg_LU_decomp_complex(A, p, &signum);
    EXPECT_EQ(result, 0);

    gsl_matrix_complex_free(A);
    gsl_permutation_free(p);
}

TEST_F(GSLLinalgTest, ComplexLUDet) {
    gsl_matrix_complex* A = gsl_matrix_complex_alloc(2, 2);
    gsl_permutation* p    = gsl_permutation_alloc(2);
    int signum;

    // Identity matrix: det = 1
    gsl_complex z1 = gsl_complex_rect(1.0, 0.0);
    gsl_complex z2 = gsl_complex_rect(0.0, 0.0);

    gsl_matrix_complex_set(A, 0, 0, z1);
    gsl_matrix_complex_set(A, 0, 1, z2);
    gsl_matrix_complex_set(A, 1, 0, z2);
    gsl_matrix_complex_set(A, 1, 1, z1);

    gsl_linalg_LU_decomp_complex(A, p, &signum);
    gsl_complex det = gsl_linalg_LU_det_complex(A, signum);

    EXPECT_NEAR(GSL_REAL(det), 1.0, 1e-10);
    EXPECT_NEAR(GSL_IMAG(det), 0.0, 1e-10);

    gsl_matrix_complex_free(A);
    gsl_permutation_free(p);
}

TEST_F(GSLLinalgTest, ComplexLUInvert) {
    gsl_matrix_complex* A    = gsl_matrix_complex_alloc(2, 2);
    gsl_matrix_complex* invA = gsl_matrix_complex_alloc(2, 2);
    gsl_permutation* p       = gsl_permutation_alloc(2);
    int signum;

    // Identity matrix
    gsl_complex z1 = gsl_complex_rect(1.0, 0.0);
    gsl_complex z2 = gsl_complex_rect(0.0, 0.0);

    gsl_matrix_complex_set(A, 0, 0, z1);
    gsl_matrix_complex_set(A, 0, 1, z2);
    gsl_matrix_complex_set(A, 1, 0, z2);
    gsl_matrix_complex_set(A, 1, 1, z1);

    gsl_linalg_LU_decomp_complex(A, p, &signum);
    gsl_linalg_LU_invert_complex(A, p, invA);

    // I^-1 = I
    gsl_complex result = gsl_matrix_complex_get(invA, 0, 0);
    EXPECT_NEAR(GSL_REAL(result), 1.0, 1e-10);
    EXPECT_NEAR(GSL_IMAG(result), 0.0, 1e-10);

    gsl_matrix_complex_free(A);
    gsl_matrix_complex_free(invA);
    gsl_permutation_free(p);
}
