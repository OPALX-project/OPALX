/**
 * \file TestGSLBLAS.cpp
 * \brief Unit tests for GSL BLAS (Basic Linear Algebra Subprograms) operations
 *
 * This file contains unit tests for GSL-compatible BLAS operations including
 * matrix-matrix and matrix-vector multiplications for both real and complex
 * data types.
 *
 * \test GSLBLASTest::DGEMM_NoTrans
 * Tests double-precision general matrix-matrix multiplication without
 * transposition: \f$\mathbf{C} = \alpha \mathbf{A} \mathbf{B} + \beta \mathbf{C}\f$
 * with \f$\alpha = 1\f$ and \f$\beta = 0\f$. Verifies correct computation
 * of matrix product for non-square matrices.
 *
 * \test GSLBLASTest::DGEMM_WithBeta
 * Tests matrix-matrix multiplication with non-zero \f$\beta\f$ coefficient.
 * Verifies that the operation \f$\mathbf{C} = \alpha \mathbf{A} \mathbf{B} + \beta \mathbf{C}\f$
 * correctly combines the matrix product with the existing matrix \f$\mathbf{C}\f$.
 *
 * \test GSLBLASTest::DGEMM_Transpose
 * Tests matrix-matrix multiplication with transposed matrices. Verifies that
 * \f$\mathbf{C} = \mathbf{A}^T \mathbf{B}\f$ is computed correctly, including
 * the case where \f$\mathbf{A}^T\f$ and \f$\mathbf{B}\f$ have compatible
 * dimensions.
 *
 * \test GSLBLASTest::DGEMV
 * Tests double-precision general matrix-vector multiplication:
 * \f$\mathbf{y} = \alpha \mathbf{A} \mathbf{x} + \beta \mathbf{y}\f$.
 * Verifies correct computation of matrix-vector product and combination with
 * existing vector \f$\mathbf{y}\f$.
 *
 * \test GSLBLASTest::ZGEMM_Complex
 * Tests complex matrix-matrix multiplication using double-precision complex
 * numbers. Verifies that complex matrix products are computed correctly,
 * including the identity matrix case \f$\mathbf{I} \mathbf{I} = \mathbf{I}\f$.
 */

#include <gtest/gtest.h>
#include <cmath>
#include "Utilities/GSLBLAS.h"
#include "Utilities/GSLComplex.h"
#include "Utilities/GSLMatrix.h"

class GSLBLASTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Test setup
    }
};

TEST_F(GSLBLASTest, DGEMM_NoTrans) {
    // Test C = A * B
    gsl_matrix* A = gsl_matrix_alloc(2, 3);
    gsl_matrix* B = gsl_matrix_alloc(3, 2);
    gsl_matrix* C = gsl_matrix_alloc(2, 2);

    // A = [1 2 3]
    //     [4 5 6]
    gsl_matrix_set(A, 0, 0, 1.0);
    gsl_matrix_set(A, 0, 1, 2.0);
    gsl_matrix_set(A, 0, 2, 3.0);
    gsl_matrix_set(A, 1, 0, 4.0);
    gsl_matrix_set(A, 1, 1, 5.0);
    gsl_matrix_set(A, 1, 2, 6.0);

    // B = [1 2]
    //     [3 4]
    //     [5 6]
    gsl_matrix_set(B, 0, 0, 1.0);
    gsl_matrix_set(B, 0, 1, 2.0);
    gsl_matrix_set(B, 1, 0, 3.0);
    gsl_matrix_set(B, 1, 1, 4.0);
    gsl_matrix_set(B, 2, 0, 5.0);
    gsl_matrix_set(B, 2, 1, 6.0);

    gsl_matrix_set_zero(C);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, C);

    // Expected: C = [22 28]
    //                [49 64]
    EXPECT_NEAR(gsl_matrix_get(C, 0, 0), 22.0, 1e-10);
    EXPECT_NEAR(gsl_matrix_get(C, 0, 1), 28.0, 1e-10);
    EXPECT_NEAR(gsl_matrix_get(C, 1, 0), 49.0, 1e-10);
    EXPECT_NEAR(gsl_matrix_get(C, 1, 1), 64.0, 1e-10);

    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_matrix_free(C);
}

TEST_F(GSLBLASTest, DGEMM_WithBeta) {
    // Test C = alpha*A*B + beta*C
    gsl_matrix* A = gsl_matrix_alloc(2, 2);
    gsl_matrix* B = gsl_matrix_alloc(2, 2);
    gsl_matrix* C = gsl_matrix_alloc(2, 2);

    gsl_matrix_set_identity(A);
    gsl_matrix_set_identity(B);
    gsl_matrix_set_all(C, 1.0);

    // C = 2.0 * I * I + 3.0 * C = 2*I + 3*C = 2*I + 3*ones = [5 3]
    //                                                          [3 5]
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 2.0, A, B, 3.0, C);

    EXPECT_NEAR(gsl_matrix_get(C, 0, 0), 5.0, 1e-10);
    EXPECT_NEAR(gsl_matrix_get(C, 0, 1), 3.0, 1e-10);
    EXPECT_NEAR(gsl_matrix_get(C, 1, 0), 3.0, 1e-10);
    EXPECT_NEAR(gsl_matrix_get(C, 1, 1), 5.0, 1e-10);

    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_matrix_free(C);
}

TEST_F(GSLBLASTest, DGEMM_Transpose) {
    // Test C = A^T * B
    gsl_matrix* A = gsl_matrix_alloc(3, 2);
    gsl_matrix* B = gsl_matrix_alloc(3, 2);
    gsl_matrix* C = gsl_matrix_alloc(2, 2);

    gsl_matrix_set_identity(A);
    gsl_matrix_set_identity(B);

    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, B, 0.0, C);

    // A^T * B = I^T * I = I * I = I
    EXPECT_NEAR(gsl_matrix_get(C, 0, 0), 1.0, 1e-10);
    EXPECT_NEAR(gsl_matrix_get(C, 0, 1), 0.0, 1e-10);
    EXPECT_NEAR(gsl_matrix_get(C, 1, 0), 0.0, 1e-10);
    EXPECT_NEAR(gsl_matrix_get(C, 1, 1), 1.0, 1e-10);

    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_matrix_free(C);
}

TEST_F(GSLBLASTest, DGEMV) {
    // Test y = alpha*A*x + beta*y
    gsl_matrix* A = gsl_matrix_alloc(2, 3);
    gsl_vector* x = gsl_vector_alloc(3);
    gsl_vector* y = gsl_vector_alloc(2);

    gsl_matrix_set_identity(A);
    gsl_vector_set_all(x, 1.0);
    gsl_vector_set_all(y, 2.0);

    // y = 1.0 * A * x + 1.0 * y = [1, 1, 0]^T + [2, 2]^T = [3, 2]^T
    // (Note: A is 2x3 identity-like, so first two elements of x are used)
    gsl_matrix_set(A, 0, 0, 1.0);
    gsl_matrix_set(A, 0, 1, 0.0);
    gsl_matrix_set(A, 0, 2, 0.0);
    gsl_matrix_set(A, 1, 0, 0.0);
    gsl_matrix_set(A, 1, 1, 1.0);
    gsl_matrix_set(A, 1, 2, 0.0);

    gsl_blas_dgemv(CblasNoTrans, 1.0, A, x, 1.0, y);

    EXPECT_NEAR(gsl_vector_get(y, 0), 3.0, 1e-10);
    EXPECT_NEAR(gsl_vector_get(y, 1), 3.0, 1e-10);

    gsl_matrix_free(A);
    gsl_vector_free(x);
    gsl_vector_free(y);
}

TEST_F(GSLBLASTest, ZGEMM_Complex) {
    gsl_matrix_complex* A = gsl_matrix_complex_alloc(2, 2);
    gsl_matrix_complex* B = gsl_matrix_complex_alloc(2, 2);
    gsl_matrix_complex* C = gsl_matrix_complex_alloc(2, 2);

    gsl_complex z1 = gsl_complex_rect(1.0, 0.0);
    gsl_complex z2 = gsl_complex_rect(0.0, 0.0);

    gsl_matrix_complex_set(A, 0, 0, z1);
    gsl_matrix_complex_set(A, 0, 1, z2);
    gsl_matrix_complex_set(A, 1, 0, z2);
    gsl_matrix_complex_set(A, 1, 1, z1);

    gsl_matrix_complex_set(B, 0, 0, z1);
    gsl_matrix_complex_set(B, 0, 1, z2);
    gsl_matrix_complex_set(B, 1, 0, z2);
    gsl_matrix_complex_set(B, 1, 1, z1);

    gsl_matrix_complex_set_zero(C);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, z1, A, B, z2, C);

    // I * I = I
    gsl_complex result = gsl_matrix_complex_get(C, 0, 0);
    EXPECT_NEAR(GSL_REAL(result), 1.0, 1e-10);
    EXPECT_NEAR(GSL_IMAG(result), 0.0, 1e-10);

    gsl_matrix_complex_free(A);
    gsl_matrix_complex_free(B);
    gsl_matrix_complex_free(C);
}
