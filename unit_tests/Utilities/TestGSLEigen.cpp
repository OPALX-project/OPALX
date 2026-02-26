/**
 * \file TestGSLEigen.cpp
 * \brief Unit tests for GSL eigenvalue computation
 *
 * This file contains unit tests for GSL-compatible eigenvalue and eigenvector
 * computation routines. The implementation uses a simplified QR algorithm for
 * computing eigenvalues of non-symmetric matrices.
 *
 * \test GSLEigenTest::EigenNonsymm_Identity
 * Tests eigenvalue computation for the identity matrix. Verifies that all
 * eigenvalues are 1.0 (with zero imaginary part) for the identity matrix,
 * confirming that the eigenvalue algorithm correctly handles this special case.
 *
 * \test GSLEigenTest::EigenNonsymm_Diagonal
 * Tests eigenvalue computation for a diagonal matrix. Verifies that the
 * computed eigenvalues match the diagonal elements (2, 3, 4), confirming
 * that the algorithm correctly identifies eigenvalues for diagonalizable
 * matrices.
 *
 * \test GSLEigenTest::EigenNonsymmv_Identity
 * Tests computation of both eigenvalues and eigenvectors for the identity
 * matrix. Verifies that eigenvalues are correct and that eigenvectors are
 * computed (though the implementation may use simplified eigenvector computation).
 *
 * \test GSLEigenTest::EigenNonsymmv_Diagonal
 * Tests eigenvalue and eigenvector computation for a diagonal matrix. Verifies
 * that eigenvalues match diagonal elements and that the algorithm completes
 * successfully.
 *
 * \test GSLEigenTest::EigenNonsymm_ErrorHandling
 * Tests error handling for invalid input (non-square matrix). Verifies that
 * the function returns an error code when the input matrix dimensions are
 * incompatible with eigenvalue computation.
 *
 * \test GSLEigenTest::WorkspaceAllocation
 * Tests allocation and deallocation of eigenvalue computation workspaces.
 * Verifies that workspaces are correctly initialized with the specified size
 * and that memory is properly managed.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "Utilities/GSLComplex.h"
#include "Utilities/GSLEigen.h"
#include "Utilities/GSLMatrix.h"

class GSLEigenTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Test setup
    }
};

TEST_F(GSLEigenTest, EigenNonsymm_Identity) {
    gsl_matrix* A                  = gsl_matrix_alloc(2, 2);
    gsl_vector_complex* eval       = gsl_vector_complex_alloc(2);
    gsl_eigen_nonsymm_workspace* w = gsl_eigen_nonsymm_alloc(2);

    // Identity matrix: eigenvalues should be 1, 1
    gsl_matrix_set_identity(A);

    int result = gsl_eigen_nonsymm(A, eval, w);
    EXPECT_EQ(result, 0);

    // Check eigenvalues (should be 1.0)
    gsl_complex e1 = gsl_vector_complex_get(eval, 0);
    gsl_complex e2 = gsl_vector_complex_get(eval, 1);

    EXPECT_NEAR(GSL_REAL(e1), 1.0, 1e-6);
    EXPECT_NEAR(GSL_IMAG(e1), 0.0, 1e-6);
    EXPECT_NEAR(GSL_REAL(e2), 1.0, 1e-6);
    EXPECT_NEAR(GSL_IMAG(e2), 0.0, 1e-6);

    gsl_matrix_free(A);
    gsl_vector_complex_free(eval);
    gsl_eigen_nonsymm_free(w);
}

TEST_F(GSLEigenTest, EigenNonsymm_Diagonal) {
    gsl_matrix* A                  = gsl_matrix_alloc(3, 3);
    gsl_vector_complex* eval       = gsl_vector_complex_alloc(3);
    gsl_eigen_nonsymm_workspace* w = gsl_eigen_nonsymm_alloc(3);

    // Diagonal matrix: [2, 0, 0; 0, 3, 0; 0, 0, 4]
    gsl_matrix_set_zero(A);
    gsl_matrix_set(A, 0, 0, 2.0);
    gsl_matrix_set(A, 1, 1, 3.0);
    gsl_matrix_set(A, 2, 2, 4.0);

    int result = gsl_eigen_nonsymm(A, eval, w);
    EXPECT_EQ(result, 0);

    // Eigenvalues should be 2, 3, 4 (in some order)
    std::vector<double> eigenvalues(3);
    for (size_t i = 0; i < 3; ++i) {
        gsl_complex e  = gsl_vector_complex_get(eval, i);
        eigenvalues[i] = GSL_REAL(e);
    }

    std::sort(eigenvalues.begin(), eigenvalues.end());
    EXPECT_NEAR(eigenvalues[0], 2.0, 1e-6);
    EXPECT_NEAR(eigenvalues[1], 3.0, 1e-6);
    EXPECT_NEAR(eigenvalues[2], 4.0, 1e-6);

    gsl_matrix_free(A);
    gsl_vector_complex_free(eval);
    gsl_eigen_nonsymm_free(w);
}

TEST_F(GSLEigenTest, EigenNonsymmv_Identity) {
    gsl_matrix* A                   = gsl_matrix_alloc(2, 2);
    gsl_vector_complex* eval        = gsl_vector_complex_alloc(2);
    gsl_matrix_complex* evec        = gsl_matrix_complex_alloc(2, 2);
    gsl_eigen_nonsymmv_workspace* w = gsl_eigen_nonsymmv_alloc(2);

    gsl_matrix_set_identity(A);

    int result = gsl_eigen_nonsymmv(A, eval, evec, w);
    EXPECT_EQ(result, 0);

    // Check eigenvalues
    gsl_complex e1 = gsl_vector_complex_get(eval, 0);
    gsl_complex e2 = gsl_vector_complex_get(eval, 1);

    EXPECT_NEAR(GSL_REAL(e1), 1.0, 1e-6);
    EXPECT_NEAR(GSL_IMAG(e1), 0.0, 1e-6);
    EXPECT_NEAR(GSL_REAL(e2), 1.0, 1e-6);
    EXPECT_NEAR(GSL_IMAG(e2), 0.0, 1e-6);

    gsl_matrix_free(A);
    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
    gsl_eigen_nonsymmv_free(w);
}

TEST_F(GSLEigenTest, EigenNonsymmv_Diagonal) {
    gsl_matrix* A                   = gsl_matrix_alloc(2, 2);
    gsl_vector_complex* eval        = gsl_vector_complex_alloc(2);
    gsl_matrix_complex* evec        = gsl_matrix_complex_alloc(2, 2);
    gsl_eigen_nonsymmv_workspace* w = gsl_eigen_nonsymmv_alloc(2);

    // Diagonal matrix: [2, 0; 0, 3]
    gsl_matrix_set_zero(A);
    gsl_matrix_set(A, 0, 0, 2.0);
    gsl_matrix_set(A, 1, 1, 3.0);

    int result = gsl_eigen_nonsymmv(A, eval, evec, w);
    EXPECT_EQ(result, 0);

    // Check eigenvalues
    std::vector<double> eigenvalues(2);
    for (size_t i = 0; i < 2; ++i) {
        gsl_complex e  = gsl_vector_complex_get(eval, i);
        eigenvalues[i] = GSL_REAL(e);
    }

    std::sort(eigenvalues.begin(), eigenvalues.end());
    EXPECT_NEAR(eigenvalues[0], 2.0, 1e-6);
    EXPECT_NEAR(eigenvalues[1], 3.0, 1e-6);

    gsl_matrix_free(A);
    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
    gsl_eigen_nonsymmv_free(w);
}

TEST_F(GSLEigenTest, EigenNonsymm_ErrorHandling) {
    gsl_matrix* A                  = gsl_matrix_alloc(2, 3);  // Not square
    gsl_vector_complex* eval       = gsl_vector_complex_alloc(2);
    gsl_eigen_nonsymm_workspace* w = gsl_eigen_nonsymm_alloc(2);

    int result = gsl_eigen_nonsymm(A, eval, w);
    EXPECT_NE(result, 0);  // Should return error

    gsl_matrix_free(A);
    gsl_vector_complex_free(eval);
    gsl_eigen_nonsymm_free(w);
}

TEST_F(GSLEigenTest, WorkspaceAllocation) {
    gsl_eigen_nonsymm_workspace* w1 = gsl_eigen_nonsymm_alloc(5);
    EXPECT_NE(w1, nullptr);
    EXPECT_EQ(w1->n, 5);
    gsl_eigen_nonsymm_free(w1);

    gsl_eigen_nonsymmv_workspace* w2 = gsl_eigen_nonsymmv_alloc(5);
    EXPECT_NE(w2, nullptr);
    EXPECT_EQ(w2->n, 5);
    gsl_eigen_nonsymmv_free(w2);
}
