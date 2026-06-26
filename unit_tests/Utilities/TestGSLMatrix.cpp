/**
 * \file TestGSLMatrix.cpp
 * \brief Unit tests for GSL matrix and vector operations
 *
 * This file contains unit tests for GSL-compatible matrix and vector data
 * structures and operations. The implementation provides wrappers around
 * standard containers with GSL-compatible APIs.
 *
 * \test GSLMatrixTest::MatrixAllocation
 * Tests matrix memory allocation and deallocation. Verifies that matrices
 * are correctly allocated with specified dimensions and that memory is
 * properly managed.
 *
 * \test GSLMatrixTest::MatrixSetGet
 * Tests setting and getting individual matrix elements. Verifies that values
 * can be stored and retrieved correctly using row and column indices.
 *
 * \test GSLMatrixTest::MatrixSetAll
 * Tests setting all matrix elements to a constant value. Verifies that
 * \c gsl_matrix_set_all() correctly initializes all elements.
 *
 * \test GSLMatrixTest::MatrixSetZero
 * Tests zeroing a matrix. Verifies that \c gsl_matrix_set_zero() sets all
 * elements to zero, even if the matrix was previously initialized with
 * non-zero values.
 *
 * \test GSLMatrixTest::MatrixSetIdentity
 * Tests creating an identity matrix. Verifies that \c gsl_matrix_set_identity()
 * sets diagonal elements to 1 and off-diagonal elements to 0.
 *
 * \test GSLMatrixTest::VectorAllocation
 * Tests vector memory allocation and deallocation. Verifies correct
 * initialization of vector size and data pointer.
 *
 * \test GSLMatrixTest::VectorSetGet
 * Tests setting and getting individual vector elements. Verifies that values
 * can be stored and retrieved correctly using element indices.
 *
 * \test GSLMatrixTest::VectorSetAll
 * Tests setting all vector elements to a constant value. Verifies uniform
 * initialization of all elements.
 *
 * \test GSLMatrixTest::VectorScale
 * Tests vector scaling operation: \f$\mathbf{v} \leftarrow \alpha \mathbf{v}\f$.
 * Verifies that all elements are multiplied by the scaling factor.
 *
 * \test GSLMatrixTest::VectorAdd
 * Tests vector addition: \f$\mathbf{v}_1 \leftarrow \mathbf{v}_1 + \mathbf{v}_2\f$.
 * Verifies element-wise addition of two vectors.
 *
 * \test GSLMatrixTest::ComplexMatrix
 * Tests complex matrix operations. Verifies that complex matrices can be
 * allocated, elements can be set and retrieved, and complex number operations
 * work correctly.
 *
 * \test GSLMatrixTest::ComplexVector
 * Tests complex vector operations including scaling with complex scalars.
 * Verifies that complex vectors can be manipulated correctly.
 */

#include <gtest/gtest.h>
#include <cmath>
#include "Utilities/GSLComplex.h"
#include "Utilities/GSLMatrix.h"

class GSLMatrixTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Test setup
    }
};

TEST_F(GSLMatrixTest, MatrixAllocation) {
    gsl_matrix* m = gsl_matrix_alloc(3, 4);
    EXPECT_EQ(m->size1, 3);
    EXPECT_EQ(m->size2, 4);
    EXPECT_NE(m->data, nullptr);

    gsl_matrix_free(m);
}

TEST_F(GSLMatrixTest, MatrixSetGet) {
    gsl_matrix* m = gsl_matrix_alloc(3, 3);

    gsl_matrix_set(m, 0, 0, 1.0);
    gsl_matrix_set(m, 0, 1, 2.0);
    gsl_matrix_set(m, 1, 0, 3.0);
    gsl_matrix_set(m, 1, 1, 4.0);

    EXPECT_DOUBLE_EQ(gsl_matrix_get(m, 0, 0), 1.0);
    EXPECT_DOUBLE_EQ(gsl_matrix_get(m, 0, 1), 2.0);
    EXPECT_DOUBLE_EQ(gsl_matrix_get(m, 1, 0), 3.0);
    EXPECT_DOUBLE_EQ(gsl_matrix_get(m, 1, 1), 4.0);

    gsl_matrix_free(m);
}

TEST_F(GSLMatrixTest, MatrixSetAll) {
    gsl_matrix* m = gsl_matrix_alloc(2, 2);
    gsl_matrix_set_all(m, 5.0);

    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(gsl_matrix_get(m, i, j), 5.0);
        }
    }

    gsl_matrix_free(m);
}

TEST_F(GSLMatrixTest, MatrixSetZero) {
    gsl_matrix* m = gsl_matrix_alloc(2, 2);
    gsl_matrix_set_all(m, 1.0);
    gsl_matrix_set_zero(m);

    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(gsl_matrix_get(m, i, j), 0.0);
        }
    }

    gsl_matrix_free(m);
}

TEST_F(GSLMatrixTest, MatrixSetIdentity) {
    gsl_matrix* m = gsl_matrix_alloc(3, 3);
    gsl_matrix_set_identity(m);

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_DOUBLE_EQ(gsl_matrix_get(m, i, j), expected);
        }
    }

    gsl_matrix_free(m);
}

TEST_F(GSLMatrixTest, VectorAllocation) {
    gsl_vector* v = gsl_vector_alloc(5);
    EXPECT_EQ(v->size, 5);
    EXPECT_NE(v->data, nullptr);

    gsl_vector_free(v);
}

TEST_F(GSLMatrixTest, VectorSetGet) {
    gsl_vector* v = gsl_vector_alloc(5);

    gsl_vector_set(v, 0, 1.0);
    gsl_vector_set(v, 2, 3.0);
    gsl_vector_set(v, 4, 5.0);

    EXPECT_DOUBLE_EQ(gsl_vector_get(v, 0), 1.0);
    EXPECT_DOUBLE_EQ(gsl_vector_get(v, 2), 3.0);
    EXPECT_DOUBLE_EQ(gsl_vector_get(v, 4), 5.0);

    gsl_vector_free(v);
}

TEST_F(GSLMatrixTest, VectorSetAll) {
    gsl_vector* v = gsl_vector_alloc(4);
    gsl_vector_set_all(v, 2.5);

    for (size_t i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ(gsl_vector_get(v, i), 2.5);
    }

    gsl_vector_free(v);
}

TEST_F(GSLMatrixTest, VectorScale) {
    gsl_vector* v = gsl_vector_alloc(3);
    gsl_vector_set(v, 0, 1.0);
    gsl_vector_set(v, 1, 2.0);
    gsl_vector_set(v, 2, 3.0);

    gsl_vector_scale(v, 2.0);

    EXPECT_DOUBLE_EQ(gsl_vector_get(v, 0), 2.0);
    EXPECT_DOUBLE_EQ(gsl_vector_get(v, 1), 4.0);
    EXPECT_DOUBLE_EQ(gsl_vector_get(v, 2), 6.0);

    gsl_vector_free(v);
}

TEST_F(GSLMatrixTest, VectorAdd) {
    gsl_vector* v1 = gsl_vector_alloc(3);
    gsl_vector* v2 = gsl_vector_alloc(3);

    gsl_vector_set(v1, 0, 1.0);
    gsl_vector_set(v1, 1, 2.0);
    gsl_vector_set(v1, 2, 3.0);

    gsl_vector_set(v2, 0, 4.0);
    gsl_vector_set(v2, 1, 5.0);
    gsl_vector_set(v2, 2, 6.0);

    gsl_vector_add(v1, v2);

    EXPECT_DOUBLE_EQ(gsl_vector_get(v1, 0), 5.0);
    EXPECT_DOUBLE_EQ(gsl_vector_get(v1, 1), 7.0);
    EXPECT_DOUBLE_EQ(gsl_vector_get(v1, 2), 9.0);

    gsl_vector_free(v1);
    gsl_vector_free(v2);
}

TEST_F(GSLMatrixTest, ComplexMatrix) {
    gsl_matrix_complex* m = gsl_matrix_complex_alloc(2, 2);

    gsl_complex z1 = gsl_complex_rect(1.0, 2.0);
    gsl_complex z2 = gsl_complex_rect(3.0, 4.0);

    gsl_matrix_complex_set(m, 0, 0, z1);
    gsl_matrix_complex_set(m, 0, 1, z2);

    gsl_complex result = gsl_matrix_complex_get(m, 0, 0);
    EXPECT_DOUBLE_EQ(GSL_REAL(result), 1.0);
    EXPECT_DOUBLE_EQ(GSL_IMAG(result), 2.0);

    result = gsl_matrix_complex_get(m, 0, 1);
    EXPECT_DOUBLE_EQ(GSL_REAL(result), 3.0);
    EXPECT_DOUBLE_EQ(GSL_IMAG(result), 4.0);

    gsl_matrix_complex_free(m);
}

TEST_F(GSLMatrixTest, ComplexVector) {
    gsl_vector_complex* v = gsl_vector_complex_alloc(3);

    gsl_complex z = gsl_complex_rect(1.0, 2.0);
    gsl_vector_complex_set(v, 1, z);

    gsl_complex result = gsl_vector_complex_get(v, 1);
    EXPECT_DOUBLE_EQ(GSL_REAL(result), 1.0);
    EXPECT_DOUBLE_EQ(GSL_IMAG(result), 2.0);

    gsl_vector_complex_scale(v, gsl_complex_rect(2.0, 0.0));
    result = gsl_vector_complex_get(v, 1);
    EXPECT_DOUBLE_EQ(GSL_REAL(result), 2.0);
    EXPECT_DOUBLE_EQ(GSL_IMAG(result), 4.0);

    gsl_vector_complex_free(v);
}
