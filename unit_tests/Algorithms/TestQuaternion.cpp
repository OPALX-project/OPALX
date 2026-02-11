#include "gtest/gtest.h"

#include "Algorithms/Quaternion.hpp"
#include "Algorithms/Matrix.h"
#include "opal_test_utilities/SilenceTest.h"
#include "Ippl.h"

#include <cmath>

class QuaternionTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc = 0;
        char** argv = nullptr;

        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() {
        ippl::finalize();
    }

    void SetUp() override {
        // nothing special
    }

    void TearDown() override {
        // nothing special
    }
};

TEST_F(QuaternionTest, DefaultConstructor) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q;
    EXPECT_DOUBLE_EQ(q.real(), 1.0);
    EXPECT_DOUBLE_EQ(q(0), 1.0);
    EXPECT_DOUBLE_EQ(q(1), 0.0);
    EXPECT_DOUBLE_EQ(q(2), 0.0);
    EXPECT_DOUBLE_EQ(q(3), 0.0);
}

TEST_F(QuaternionTest, ComponentConstructor) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q(1.0, 2.0, 3.0, 4.0);
    EXPECT_DOUBLE_EQ(q.real(), 1.0);
    EXPECT_DOUBLE_EQ(q(0), 1.0);
    EXPECT_DOUBLE_EQ(q(1), 2.0);
    EXPECT_DOUBLE_EQ(q(2), 3.0);
    EXPECT_DOUBLE_EQ(q(3), 4.0);
}

TEST_F(QuaternionTest, VectorConstructor) {
    OpalTestUtilities::SilenceTest silencer;

    ippl::Vector<double, 3> vec(1.0, 2.0, 3.0);
    Quaternion q(vec);
    EXPECT_DOUBLE_EQ(q.real(), 0.0);
    EXPECT_DOUBLE_EQ(q(1), 1.0);
    EXPECT_DOUBLE_EQ(q(2), 2.0);
    EXPECT_DOUBLE_EQ(q(3), 3.0);
    EXPECT_TRUE(q.isPure());
}

TEST_F(QuaternionTest, RealPlusVectorConstructor) {
    OpalTestUtilities::SilenceTest silencer;

    ippl::Vector<double, 3> vec(1.0, 2.0, 3.0);
    Quaternion q(0.5, vec);
    EXPECT_DOUBLE_EQ(q.real(), 0.5);
    EXPECT_DOUBLE_EQ(q(1), 1.0);
    EXPECT_DOUBLE_EQ(q(2), 2.0);
    EXPECT_DOUBLE_EQ(q(3), 3.0);
}

TEST_F(QuaternionTest, CopyConstructor) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q1(1.0, 2.0, 3.0, 4.0);
    Quaternion q2(q1);
    EXPECT_DOUBLE_EQ(q2.real(), 1.0);
    EXPECT_DOUBLE_EQ(q2(1), 2.0);
    EXPECT_DOUBLE_EQ(q2(2), 3.0);
    EXPECT_DOUBLE_EQ(q2(3), 4.0);
}

TEST_F(QuaternionTest, RealAndImagParts) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q(1.0, 2.0, 3.0, 4.0);
    EXPECT_DOUBLE_EQ(q.real(), 1.0);
    
    ippl::Vector<double, 3> imag = q.imag();
    EXPECT_DOUBLE_EQ(imag(0), 2.0);
    EXPECT_DOUBLE_EQ(imag(1), 3.0);
    EXPECT_DOUBLE_EQ(imag(2), 4.0);
}

TEST_F(QuaternionTest, NormAndLength) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q(1.0, 2.0, 3.0, 4.0);
    // Norm = w^2 + x^2 + y^2 + z^2 = 1 + 4 + 9 + 16 = 30
    EXPECT_DOUBLE_EQ(q.Norm(), 30.0);
    EXPECT_DOUBLE_EQ(q.length(), std::sqrt(30.0));
}

TEST_F(QuaternionTest, IsUnit) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q1(1.0, 0.0, 0.0, 0.0);
    EXPECT_TRUE(q1.isUnit());

    Quaternion q2(0.5, 0.5, 0.5, 0.5);
    EXPECT_TRUE(q2.isUnit());

    Quaternion q3(1.0, 1.0, 0.0, 0.0);
    EXPECT_FALSE(q3.isUnit());
}

TEST_F(QuaternionTest, IsPure) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q1(0.0, 1.0, 2.0, 3.0);
    EXPECT_TRUE(q1.isPure());

    Quaternion q2(1.0, 1.0, 2.0, 3.0);
    EXPECT_FALSE(q2.isPure());
}

TEST_F(QuaternionTest, IsPureUnit) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q1(0.0, 1.0, 0.0, 0.0);
    EXPECT_TRUE(q1.isPureUnit());

    Quaternion q2(0.0, 0.5, 0.5, 0.5);
    EXPECT_FALSE(q2.isPureUnit());

    Quaternion q3(0.5, 0.5, 0.5, 0.5);
    EXPECT_FALSE(q3.isPureUnit());
}

TEST_F(QuaternionTest, Conjugate) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q(1.0, 2.0, 3.0, 4.0);
    Quaternion qConj = q.conjugate();
    
    EXPECT_DOUBLE_EQ(qConj.real(), 1.0);
    EXPECT_DOUBLE_EQ(qConj(1), -2.0);
    EXPECT_DOUBLE_EQ(qConj(2), -3.0);
    EXPECT_DOUBLE_EQ(qConj(3), -4.0);
}

TEST_F(QuaternionTest, Normalize) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q(1.0, 2.0, 3.0, 4.0);
    q.normalize();
    
    EXPECT_TRUE(q.isUnit());
    EXPECT_NEAR(q.Norm(), 1.0, 1e-12);
}

TEST_F(QuaternionTest, Inverse) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q(1.0, 2.0, 3.0, 4.0);
    Quaternion qInv = q.inverse();
    
    // q * q^-1 should be identity (1, 0, 0, 0)
    Quaternion identity = q * qInv;
    
    EXPECT_NEAR(identity.real(), 1.0, 1e-12);
    EXPECT_NEAR(identity(1), 0.0, 1e-12);
    EXPECT_NEAR(identity(2), 0.0, 1e-12);
    EXPECT_NEAR(identity(3), 0.0, 1e-12);
}

TEST_F(QuaternionTest, ScalarMultiplication) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q(1.0, 2.0, 3.0, 4.0);
    Quaternion q2 = q * 2.0;
    
    EXPECT_DOUBLE_EQ(q2.real(), 2.0);
    EXPECT_DOUBLE_EQ(q2(1), 4.0);
    EXPECT_DOUBLE_EQ(q2(2), 6.0);
    EXPECT_DOUBLE_EQ(q2(3), 8.0);
}

TEST_F(QuaternionTest, ScalarDivision) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q(2.0, 4.0, 6.0, 8.0);
    Quaternion q2 = q / 2.0;
    
    EXPECT_DOUBLE_EQ(q2.real(), 1.0);
    EXPECT_DOUBLE_EQ(q2(1), 2.0);
    EXPECT_DOUBLE_EQ(q2(2), 3.0);
    EXPECT_DOUBLE_EQ(q2(3), 4.0);
}

TEST_F(QuaternionTest, QuaternionMultiplication) {
    OpalTestUtilities::SilenceTest silencer;

    // Test with simple unit quaternions
    Quaternion q1(1.0, 0.0, 0.0, 0.0);
    Quaternion q2(0.0, 1.0, 0.0, 0.0);
    
    Quaternion result = q1 * q2;
    EXPECT_DOUBLE_EQ(result.real(), 0.0);
    EXPECT_DOUBLE_EQ(result(1), 1.0);
    EXPECT_DOUBLE_EQ(result(2), 0.0);
    EXPECT_DOUBLE_EQ(result(3), 0.0);
}

TEST_F(QuaternionTest, QuaternionMultiplicationAssociativity) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q1(1.0, 2.0, 3.0, 4.0);
    Quaternion q2(5.0, 6.0, 7.0, 8.0);
    Quaternion q3(9.0, 10.0, 11.0, 12.0);
    
    Quaternion result1 = (q1 * q2) * q3;
    Quaternion result2 = q1 * (q2 * q3);
    
    EXPECT_NEAR(result1.real(), result2.real(), 1e-10);
    EXPECT_NEAR(result1(1), result2(1), 1e-10);
    EXPECT_NEAR(result1(2), result2(2), 1e-10);
    EXPECT_NEAR(result1(3), result2(3), 1e-10);
}

TEST_F(QuaternionTest, RotateVector90DegreesAroundZ) {
    OpalTestUtilities::SilenceTest silencer;

    // 90 degree rotation around z-axis
    double angle = M_PI / 2.0;
    ippl::Vector<double, 3> axis(0.0, 0.0, 1.0);
    Quaternion q(std::cos(angle / 2.0), std::sin(angle / 2.0) * axis);
    
    ippl::Vector<double, 3> vec(1.0, 0.0, 0.0);
    ippl::Vector<double, 3> rotated = q.rotate(vec);
    
    EXPECT_NEAR(rotated(0), 0.0, 1e-12);
    EXPECT_NEAR(rotated(1), 1.0, 1e-12);
    EXPECT_NEAR(rotated(2), 0.0, 1e-12);
}

TEST_F(QuaternionTest, RotateVector180DegreesAroundX) {
    OpalTestUtilities::SilenceTest silencer;

    // 180 degree rotation around x-axis
    double angle = M_PI;
    ippl::Vector<double, 3> axis(1.0, 0.0, 0.0);
    Quaternion q(std::cos(angle / 2.0), std::sin(angle / 2.0) * axis);
    
    ippl::Vector<double, 3> vec(0.0, 1.0, 0.0);
    ippl::Vector<double, 3> rotated = q.rotate(vec);
    
    EXPECT_NEAR(rotated(0), 0.0, 1e-12);
    EXPECT_NEAR(rotated(1), -1.0, 1e-12);
    EXPECT_NEAR(rotated(2), 0.0, 1e-12);
}

TEST_F(QuaternionTest, RotateVectorIdentity) {
    OpalTestUtilities::SilenceTest silencer;

    // Identity quaternion (no rotation)
    Quaternion q(1.0, 0.0, 0.0, 0.0);
    
    ippl::Vector<double, 3> vec(1.0, 2.0, 3.0);
    ippl::Vector<double, 3> rotated = q.rotate(vec);
    
    EXPECT_NEAR(rotated(0), vec(0), 1e-12);
    EXPECT_NEAR(rotated(1), vec(1), 1e-12);
    EXPECT_NEAR(rotated(2), vec(2), 1e-12);
}

TEST_F(QuaternionTest, GetRotationMatrix90DegreesZ) {
    OpalTestUtilities::SilenceTest silencer;

    // 90 degree rotation around z-axis
    double angle = M_PI / 2.0;
    ippl::Vector<double, 3> axis(0.0, 0.0, 1.0);
    Quaternion q(std::cos(angle / 2.0), std::sin(angle / 2.0) * axis);
    
    matrix3x3_t R = q.getRotationMatrix();
    
    // Expected rotation matrix for 90 degrees around z
    // [0 -1  0]
    // [1  0  0]
    // [0  0  1]
    EXPECT_NEAR(R(0, 0), 0.0, 1e-12);
    EXPECT_NEAR(R(0, 1), -1.0, 1e-12);
    EXPECT_NEAR(R(0, 2), 0.0, 1e-12);
    EXPECT_NEAR(R(1, 0), 1.0, 1e-12);
    EXPECT_NEAR(R(1, 1), 0.0, 1e-12);
    EXPECT_NEAR(R(1, 2), 0.0, 1e-12);
    EXPECT_NEAR(R(2, 0), 0.0, 1e-12);
    EXPECT_NEAR(R(2, 1), 0.0, 1e-12);
    EXPECT_NEAR(R(2, 2), 1.0, 1e-12);
}

TEST_F(QuaternionTest, GetRotationMatrixIdentity) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q(1.0, 0.0, 0.0, 0.0);
    matrix3x3_t R = q.getRotationMatrix();
    
    // Should be identity matrix
    EXPECT_NEAR(R(0, 0), 1.0, 1e-12);
    EXPECT_NEAR(R(0, 1), 0.0, 1e-12);
    EXPECT_NEAR(R(0, 2), 0.0, 1e-12);
    EXPECT_NEAR(R(1, 0), 0.0, 1e-12);
    EXPECT_NEAR(R(1, 1), 1.0, 1e-12);
    EXPECT_NEAR(R(1, 2), 0.0, 1e-12);
    EXPECT_NEAR(R(2, 0), 0.0, 1e-12);
    EXPECT_NEAR(R(2, 1), 0.0, 1e-12);
    EXPECT_NEAR(R(2, 2), 1.0, 1e-12);
}

TEST_F(QuaternionTest, MatrixConstructorRoundTrip) {
    OpalTestUtilities::SilenceTest silencer;

    // Create a rotation matrix for 45 degrees around x-axis
    double angle = M_PI / 4.0;
    double c = std::cos(angle);
    double s = std::sin(angle);
    
    matrix3x3_t R(0.0);
    R(0, 0) = 1.0; R(0, 1) = 0.0; R(0, 2) = 0.0;
    R(1, 0) = 0.0; R(1, 1) = c;   R(1, 2) = -s;
    R(2, 0) = 0.0; R(2, 1) = s;   R(2, 2) = c;
    
    // Convert to quaternion and back to matrix
    Quaternion q(R);
    matrix3x3_t R2 = q.getRotationMatrix();
    
    // Check all elements match
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(R(i, j), R2(i, j), 1e-10);
        }
    }
}

TEST_F(QuaternionTest, GetQuaternionSameVectors) {
    OpalTestUtilities::SilenceTest silencer;

    ippl::Vector<double, 3> u(1.0, 0.0, 0.0);
    ippl::Vector<double, 3> ref(1.0, 0.0, 0.0);
    
    Quaternion q = getQuaternion(u, ref);
    
    // Should be identity quaternion
    EXPECT_NEAR(q.real(), 1.0, 1e-12);
    EXPECT_TRUE(q.isUnit());
}

TEST_F(QuaternionTest, GetQuaternionOrthogonalVectors) {
    OpalTestUtilities::SilenceTest silencer;

    ippl::Vector<double, 3> u(1.0, 0.0, 0.0);
    ippl::Vector<double, 3> ref(0.0, 1.0, 0.0);
    
    Quaternion q = getQuaternion(u, ref);
    
    // Should be 90 degree rotation
    EXPECT_TRUE(q.isUnit());
    
    // Rotate u by q should give ref
    ippl::Vector<double, 3> rotated = q.rotate(u);
    EXPECT_NEAR(rotated(0), ref(0), 1e-12);
    EXPECT_NEAR(rotated(1), ref(1), 1e-12);
    EXPECT_NEAR(rotated(2), ref(2), 1e-12);
}

TEST_F(QuaternionTest, GetQuaternionOppositeVectors) {
    OpalTestUtilities::SilenceTest silencer;

    ippl::Vector<double, 3> u(1.0, 0.0, 0.0);
    ippl::Vector<double, 3> ref(-1.0, 0.0, 0.0);
    
    Quaternion q = getQuaternion(u, ref);
    
    // Should be 180 degree rotation
    EXPECT_TRUE(q.isUnit());
    EXPECT_NEAR(q.real(), 0.0, 1e-12);
    
    // Rotate u by q should give ref
    ippl::Vector<double, 3> rotated = q.rotate(u);
    EXPECT_NEAR(rotated(0), ref(0), 1e-10);
    EXPECT_NEAR(rotated(1), ref(1), 1e-10);
    EXPECT_NEAR(rotated(2), ref(2), 1e-10);
}

TEST_F(QuaternionTest, ConjugatePreservesNorm) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q(1.0, 2.0, 3.0, 4.0);
    Quaternion qConj = q.conjugate();
    
    EXPECT_DOUBLE_EQ(q.Norm(), qConj.Norm());
}

TEST_F(QuaternionTest, MultiplicationByConjugateGivesNorm) {
    OpalTestUtilities::SilenceTest silencer;

    Quaternion q(1.0, 2.0, 3.0, 4.0);
    Quaternion result = q * q.conjugate();
    
    EXPECT_NEAR(result.real(), q.Norm(), 1e-12);
    EXPECT_NEAR(result(1), 0.0, 1e-12);
    EXPECT_NEAR(result(2), 0.0, 1e-12);
    EXPECT_NEAR(result(3), 0.0, 1e-12);
}
