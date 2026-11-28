#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include "Algorithms/TpsLinearMapAlgebra.h"
#include "Ippl.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

namespace {
    using Matrix = matrix6x6_t;
    Matrix identity() {
        Matrix result(0.);
        for (int i = 0; i < 6; ++i)
            result(i, i) = 1;
        return result;
    }
    Matrix multiply(const Matrix& a, const Matrix& b) {
        Matrix result(0.);
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                for (int k = 0; k < 6; ++k)
                    result(i, j) += a(i, k) * b(k, j);
        return result;
    }
    void expectMatrix(const Matrix& a, const Matrix& b, double tolerance = 3e-13) {
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                EXPECT_NEAR(a(i, j), b(i, j), tolerance * std::max(1., std::abs(b(i, j))))
                        << i << ',' << j;
    }
    class TpsLinearMapAlgebraTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }
        static void TearDownTestSuite() { ippl::finalize(); }
    };
}  // namespace

TEST_F(TpsLinearMapAlgebraTest, NonsymmetricRightSolveWithPivotsAndMixedCoordinateScales) {
    std::mt19937 generator(692);
    std::uniform_real_distribution<double> random(-.3, .3);
    for (int trial = 0; trial < 12; ++trial) {
        Matrix input(0.), expected(0.);
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) {
                input(i, j)    = (i == j ? 2. : random(generator)) * std::pow(10., -4 + j);
                expected(i, j) = random(generator);
            }
        // Force row interchanges, preserving nonsymmetry and mixed column scales.
        for (int j = 0; j < 6; ++j)
            std::swap(input(0, j), input(4, j));
        const auto output     = multiply(expected, input);
        const auto savedInput = input, savedOutput = output;
        double condition  = 0;
        const auto result = opalx::maps::solveInputDifferences(input, output, condition);
        expectMatrix(result, expected);
        expectMatrix(multiply(result, input), output);
        expectMatrix(input, savedInput, 0);
        expectMatrix(output, savedOutput, 0);
        EXPECT_GT(condition, 1e4);
        EXPECT_LT(condition, 1e7);
    }
}

TEST_F(TpsLinearMapAlgebraTest, ConditionUsesOriginalInfinityNormNotTransposeNorm) {
    auto input  = identity();
    input(0, 1) = 2;
    input(0, 2) = 3;
    // A=I+N, N^2=0: inverse=I-N. Row norms give6*6=36; column norms give4*4=16.
    double condition = 0;
    expectMatrix(opalx::maps::solveInputDifferences(input, input, condition), identity());
    EXPECT_DOUBLE_EQ(condition, 36);
}

TEST_F(TpsLinearMapAlgebraTest, PivotFloorAndConditionLimitRemainIndependent) {
    double condition = 0;
    for (double scale : {std::nextafter(1e-18, 1.), 1e-18, 1e-10}) {
        Matrix input(0.);
        for (int i = 0; i < 6; ++i)
            input(i, i) = scale;
        expectMatrix(opalx::maps::solveInputDifferences(input, input, condition), identity());
        EXPECT_NEAR(condition, 1., 3e-16);
    }
    Matrix tiny(0.);
    for (int i = 0; i < 6; ++i)
        tiny(i, i) = std::nextafter(1e-18, 0.);
    EXPECT_THROW(opalx::maps::solveInputDifferences(tiny, tiny, condition), OpalException);
    auto mixed  = identity();
    mixed(0, 0) = 1e-12;
    expectMatrix(opalx::maps::solveInputDifferences(mixed, mixed, condition), identity());
    EXPECT_DOUBLE_EQ(condition, 1e12);
    mixed(0, 0) = 0.5e-12;
    EXPECT_THROW(opalx::maps::solveInputDifferences(mixed, mixed, condition), OpalException);
    EXPECT_DOUBLE_EQ(condition, 2e12);
    auto singular = identity();
    for (int j = 0; j < 6; ++j)
        singular(3, j) = singular(2, j);
    EXPECT_THROW(opalx::maps::solveInputDifferences(singular, singular, condition), OpalException);
}

TEST_F(TpsLinearMapAlgebraTest, PivotFloorAppliesToOriginalMatrixNotItsTranspose) {
    Matrix input(0.);
    for (int i = 0; i < 6; ++i)
        input(i, i) = 2e-18;
    input(0, 1) = 100e-18;
    input(1, 0) = 1e-18;
    input(1, 1) = 1e-18;
    // Original first two U pivots: 2e-18,-49e-18. Transposed:100e-18,.98e-18.
    // Applying the floor to a new factorization of A^T would incorrectly reject it.
    double condition = 0;
    expectMatrix(opalx::maps::solveInputDifferences(input, input, condition), identity());
    EXPECT_LT(condition, 1e3);
}

TEST_F(TpsLinearMapAlgebraTest, DeterminantOfCoupledBlockTriangularMap) {
    Matrix map(0.);
    const double blocks[3][4] = {{.1, 2, 3, 4}, {2, 1, 4, 3}, {1, -1, 2, 5}};
    double expected           = 1;
    for (int pair = 0; pair < 3; ++pair) {
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                map(2 * pair + i, 2 * pair + j) = blocks[pair][2 * i + j];
        expected *= blocks[pair][0] * blocks[pair][3] - blocks[pair][1] * blocks[pair][2];
    }
    map(0, 3)    = .73;
    map(1, 4)    = -1.1;
    map(2, 5)    = 2.4;
    const auto d = opalx::maps::diagnoseLinearMap(map);
    EXPECT_NEAR(d.determinantResidual, std::abs(expected - 1), 3e-13);
}

TEST_F(TpsLinearMapAlgebraTest, SingularDeterminantIsDiagnosticNotFatal) {
    auto map               = identity();
    map(0, 0)              = 0;
    const auto diagnostics = opalx::maps::diagnoseLinearMap(map);
    EXPECT_DOUBLE_EQ(diagnostics.determinantResidual, 1);
    EXPECT_DOUBLE_EQ(diagnostics.symplecticResidual, 1);
    for (int shift : {1, 2, 3, 5}) {
        Matrix permutation(0.);
        for (int i = 0; i < 6; ++i)
            permutation(i, (i + shift) % 6) = 1;
        const auto d = opalx::maps::diagnoseLinearMap(permutation);
        EXPECT_DOUBLE_EQ(d.determinantResidual, shift == 2 ? 0 : 2);
    }
}

TEST_F(TpsLinearMapAlgebraTest, CanonicalFormUsesMaximumEntryAndIndependentDenseOracle) {
    Matrix j(0.);
    for (int pair = 0; pair < 3; ++pair) {
        j(2 * pair, 2 * pair + 1) = 1;
        j(2 * pair + 1, 2 * pair) = -1;
    }
    auto map  = identity();
    map(0, 0) = 2;
    map(1, 1) = .5;
    map(2, 3) = .4;
    map(4, 5) = -.3;
    auto d    = opalx::maps::diagnoseLinearMap(map);
    EXPECT_NEAR(d.determinantResidual, 0, 1e-15);
    EXPECT_NEAR(d.symplecticResidual, 0, 1e-15);
    map(0, 3) = .7;
    map(5, 1) = -.2;
    Matrix transpose(0.);
    for (int i = 0; i < 6; ++i)
        for (int k = 0; k < 6; ++k)
            transpose(i, k) = map(k, i);
    const auto oracle = multiply(transpose, multiply(j, map));
    double maximum = 0, sumSquares = 0;
    for (int i = 0; i < 6; ++i)
        for (int k = 0; k < 6; ++k) {
            const auto residual = oracle(i, k) - j(i, k);
            maximum             = std::max(maximum, std::abs(residual));
            sumSquares += residual * residual;
        }
    d = opalx::maps::diagnoseLinearMap(map);
    EXPECT_NEAR(d.symplecticResidual, maximum, 1e-15);
    EXPECT_GT(std::sqrt(sumSquares), maximum);
}

TEST_F(TpsLinearMapAlgebraTest, RejectsNonfiniteInputsAndOverflow) {
    double condition = 0;
    for (double value :
         {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        auto invalid  = identity();
        invalid(1, 2) = value;
        EXPECT_THROW(
                opalx::maps::solveInputDifferences(invalid, identity(), condition), OpalException);
        EXPECT_THROW(
                opalx::maps::solveInputDifferences(identity(), invalid, condition), OpalException);
        EXPECT_THROW(opalx::maps::diagnoseLinearMap(invalid), OpalException);
    }
    auto huge = identity();
    for (int i = 0; i < 6; ++i)
        huge(i, i) = 1e200;
    EXPECT_THROW(opalx::maps::diagnoseLinearMap(huge), OpalException);
    for (int i = 0; i < 6; ++i)
        huge(i, i) = 1e60;  // Finite canonical residual, determinant overflows.
    EXPECT_THROW(opalx::maps::diagnoseLinearMap(huge), OpalException);
    Matrix tiny(0.);
    for (int i = 0; i < 6; ++i) {
        tiny(i, i) = 1e-18;
        huge(i, i) = 1e308;
    }
    EXPECT_THROW(opalx::maps::solveInputDifferences(tiny, huge, condition), OpalException);
}
