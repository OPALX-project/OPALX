#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "Algorithms/MatchedCovariance.h"
#include "Ippl.h"
#include "gtest/gtest.h"

namespace {
    using Matrix = MatchedCovariance::Matrix;
    Matrix identity() {
        Matrix a{};
        for (unsigned i = 0; i < 4; ++i)
            a[i][i] = 1;
        return a;
    }
    // Independent fixed-size reference products, separate from TPS/Kokkos kernels.
    Matrix product(const Matrix& a, const Matrix& b) {
        Matrix r{};
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                for (unsigned k = 0; k < 4; ++k)
                    r[i][j] += a[i][k] * b[k][j];
        return r;
    }
    Matrix transpose(const Matrix& a) {
        Matrix r{};
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                r[i][j] = a[j][i];
        return r;
    }
    Matrix rotations(double q1 = .13, double q2 = .31) {
        Matrix a{};
        for (unsigned k = 0; k < 2; ++k) {
            const unsigned i   = 2 * k;
            const double phase = 2 * std::acos(-1.) * (k ? q2 : q1);
            a[i][i] = a[i + 1][i + 1] = std::cos(phase);
            a[i][i + 1]               = std::sin(phase);
            a[i + 1][i]               = -std::sin(phase);
        }
        return a;
    }
    ClosedOrbitOptics optics(const Matrix& a) {
        ClosedOrbitOptics o;
        o.matrix          = a;
        o.launch.momentum = ippl::Vector<double, 3>(0, 0, 2);
        o.launch.chargeE  = 1;
        o.launch.massEV   = 938272088.16;
        return o;
    }
    MatchedCovariance::Result matched(
            const ClosedOrbitOptics& optics, const std::array<double, 2>& eps) {
        try {
            return MatchedCovariance::compute(optics, eps);
        } catch (const OpalException& ex) {
            throw std::runtime_error(ex.where() + ": " + ex.what());
        }
    }
    void expectMatrix(const Matrix& actual, const Matrix& expected, double tolerance = 2e-12) {
        double scale = 0;
        for (const auto& row : expected)
            for (double x : row)
                scale = std::max(scale, std::abs(x));
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                EXPECT_NEAR(actual[i][j], expected[i][j], tolerance * scale) << i << ',' << j;
    }
    class MatchedCovarianceTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }
        static void TearDownTestSuite() { ippl::finalize(); }
    };
}  // namespace

TEST_F(MatchedCovarianceTest, AnalyticTwissNormalizedRmsWithoutPiAndIndependentPeriodicity) {
    Matrix a{}, expected{};
    const std::array<double, 2> eps{2e-6, 7e-6}, beta{3.2, .8}, alpha{-.7, 1.1};
    for (unsigned k = 0; k < 2; ++k) {
        const unsigned i   = 2 * k;
        const double phase = 2 * std::acos(-1.) * (k ? .31 : .13), c = std::cos(phase),
                     s     = std::sin(phase);
        const double gamma = (1 + alpha[k] * alpha[k]) / beta[k];
        a[i][i]            = c + alpha[k] * s;
        a[i][i + 1]        = beta[k] * s;
        a[i + 1][i]        = -gamma * s;
        a[i + 1][i + 1]    = c - alpha[k] * s;
        expected[i][i]     = eps[k] * beta[k];
        expected[i][i + 1] = expected[i + 1][i] = -eps[k] * alpha[k];
        expected[i + 1][i + 1]                  = eps[k] * gamma;
    }
    const auto actual = matched(optics(a), eps);
    expectMatrix(actual.covariance, expected);
    expectMatrix(product(product(a, actual.covariance), transpose(a)), expected);
    for (unsigned k = 0; k < 2; ++k) {
        const unsigned i  = 2 * k;
        const auto& sigma = actual.covariance;
        EXPECT_NEAR(
                std::sqrt(sigma[i][i] * sigma[i + 1][i + 1] - sigma[i][i + 1] * sigma[i + 1][i]),
                eps[k], eps[k] * 1e-12);
    }
    EXPECT_LT(actual.symplecticResidual, 1e-13);
    EXPECT_LT(actual.periodicResidual, 1e-13);
}

TEST_F(MatchedCovarianceTest, CoupledModesUndoScalesAndOrderByPositivePhase) {
    const double c = std::cos(.6), s = std::sin(.6);
    const Matrix coupling{{{c, 0, s, 0}, {0, c, 0, s}, {-s, 0, c, 0}, {0, -s, 0, c}}};
    // First physical mode has phase .83; positive-imaginary representative is .17.
    // Second is .11, so prescribed emittances are intentionally reversed by plane.
    auto source   = optics(product(product(coupling, rotations(.83, .11)), transpose(coupling)));
    source.scales = {.01, .4, .02, .2};
    Matrix diagonal{};
    diagonal[0][0] = diagonal[1][1] = 7e-6;
    diagonal[2][2] = diagonal[3][3] = 2e-6;
    const auto result               = matched(source, {2e-6, 7e-6});
    expectMatrix(result.covariance, product(product(coupling, diagonal), transpose(coupling)));
    EXPECT_NEAR(result.tunes[0], .11, 1e-14);
    EXPECT_NEAR(result.tunes[1], .17, 1e-14);
    EXPECT_GT(std::abs(result.covariance[0][2]), 1e-6);
    auto unit   = source;
    unit.scales = {1, 1, 1, 1};
    expectMatrix(matched(unit, {2e-6, 7e-6}).covariance, result.covariance);
}

TEST_F(MatchedCovarianceTest, NoncanonicalMagneticSectionChargeSignUnitsAndRotation) {
    for (double charge : {-1., 2.}) {
        auto source           = optics(rotations());
        source.launch.massEV  = 2.99792458e8;
        source.launch.chargeE = charge;
        source.launch.section = CoordinateSystemTrafo(
                ippl::Vector<double, 3>(2, -1, 3), Quaternion(std::cos(.3), std::sin(.3), 0, 0));
        const auto labField =
                source.launch.section.rotateFrom(ippl::Vector<double, 3>(.1, -.2, 3.));
        source.normalMagneticField = source.launch.section.rotateTo(labField)[2];
        const auto c               = MatchedCovariance::canonicalTransform(
                source.normalMagneticField, charge, source.launch.massEV);
        EXPECT_NEAR(c[1][2], -1.5 * charge, 1e-14);
        EXPECT_NEAR(c[3][0], 1.5 * charge, 1e-14);
        auto inverse  = identity();
        inverse[1][2] = -c[1][2];
        inverse[3][0] = -c[3][0];
        source.matrix = product(product(inverse, rotations()), c);
        Matrix diagonal{};
        diagonal[0][0] = diagonal[1][1] = 2e-6;
        diagonal[2][2] = diagonal[3][3] = 7e-6;
        const auto actual               = matched(source, {2e-6, 7e-6});
        expectMatrix(actual.covariance, product(product(inverse, diagonal), transpose(inverse)));
        expectMatrix(
                product(product(source.matrix, actual.covariance), transpose(source.matrix)),
                actual.covariance);
        source.normalMagneticField = 0;
        EXPECT_THROW(MatchedCovariance::compute(source, {2e-6, 7e-6}), OpalException);
    }
}

TEST_F(MatchedCovarianceTest, RejectsUnsupportedMapsAndInvalidInputs) {
    const auto valid = optics(rotations());
    for (double invalid :
         {0., -1., std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::quiet_NaN()}) {
        EXPECT_THROW(MatchedCovariance::compute(valid, {invalid, 1e-6}), OpalException);
        EXPECT_THROW(MatchedCovariance::compute(valid, {1e-6, invalid}), OpalException);
    }
    for (unsigned change = 0; change < 13; ++change) {
        auto o = valid;
        if (change == 0) o.residual[0] = std::numeric_limits<double>::quiet_NaN();
        if (change == 1) o.matrix = identity();
        if (change == 2) o.matrix = rotations(.2, .2);
        if (change == 3) o.matrix = rotations(.2, .2000001);
        if (change == 4) o.matrix = rotations(.5 - 1e-8, .2);
        if (change == 5)
            for (auto& row : o.matrix)
                for (double& x : row)
                    x *= 1.001;
        if (change == 6) o.matrix[0][2] = .01;  // Spectrum stable, eigenbasis not symplectic.
        if (change == 7) o.scales[1] = 0;
        if (change == 8) o.normalMagneticField = std::numeric_limits<double>::infinity();
        if (change == 9) o.launch.massEV = -1;
        if (change == 10) o.relativeEnergyDrift = 1e-7;
        if (change == 11) o.relativeEnergyDrift = std::numeric_limits<double>::quiet_NaN();
        if (change == 12) o.launch.chargeE = 0;
        EXPECT_THROW(MatchedCovariance::compute(o, {2e-6, 7e-6}), OpalException) << change;
    }
    MatchedCovariance::Settings settings;
    settings.periodicTolerance = 0;
    EXPECT_THROW(MatchedCovariance::compute(valid, {2e-6, 7e-6}, settings), OpalException);
}
