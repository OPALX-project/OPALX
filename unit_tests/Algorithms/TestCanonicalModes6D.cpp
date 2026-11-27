/** @file TestCanonicalModes6D.cpp
 * @brief Independent canonical fixtures; no OPAL input files or particle runs.
 */
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include "Algorithms/CanonicalModes6D.h"
#include "Ippl.h"
#include "Utilities/OpalException.h"

namespace {
    namespace modes = opalx::maps::modes6d;
    using Matrix    = modes::Matrix;
    using Eps       = modes::Emittances;
    using Status    = modes::Status;
    Matrix identity() {
        Matrix a{};
        for (unsigned i = 0; i < 6; ++i)
            a[i][i] = 1;
        return a;
    }
    Matrix transpose(const Matrix& a) {
        Matrix b{};
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                b[i][j] = a[j][i];
        return b;
    }
    Matrix product(const Matrix& a, const Matrix& b) {
        Matrix c{};
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                for (unsigned k = 0; k < 6; ++k)
                    c[i][j] += a[i][k] * b[k][j];
        return c;
    }
    Matrix form() {
        Matrix j{};
        for (unsigned k = 0; k < 6; k += 2) {
            j[k][k + 1] = 1;
            j[k + 1][k] = -1;
        }
        return j;
    }
    Matrix inverseCanonical(const Matrix& a) {
        auto b = product(product(form(), transpose(a)), form());
        for (auto& row : b)
            for (auto& x : row)
                x = -x;
        return b;
    }
    Matrix rotation(Eps phase = {0.13, 0.27, 0.39}, std::array<int, 3> sign = {1, 1, 1}) {
        Matrix r{};
        for (unsigned k = 0; k < 3; ++k) {
            double angle    = 2 * std::acos(-1.0) * phase[k];
            r[2 * k][2 * k] = r[2 * k + 1][2 * k + 1] = std::cos(angle);
            r[2 * k][2 * k + 1]                       = sign[k] * std::sin(angle);
            r[2 * k + 1][2 * k]                       = -sign[k] * std::sin(angle);
        }
        return r;
    }
    Matrix planeMix(double angle, unsigned a = 0, unsigned b = 1) {
        auto s = identity();
        for (unsigned side = 0; side < 2; ++side) {
            s[2 * a + side][2 * a + side] = s[2 * b + side][2 * b + side] = std::cos(angle);
            s[2 * a + side][2 * b + side]                                 = std::sin(angle);
            s[2 * b + side][2 * a + side]                                 = -std::sin(angle);
        }
        return s;
    }
    Matrix coupledBasis() {
        auto shear  = identity();
        shear[1][0] = 0.3;
        shear[1][2] = shear[3][0] = -0.16;
        shear[3][4] = shear[5][2] = 0.23;
        shear[5][4]               = -0.1;
        auto squeeze              = identity();
        for (unsigned k = 0; k < 3; ++k) {
            squeeze[2 * k][2 * k]         = 0.7 + 0.6 * k;
            squeeze[2 * k + 1][2 * k + 1] = 1 / squeeze[2 * k][2 * k];
        }
        return product(
                product(product(planeMix(0.33, 0, 2), planeMix(-0.42, 1, 2)), shear), squeeze);
    }
    Matrix covariance(const Matrix& s, Eps eps) {
        Matrix d{};
        for (unsigned i = 0; i < 6; ++i)
            d[i][i] = eps[i / 2];
        return product(product(s, d), transpose(s));
    }
    Matrix conjugate(const Matrix& s, const Matrix& r) {
        return product(product(s, r), inverseCanonical(s));
    }
    void expectMatrix(const Matrix& a, const Matrix& b, double tol = 3e-12) {
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                EXPECT_NEAR(a[i][j], b[i][j], tol * std::sqrt(std::abs(b[i][i] * b[j][j]) + 1e-100))
                        << i << ',' << j;
    }
    void expectCanonical(const Matrix& s) {
        const auto actual   = product(product(transpose(s), form()), s);
        const auto expected = form();
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                EXPECT_NEAR(actual[i][j], expected[i][j], 1e-11) << i << ',' << j;
    }
    // Independent spectral invariants: tr[(J Sigma)^(2p)]=-2 sum_k eps_k^(2p).
    void expectActions(const Matrix& sigma, Eps eps) {
        const auto a2 = product(product(form(), sigma), product(form(), sigma));
        auto power    = identity();
        for (unsigned p = 1; p <= 3; ++p) {
            power        = product(power, a2);
            double trace = 0, expected = 0;
            for (unsigned i = 0; i < 6; ++i)
                trace += power[i][i];
            for (double e : eps)
                expected += 2 * (p % 2 ? -1 : 1) * std::pow(e, 2 * p);
            EXPECT_NEAR(trace, expected, 2e-10 * std::abs(expected));
        }
    }
    tps::PhaseSpaceChart chart() {
        tps::PhaseSpaceChart c;
        c.name           = "periodic canonical section";
        c.frame.name     = "section axes";
        c.convention     = tps::CoordinateConvention::canonical;
        c.canonicalPairs = {{0, 1}, {2, 3}, {4, 5}};
        for (unsigned i = 0; i < 6; ++i)
            c.coordinates.push_back({"c" + std::to_string(i), i % 2 ? "p/(mc)" : "m", 1});
        c.reference.resize(6);
        c.referenceTime = 1.2;
        c.parameters    = {{"species_rest_energy", "eV", 8e8}, {"canonical_kappa_Bx", "1/m", 0}};
        return c;
    }
    tps::Map<double> tpsMap(
            const Matrix& physical, tps::PhaseSpaceChart in, tps::PhaseSpaceChart out,
            int order = 1, double affine = 0) {
        const int saved = tps::Tps<double>::getGlobalTruncOrder();
        tps::Tps<double>::setGlobalTruncOrder(tps::Tps<double>::EXACT);
        tps::VpsMap<double> p(6, 6);
        for (unsigned i = 0; i < 6; ++i) {
            p[i]    = tps::Tps<double>(order, 6);
            p[i][0] = affine;
            for (unsigned j = 0; j < 6; ++j)
                p[i][j + 1] = physical[i][j] * in.coordinates[j].scale / out.coordinates[i].scale;
        }
        auto result = tps::Map<double>(tps::MapContext(in, out, order, {"periodic"}), p);
        tps::Tps<double>::setGlobalTruncOrder(saved);
        return result;
    }
    class CanonicalModes6DTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }
        static void TearDownTestSuite() { ippl::finalize(); }
        void SetUp() override { order = tps::Tps<double>::getGlobalTruncOrder(); }
        void TearDown() override { tps::Tps<double>::setGlobalTruncOrder(order); }
        int order;
    };
}  // namespace

TEST_F(CanonicalModes6DTest, IndependentRotationsBothKreinSignsAndCovariance) {
    const auto a = modes::analyze(rotation({0.13, 0.27, 0.39}, {1, -1, 1}));
    ASSERT_EQ(a.status, Status::Stable);
    ASSERT_EQ(a.modes.size(), 3u);
    const Eps phase{0.13, 0.27, 0.39}, eps{2e-6, 7e-6, 1.1e-5};
    for (unsigned k = 0; k < 3; ++k) {
        EXPECT_NEAR(a.modes[k].phase, phase[k], 1e-14);
        EXPECT_EQ(a.modes[k].kreinSign, k == 1 ? -1 : 1);
    }
    expectCanonical(a.basis);
    auto sigma = modes::invariantCovariance(a, eps);
    expectMatrix(sigma, covariance(identity(), eps));
    expectActions(sigma, eps);
    EXPECT_LT(a.maximumEigenResidual, 1e-14);
    EXPECT_LT(a.normalizationResidual, 1e-13);
}

TEST_F(CanonicalModes6DTest, CoupledSymplecticSimilarityHasKnownInvariantCovariance) {
    const auto s = coupledBasis();
    expectCanonical(s);
    const auto m = conjugate(s, rotation({0.13, 0.27, 0.39}, {-1, 1, -1}));
    const auto a = modes::analyze(m);
    ASSERT_EQ(a.status, Status::Stable);
    const Eps eps{1.3e-6, 6.5e-6, 2.1e-5};
    auto sigma = modes::invariantCovariance(a, eps);
    expectMatrix(sigma, covariance(s, eps));
    expectMatrix(product(product(m, sigma), transpose(m)), sigma);
    expectActions(sigma, eps);
    expectCanonical(a.basis);
}

TEST_F(CanonicalModes6DTest, NumericalScalesCarryCanonicalFormAndExtremeCommonScale) {
    auto m     = conjugate(coupledBasis(), rotation());
    auto plain = modes::analyze(m);
    ASSERT_EQ(plain.status, Status::Stable);
    const Eps eps{1e-6, 2e-6, 5e-6};
    for (double common : {1., 1e-150, 1e150}) {
        modes::Settings settings;
        settings.scales = {common * 0.3, common * 2,   common * 4,
                           common * 0.7, common * 1.1, common * 0.8};
        const auto a    = modes::analyze(m, settings);
        ASSERT_EQ(a.status, Status::Stable);
        expectMatrix(
                modes::invariantCovariance(a, eps, settings),
                modes::invariantCovariance(plain, eps));
        const auto w = modes::williamson(covariance(coupledBasis(), eps), settings);
        for (unsigned k = 0; k < 3; ++k)
            EXPECT_NEAR(w.emittances[k], eps[k], eps[k] * 2e-12);
        expectCanonical(w.basis);
    }
}

TEST_F(CanonicalModes6DTest, HyperbolicAndComplexQuartetInstability) {
    auto m  = rotation();
    m[0][0] = 2;
    m[0][1] = m[1][0] = 0;
    m[1][1]           = 0.5;
    const auto a      = modes::analyze(m);
    EXPECT_EQ(a.status, Status::Unstable);
    EXPECT_TRUE(a.modes.empty());
    // q->exp(g)*O q, p->exp(-g)*O p gives a real symplectic complex quartet.
    m = rotation();
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j)
            m[i][j] = 0;
    for (unsigned side = 0; side < 2; ++side) {
        const double gain = std::exp(side ? -0.1 : 0.1);
        m[side][side] = m[2 + side][2 + side] = gain * std::cos(0.7);
        m[side][2 + side]                     = gain * std::sin(0.7);
        m[2 + side][side]                     = -gain * std::sin(0.7);
    }
    EXPECT_EQ(modes::analyze(m).status, Status::Unstable);
}

TEST_F(CanonicalModes6DTest, RejectsUnitCircleNonCanonicalSimilarityAndDamping) {
    auto m = rotation();
    // Mixing q and p of different pairs without the canonical partner change.
    auto s = identity(), inv = identity();
    s[0][2]      = 0.3;
    inv[0][2]    = -0.3;
    m            = product(product(s, m), inv);
    const auto a = modes::analyze(m);
    EXPECT_EQ(a.status, Status::NonSymplectic);
    EXPECT_TRUE(a.modes.empty());
    EXPECT_GT(a.symplecticResidual, 1e-3);
    m = rotation();
    for (auto& row : m)
        for (auto& x : row)
            x *= 0.9;
    EXPECT_EQ(modes::analyze(m).status, Status::NonSymplectic);
    modes::Settings loose;
    loose.symplecticTolerance = 0.2;
    EXPECT_EQ(modes::analyze(m, loose).status, Status::NonUnitCircle);
}

TEST_F(CanonicalModes6DTest, MarginalIntegerHalfIntegerAndJordanAreNotStableModes) {
    for (double phase : {0., 1e-8, 0.5 - 1e-8, 0.5}) {
        const auto a = modes::analyze(rotation({0.13, 0.27, phase}));
        EXPECT_EQ(a.status, Status::Marginal);
        EXPECT_TRUE(a.modes.empty());
        EXPECT_THROW(modes::invariantCovariance(a, {1e-6, 2e-6, 3e-6}), OpalException);
    }
    auto jordan  = identity();
    jordan[0][1] = 1;
    EXPECT_EQ(modes::analyze(jordan).status, Status::Marginal);
    EXPECT_EQ(modes::analyze(identity()).status, Status::Marginal);
}

TEST_F(CanonicalModes6DTest, UnresolvedOscillatoryModesAndIllConditionedBasis) {
    EXPECT_EQ(modes::analyze(rotation({0.13, 0.13, 0.39})).status, Status::Degenerate);
    EXPECT_EQ(modes::analyze(rotation({0.13, 0.13 + 1e-8, 0.39})).status, Status::Degenerate);
    auto squeeze  = identity();
    squeeze[0][0] = 1e6;
    squeeze[1][1] = 1e-6;
    const auto m  = conjugate(squeeze, rotation());
    EXPECT_EQ(modes::analyze(m).status, Status::IllConditioned);
    modes::Settings equilibrated;
    equilibrated.scales = {1e6, 1e-6, 1, 1, 1, 1};
    const auto a        = modes::analyze(m, equilibrated);
    ASSERT_EQ(a.status, Status::Stable);
    expectMatrix(
            modes::invariantCovariance(a, {1e-6, 2e-6, 3e-6}, equilibrated),
            covariance(squeeze, {1e-6, 2e-6, 3e-6}));
}

TEST_F(CanonicalModes6DTest, WilliamsonReconstructsCoupledCovarianceAndRestoresActions) {
    const Eps original{8e-6, 2e-6, 4e-6}, target{3e-6, 7e-6, 1.5e-5};
    const auto sigma = covariance(coupledBasis(), original);
    auto w           = modes::williamson(sigma);
    expectCanonical(w.basis);
    expectMatrix(covariance(w.basis, w.emittances), sigma);
    Eps sorted = original;
    std::sort(sorted.begin(), sorted.end());
    for (unsigned k = 0; k < 3; ++k)
        EXPECT_NEAR(w.emittances[k], sorted[k], sorted[k] * 2e-12);
    auto restored = modes::withEmittances(w, target);
    expectActions(restored, target);
    auto reset = modes::williamson(restored);
    for (unsigned k = 0; k < 3; ++k)
        EXPECT_NEAR(reset.emittances[k], target[k], target[k] * 2e-12);
    EXPECT_LT(w.reconstructionResidual, 1e-12);
}

TEST_F(CanonicalModes6DTest, EqualActionsAndCovarianceMixing) {
    const Eps eps{2e-6, 2e-6, 2e-6};
    auto sigma = covariance(coupledBasis(), eps);
    auto w     = modes::williamson(sigma);
    expectCanonical(w.basis);
    expectMatrix(covariance(w.basis, w.emittances), sigma);
    for (double e : w.emittances)
        EXPECT_NEAR(e, 2e-6, 4e-18);
    auto other = covariance(identity(), eps);
    for (unsigned i = 0; i < 6; ++i)
        for (unsigned j = 0; j < 6; ++j)
            sigma[i][j] = (sigma[i][j] + other[i][j]) / 2;
    w = modes::williamson(sigma);
    EXPECT_GT(w.emittances[2], 2.01e-6);
    expectActions(modes::withEmittances(w, eps), eps);
}

TEST_F(CanonicalModes6DTest, VeryDifferentPhysicalCovarianceScales) {
    auto s  = identity();
    s[0][0] = 1e-20;
    s[1][1] = 1e20;
    s[4][4] = 1e15;
    s[5][5] = 1e-15;
    const Eps eps{2e-6, 3e-6, 5e-6};
    auto sigma = covariance(s, eps);
    modes::Settings settings;
    settings.scales = {1e-20, 1e20, 1, 1, 1e15, 1e-15};
    auto w          = modes::williamson(sigma, settings);
    expectMatrix(covariance(w.basis, w.emittances), sigma);
    for (unsigned k = 0; k < 3; ++k)
        EXPECT_NEAR(w.emittances[k], eps[k], eps[k] * 1e-13);
    for (double factor : {1e-200, 1e200}) {
        auto tiny   = covariance(identity(), {factor, factor * 2, factor * 3});
        auto scaled = modes::williamson(tiny);
        for (unsigned k = 0; k < 3; ++k)
            EXPECT_NEAR(scaled.emittances[k], factor * (k + 1), factor * (k + 1) * 1e-13);
    }
}

TEST_F(CanonicalModes6DTest, TracksSubspacesAcrossPhaseOrderingAndRejectsAmbiguity) {
    auto previous = modes::analyze(rotation({0.13, 0.27, 0.39}));
    auto current  = modes::analyze(rotation({0.29, 0.15, 0.39}));
    ASSERT_EQ(current.status, Status::Stable);
    auto tracked = modes::trackModes(current, previous);
    EXPECT_NEAR(tracked.modes[0].phase, 0.29, 1e-14);
    EXPECT_NEAR(tracked.modes[1].phase, 0.15, 1e-14);
    const Eps eps{2e-6, 5e-6, 9e-6};
    expectMatrix(modes::invariantCovariance(tracked, eps), covariance(identity(), eps));
    auto mixed = modes::analyze(conjugate(planeMix(std::acos(-1.) / 4), rotation()));
    ASSERT_EQ(mixed.status, Status::Stable);
    EXPECT_THROW(modes::trackModes(mixed, previous), OpalException);
    auto reversed = modes::analyze(rotation({0.13, 0.27, 0.39}, {1, -1, 1}));
    EXPECT_THROW(modes::trackModes(reversed, previous), OpalException);
}

TEST_F(CanonicalModes6DTest, TpsCanonicalAdapterRestoresLocalScalesAndPreservesAmbientOrder) {
    auto m = conjugate(coupledBasis(), rotation());
    auto c = chart();
    for (unsigned i = 0; i < 6; ++i)
        c.coordinates[i].scale = 0.3 + 0.4 * i;
    auto map = tpsMap(m, c, c);
    tps::Tps<double>::setGlobalTruncOrder(0);
    const auto a = modes::analyze(map);
    EXPECT_EQ(tps::Tps<double>::getGlobalTruncOrder(), 0);
    ASSERT_EQ(a.status, Status::Stable);
    expectMatrix(a.map, m);
    expectMatrix(
            modes::invariantCovariance(a, {1e-6, 2e-6, 4e-6}),
            covariance(coupledBasis(), {1e-6, 2e-6, 4e-6}));
}

TEST_F(CanonicalModes6DTest, TpsAdapterRejectsBoundaryTimeGaugeSpeciesAndInvalidMaps) {
    const auto c      = chart();
    auto out          = c;
    out.referenceTime = *c.referenceTime + 1e-6;
    EXPECT_THROW(modes::analyze(tpsMap(rotation(), c, out)), OpalException);
    out = c;
    out.parameters[0].value *= 2;
    EXPECT_THROW(modes::analyze(tpsMap(rotation(), c, out)), OpalException);
    out                     = c;
    out.parameters[1].value = 1;
    EXPECT_THROW(modes::analyze(tpsMap(rotation(), c, out)), OpalException);
    out              = c;
    out.reference[0] = 0.01;
    EXPECT_THROW(modes::analyze(tpsMap(rotation(), c, out)), OpalException);
    auto mechanical       = c;
    mechanical.convention = tps::CoordinateConvention::mechanical;
    mechanical.canonicalPairs.clear();
    EXPECT_THROW(modes::analyze(tpsMap(rotation(), mechanical, mechanical)), OpalException);
    auto wrong                = c;
    wrong.coordinates[1].unit = "1";
    EXPECT_THROW(modes::analyze(tpsMap(rotation(), wrong, wrong)), OpalException);
    EXPECT_THROW(modes::analyze(tpsMap(rotation(), c, c, 1, 0.01)), OpalException);
    EXPECT_THROW(modes::analyze(tpsMap(rotation(), c, c, 2)), OpalException);
    auto nonfinite  = rotation();
    nonfinite[0][0] = std::numeric_limits<double>::infinity();
    EXPECT_THROW(modes::analyze(tpsMap(nonfinite, c, c)), OpalException);
}

TEST_F(CanonicalModes6DTest, InvalidInputsSettingsAndPublicSnapshotsFailExplicitly) {
    EXPECT_EQ(modes::analyze(Matrix{}).status, Status::NonSymplectic);
    auto m  = rotation();
    m[0][0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(modes::analyze(m), OpalException);
    modes::Settings bad;
    bad.scales[0] = 0;
    EXPECT_THROW(modes::analyze(rotation(), bad), OpalException);
    bad                   = {};
    bad.minimumBasisRcond = -1;
    EXPECT_THROW(modes::analyze(rotation(), bad), OpalException);
    bad                      = {};
    bad.minimumAssignmentGap = std::numeric_limits<double>::infinity();
    EXPECT_THROW(modes::williamson(identity(), bad), OpalException);
    auto sigma  = identity();
    sigma[0][0] = -1;
    EXPECT_THROW(modes::williamson(sigma), OpalException);
    sigma       = identity();
    sigma[0][1] = 0.1;
    EXPECT_THROW(modes::williamson(sigma), OpalException);
    sigma       = identity();
    sigma[0][1] = sigma[1][0] = 2;
    EXPECT_THROW(modes::williamson(sigma), OpalException);
    const auto a = modes::analyze(rotation());
    EXPECT_THROW(modes::invariantCovariance(a, {0, 1, 2}), OpalException);
    auto forged = a;
    for (auto& row : forged.basis)
        row[0] *= 2;
    EXPECT_THROW(modes::invariantCovariance(forged, {1, 2, 3}), OpalException);
    forged = a;
    forged.map[0][0] += 0.1;
    EXPECT_THROW(modes::invariantCovariance(forged, {1, 2, 3}), OpalException);
    forged                    = a;
    forged.modes[0].kreinSign = 0;
    EXPECT_THROW(modes::trackModes(forged, a), OpalException);
    forged = a;
    forged.modes[0].phase += 0.01;
    EXPECT_THROW(modes::trackModes(forged, a), OpalException);
    auto w        = modes::williamson(identity());
    w.basis[0][0] = 0;
    EXPECT_THROW(modes::withEmittances(w, {1, 2, 3}), OpalException);
}

TEST_F(CanonicalModes6DTest, AllRanksObtainIdenticalModesNormalizationAndCovariance) {
    auto a = modes::analyze(conjugate(coupledBasis(), rotation()));
    ASSERT_EQ(a.status, Status::Stable);
    const auto sigma = modes::invariantCovariance(a, {1e-6, 2e-6, 4e-6});
    const auto w     = modes::williamson(sigma);
    std::array<double, 81> values{}, minimum{}, maximum{};
    unsigned index = 0;
    for (const auto& row : a.basis)
        for (double x : row)
            values[index++] = x;
    for (const auto& row : sigma)
        for (double x : row)
            values[index++] = x;
    for (unsigned k = 0; k < 3; ++k) {
        values[index++] = a.modes[k].phase;
        values[index++] = a.modes[k].kreinSign;
        values[index++] = w.emittances[k];
    }
    MPI_Allreduce(
            values.data(), minimum.data(), values.size(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(
            values.data(), maximum.data(), values.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    EXPECT_EQ(minimum, maximum);
}

TEST_F(CanonicalModes6DTest, FamilyOfFullyCoupledMapsAndCovariancesHasKnownSolutions) {
    // Deterministic parameter sweep: independent symplectic factors, all three planes,
    // distinct action ratios and both rotation orientations. No fitted expected values.
    for (unsigned sample = 0; sample < 24; ++sample) {
        SCOPED_TRACE(sample);
        const double f = (sample + 1.) / 24;
        auto squeeze = identity(), shear = identity();
        for (unsigned k = 0; k < 3; ++k) {
            squeeze[2 * k][2 * k]         = std::exp((k - 1.) * 2 * f);
            squeeze[2 * k + 1][2 * k + 1] = 1 / squeeze[2 * k][2 * k];
            for (unsigned l = 0; l < 3; ++l)
                shear[2 * k + 1][2 * l] = (k + l + 1.) * 0.07 * f;
        }
        const auto s =
                product(product(product(planeMix(0.6 * f, 0, 1), planeMix(-0.7 * f, 0, 2)), shear),
                        squeeze);
        const Eps phase{0.08 + 0.04 * f, 0.22 + 0.03 * f, 0.39 - 0.02 * f};
        const Eps eps{2e-6 * (1 + f), 7e-6 * (1 + 0.5 * f), 1e-5 * (2 - f)};
        const auto m = conjugate(s, rotation(phase, {sample % 2 ? -1 : 1, 1, sample % 3 ? -1 : 1}));
        modes::Settings settings;
        settings.scales = {0.8 + f, 2 + f, 0.4 + f, 1.1 + f, 2.3 - f, 1.7 - f};
        const auto a    = modes::analyze(m, settings);
        ASSERT_EQ(a.status, Status::Stable);
        const auto sigma = modes::invariantCovariance(a, eps, settings);
        expectMatrix(sigma, covariance(s, eps), 1e-11);
        expectActions(sigma, eps);
        const auto w = modes::williamson(sigma, settings);
        auto sorted  = eps;
        std::sort(sorted.begin(), sorted.end());
        for (unsigned k = 0; k < 3; ++k)
            EXPECT_NEAR(w.emittances[k], sorted[k], sorted[k] * 1e-11);
        expectCanonical(w.basis);
    }
}

TEST_F(CanonicalModes6DTest, FrozenAnalyticSpaceChargeKicksPreserveWilliamsonActions) {
    using namespace opalx::maps;
    SpaceChargeReference ref;
    ref.massEV           = 9.38e8;
    ref.chargeE          = 1;
    ref.orbit.momentum   = {0.1, -0.05, 2};
    ref.magneticFieldLab = {0.2, -0.1, 0.4};
    const SpaceChargeCoordinates coordinates(ref);
    const Eps eps{2e-6, 5e-6, 9e-6};
    const auto sigma = covariance(coupledBasis(), eps);
    for (auto profile :
         {opalx::spacecharge::analytic3d::Profile::UniformEllipsoid,
          opalx::spacecharge::analytic3d::Profile::GaussianRms}) {
        const auto kick = makeSpaceChargeKick(
                coordinates, sigma, SpaceChargeChart::CommonTimeCanonical, profile, 1e-10, 0.03,
                "kick");
        const auto transported = transportSpaceChargeCovariance(kick.polynomial, sigma);
        const auto w           = modes::williamson(transported);
        for (unsigned k = 0; k < 3; ++k)
            EXPECT_NEAR(w.emittances[k], eps[k], eps[k] * 3e-12);
        expectActions(transported, eps);
    }
}
