/**
 * @file TestBeamMoments.cpp
 * @brief Unit tests for opalx::moments::compute() and BeamMoments.
 *
 * Compares every group against a two-pass host reference in long double, checks the request
 * contract (checked accessors, dependencies, request independence of shared sums), the
 * degenerate conventions, the union over several containers, and the rest-frame temperature of a
 * boosted set. Particles are dealt round-robin over the ranks, so `mpirun -n 2` also checks the
 * MPI reductions; CTest runs one rank.
 */

#include <gtest/gtest.h>
#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <random>
#include <vector>

#include "Algorithms/PartData.h"
#include "Ippl.h"
#include "Moments/BeamMoments.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Structure/DataSink.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

namespace {

    using opalx::moments::BeamMoments;
    using opalx::moments::compute;
    using opalx::moments::Moment;
    using opalx::moments::Request;
    using PC_t = ParticleContainer<double, 3>;
    using Vec3 = ippl::Vector<double, 3>;

    static_assert((Moment::Bounds | Moment::Spread).contains(Moment::Spread));
    static_assert(!Request(Moment::Bounds).contains(Moment::MeanPosition));
    static_assert(
            opalx::moments::allMoments.without(Moment::Polarization).contains(Moment::Spread));
    static_assert((Request(Moment::Bounds) | Moment::Spread) == (Moment::Spread | Moment::Bounds));

    struct Particle {
        Vec3 r;
        Vec3 p;
        double massGeV;
        Vec3 pol{0.0};
    };

    /// Two-pass reference in long double, with the same definitions as compute().
    struct Reference {
        std::array<long double, 3> minR{}, maxR{}, meanR{}, meanP{};
        long double meanGamma = 0, meanEkin = 0, stdEkin = 0;
        std::array<std::array<long double, 6>, 6> cov{};
        std::array<long double, 3> temperature{}, meanPol{}, rmsPol{};
    };

    Reference reference(const std::vector<Particle>& set) {
        Reference ref;
        const long double n = static_cast<long double>(std::max<std::size_t>(set.size(), 1));
        for (int d = 0; d < 3; ++d) {
            ref.minR[d] = set.empty() ? 0 : set[0].r[d];
            ref.maxR[d] = ref.minR[d];
        }
        auto gammaOf = [](const Vec3& p) {
            return std::sqrt(
                    1.0L + (long double)p[0] * p[0] + (long double)p[1] * p[1]
                    + (long double)p[2] * p[2]);
        };
        for (const Particle& q : set) {
            const long double gamma = gammaOf(q.p);
            for (int d = 0; d < 3; ++d) {
                ref.minR[d] = std::min<long double>(ref.minR[d], q.r[d]);
                ref.maxR[d] = std::max<long double>(ref.maxR[d], q.r[d]);
                ref.meanR[d] += q.r[d] / n;
                ref.meanP[d] += q.p[d] / n;
                ref.meanPol[d] += q.pol[d] / n;
            }
            ref.meanGamma += gamma / n;
            ref.meanEkin += (gamma - 1) * q.massGeV * Units::GeV2MeV / n;
        }
        std::array<long double, 3> betaC{};
        long double betaC2 = 0;
        for (int d = 0; d < 3; ++d) {
            betaC[d] = ref.meanP[d] / ref.meanGamma;
            betaC2 += betaC[d] * betaC[d];
        }
        const long double gammaC = 1 / std::sqrt(1 - betaC2);
        for (const Particle& q : set) {
            const long double gamma = gammaOf(q.p);
            long double c[6];
            for (int d = 0; d < 3; ++d) {
                c[2 * d]     = q.r[d] - ref.meanR[d];
                c[2 * d + 1] = q.p[d] - ref.meanP[d];
            }
            for (int a = 0; a < 6; ++a) {
                for (int b = 0; b < 6; ++b) {
                    ref.cov[a][b] += c[a] * c[b] / n;
                }
            }
            const long double ekin = (gamma - 1) * q.massGeV * Units::GeV2MeV;
            ref.stdEkin += (ekin - ref.meanEkin) * (ekin - ref.meanEkin) / n;
            const long double betaDotDu = betaC[0] * c[1] + betaC[1] * c[3] + betaC[2] * c[5];
            const long double shift =
                    gammaC * gammaC / (gammaC + 1) * betaDotDu - gammaC * (gamma - ref.meanGamma);
            long double uRest[3], u2 = 0;
            for (int d = 0; d < 3; ++d) {
                uRest[d] = c[2 * d + 1] + shift * betaC[d];
                u2 += uRest[d] * uRest[d];
            }
            for (int d = 0; d < 3; ++d) {
                ref.temperature[d] += q.massGeV * uRest[d] * uRest[d] / std::sqrt(1 + u2)
                                      * Units::GeV2eV / n / Physics::kB;
                ref.rmsPol[d] += (long double)q.pol[d] * q.pol[d] / n;
            }
        }
        ref.stdEkin = std::sqrt(ref.stdEkin);
        for (int d = 0; d < 3; ++d) {
            ref.rmsPol[d] =
                    std::sqrt(std::max(ref.rmsPol[d] - ref.meanPol[d] * ref.meanPol[d], 0.0L));
        }
        return ref;
    }

    void expectRelative(double value, long double expected, double tolerance, const char* what) {
        const long double scale = std::max<long double>(std::fabs(expected), 1e-300L);
        EXPECT_LE(std::fabs(value - expected) / scale, tolerance)
                << what << ": " << value << " vs " << static_cast<double>(expected);
    }

    std::vector<Particle> correlatedSet(std::size_t n, unsigned seed, double massGeV) {
        std::mt19937_64 engine(seed);
        std::normal_distribution<double> normal;
        std::vector<Particle> set(n);
        for (Particle& q : set) {
            const double x = normal(engine), y = normal(engine), z = normal(engine);
            q.r = Vec3(1e-3 * x + 2e-4, 2e-3 * y - 1e-4, 5e-4 * z + 3e-3);
            q.p =
                    Vec3(0.02 * x + 0.01 * normal(engine), 0.03 * normal(engine) - 0.002,
                         5.0 + 0.4 * z + 0.05 * normal(engine));
            q.massGeV = massGeV;
            q.pol = Vec3(0.3 + 0.1 * normal(engine), -0.2 + 0.02 * normal(engine), 0.5 + 0.05 * x);
        }
        return set;
    }

    class BeamMomentsTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }

        void SetUp() override { Options::useQMAttributes = false; }
        void TearDown() override { Options::useQMAttributes = false; }

        /// Indices of @p n particles that this rank holds (round-robin over the ranks).
        static std::vector<std::size_t> localIndices(std::size_t n) {
            std::vector<std::size_t> indices;
            for (std::size_t i = static_cast<std::size_t>(ippl::Comm->rank()); i < n;
                 i += static_cast<std::size_t>(ippl::Comm->size())) {
                indices.push_back(i);
            }
            return indices;
        }

        /// Container on an 8^3 mesh over [-4, 4]^3 holding @p set, dealt over the ranks.
        std::shared_ptr<PC_t> makeContainer(
                const std::vector<Particle>& set, bool spin = false, double scalarMassGeV = 1.0) {
            ippl::Vector<int, 3> nr        = 8;
            ippl::Vector<double, 3> origin = -4.0;
            ippl::Vector<double, 3> hr     = 1.0;
            ippl::NDIndex<3> domain;
            for (unsigned i = 0; i < 3; i++) {
                domain[i] = ippl::Index(nr[i]);
            }
            Mesh_t<3> mesh(domain, hr, origin);
            FieldLayout_t<3> layout(MPI_COMM_WORLD, domain, {true, true, true}, true);
            auto pc = std::make_shared<PC_t>(
                    mesh, layout, spin, PC_t::LayoutType::Spatial, 0.0, ippl::BC::NO);
            pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());

            const auto local = localIndices(set.size());
            pc->createParticles(local.size());
            auto r = pc->R.getHostMirror();
            auto p = pc->P.getHostMirror();
            for (std::size_t i = 0; i < local.size(); ++i) {
                r(i) = set[local[i]].r;
                p(i) = set[local[i]].p;
            }
            Kokkos::deep_copy(pc->R.getView(), r);
            Kokkos::deep_copy(pc->P.getView(), p);
            if (spin) {
                auto pol = pc->Pol.getHostMirror();
                for (std::size_t i = 0; i < local.size(); ++i) {
                    pol(i) = set[local[i]].pol;
                }
                Kokkos::deep_copy(pc->Pol.getView(), pol);
            }
            pc->setM(scalarMassGeV);
            Kokkos::fence();
            return pc;
        }
    };

    void expectMatchesReference(const BeamMoments& m, const Reference& ref, double tolerance) {
        for (int d = 0; d < 3; ++d) {
            EXPECT_EQ(m.minR()[d], static_cast<double>(ref.minR[d]));
            EXPECT_EQ(m.maxR()[d], static_cast<double>(ref.maxR[d]));
            expectRelative(m.meanR()[d], ref.meanR[d], tolerance, "meanR");
            expectRelative(m.meanP()[d], ref.meanP[d], tolerance, "meanP");
            expectRelative(m.restFrameTemperature()[d], ref.temperature[d], 1e-10, "temperature");
            expectRelative(m.rmsR()[d], std::sqrt(ref.cov[2 * d][2 * d]), tolerance, "rmsR");
            expectRelative(
                    m.rmsP()[d], std::sqrt(ref.cov[2 * d + 1][2 * d + 1]), tolerance, "rmsP");
            const long double sx = std::sqrt(ref.cov[2 * d][2 * d]);
            const long double sp = std::sqrt(ref.cov[2 * d + 1][2 * d + 1]);
            const long double xp = ref.cov[2 * d][2 * d + 1];
            expectRelative(m.correlationRP()[d], xp / (sx * sp), tolerance, "correlationRP");
            const long double emittance = std::sqrt(sx * sx * sp * sp - xp * xp);
            expectRelative(m.normalizedEmittance()[d], emittance, 1e-10, "emittance");
            expectRelative(
                    m.geometricEmittance()[d],
                    emittance / std::sqrt(ref.meanGamma * ref.meanGamma - 1), 1e-10,
                    "geometricEmittance");
        }
        for (unsigned a = 0; a < 6; ++a) {
            for (unsigned b = 0; b < 6; ++b) {
                expectRelative(m.covariance(a, b), ref.cov[a][b], 1e-11, "covariance");
                EXPECT_EQ(m.covariance(a, b), m.covariance(b, a));
            }
        }
        expectRelative(m.meanGamma(), ref.meanGamma, tolerance, "meanGamma");
        expectRelative(
                m.meanGammaZ(), std::sqrt(ref.meanP[2] * ref.meanP[2] + 1), tolerance,
                "meanGammaZ");
        expectRelative(m.meanKineticEnergy(), ref.meanEkin, tolerance, "meanEkin");
        expectRelative(m.stdKineticEnergy(), ref.stdEkin, 1e-10, "stdEkin");
    }

}  // namespace

TEST_F(BeamMomentsTest, AllGroupsMatchHostReference) {
    // Enough particles for many work items, so GPU builds also reduce over several blocks.
    const auto set = correlatedSet(20000, 7, 0.938);
    auto pc        = makeContainer(set, false, 0.938);

    const BeamMoments m = compute(*pc, opalx::moments::allMoments);
    EXPECT_EQ(m.count(), set.size());
    EXPECT_TRUE(m.available().contains(opalx::moments::allMoments.without(Moment::Polarization)));
    expectMatchesReference(m, reference(set), 1e-12);
}

TEST_F(BeamMomentsTest, KineticEnergyUsesReferenceMassWhenSet) {
    const double referenceMassGeV = Physics::m_e;
    auto set                      = correlatedSet(50, 11, referenceMassGeV);
    // The macro-particle mass differs from the species mass; the reference mass must win.
    auto pc = makeContainer(set, false, 1.0e3 * referenceMassGeV);
    const PartData species(1.0, referenceMassGeV * Units::GeV2eV, 1.0);
    pc->setReference(&species);

    const BeamMoments m = compute(*pc, Moment::MeanMomentum | Moment::Spread);
    const Reference ref = reference(set);
    expectRelative(m.meanKineticEnergy(), ref.meanEkin, 1e-12, "meanEkin");
    expectRelative(m.stdKineticEnergy(), ref.stdEkin, 1e-10, "stdEkin");
}

TEST_F(BeamMomentsTest, PerParticleMassWithoutReference) {
    Options::useQMAttributes = true;
    auto set                 = correlatedSet(40, 3, 1.0);
    for (std::size_t i = 0; i < set.size(); ++i) {
        set[i].massGeV = 0.5 + 0.01 * i;
    }
    auto pc         = makeContainer(set);
    const auto mine = localIndices(set.size());
    auto mass       = pc->getMView();
    auto host       = Kokkos::create_mirror_view(mass);
    for (std::size_t i = 0; i < mine.size(); ++i) {
        host(i) = set[mine[i]].massGeV;
    }
    Kokkos::deep_copy(mass, host);

    const BeamMoments m = compute(*pc, Moment::Spread);
    const Reference ref = reference(set);
    expectRelative(m.meanKineticEnergy(), ref.meanEkin, 1e-12, "meanEkin");
    expectRelative(m.stdKineticEnergy(), ref.stdEkin, 1e-10, "stdEkin");
}

TEST_F(BeamMomentsTest, RequestContract) {
    auto pc = makeContainer(correlatedSet(20, 5, 1.0));

    const BeamMoments bounds = compute(*pc, Moment::Bounds);
    EXPECT_EQ(bounds.available(), Request(Moment::Bounds));
    EXPECT_THROW(static_cast<void>(bounds.meanP()), OpalException);
    EXPECT_THROW(static_cast<void>(bounds.rmsR()), OpalException);

    // Spread is centred on both means, so they become available too.
    const BeamMoments spread = compute(*pc, Moment::Spread);
    EXPECT_TRUE(spread.available().contains(Moment::MeanPosition | Moment::MeanMomentum));
    EXPECT_FALSE(spread.available().contains(Moment::Bounds));
    EXPECT_THROW(static_cast<void>(spread.covariance(6, 0)), OpalException);

    // Without spin there is no polarization to report.
    const BeamMoments all = compute(*pc, opalx::moments::allMoments);
    EXPECT_FALSE(all.available().contains(Moment::Polarization));
    EXPECT_THROW(static_cast<void>(all.meanPolarization()), OpalException);
}

TEST_F(BeamMomentsTest, OutputRequestsAreWhatComputeProvides) {
    // The tracker reuses one value for every output of a step only if the value holds the full
    // request, so an output must not ask for a group compute() leaves out (polarization without
    // spin).
    for (const bool spin : {false, true}) {
        auto pc = makeContainer(correlatedSet(20, 31, 1.0), spin);
        for (const bool psDump : {false, true}) {
            for (const bool statDump : {false, true}) {
                const Request request = DataSink::requiredMoments(psDump, statDump, pc->hasSpin());
                EXPECT_TRUE(compute(*pc, request).available().contains(request))
                        << "spin " << spin << ", ps " << psDump << ", stat " << statDump;
            }
        }
    }
}

TEST_F(BeamMomentsTest, SharedSumsDoNotDependOnTheRequest) {
    // Output steps request more groups than ordinary steps; the quantities the step itself
    // uses must come out bitwise the same either way.
    auto pc                   = makeContainer(correlatedSet(300, 13, 1.0));
    const BeamMoments narrow  = compute(*pc, Moment::MeanMomentum);
    const BeamMoments wide    = compute(*pc, opalx::moments::allMoments);
    const BeamMoments bounded = compute(*pc, Moment::Bounds);
    for (int d = 0; d < 3; ++d) {
        EXPECT_EQ(narrow.meanP()[d], wide.meanP()[d]);
        EXPECT_EQ(bounded.minR()[d], wide.minR()[d]);
        EXPECT_EQ(bounded.maxR()[d], wide.maxR()[d]);
    }
    EXPECT_EQ(narrow.meanGamma(), wide.meanGamma());
}

TEST_F(BeamMomentsTest, EmptyContainerGivesZerosAndFiniteDerivedValues) {
    auto pc             = makeContainer({});
    const BeamMoments m = compute(*pc, opalx::moments::allMoments);
    EXPECT_EQ(m.count(), 0u);
    for (int d = 0; d < 3; ++d) {
        EXPECT_EQ(m.minR()[d], 0.0);
        EXPECT_EQ(m.maxR()[d], 0.0);
        EXPECT_EQ(m.meanR()[d], 0.0);
        EXPECT_EQ(m.rmsR()[d], 0.0);
        EXPECT_EQ(m.correlationRP()[d], 0.0);
        EXPECT_EQ(m.normalizedEmittance()[d], 0.0);
        EXPECT_EQ(m.geometricEmittance()[d], 0.0);
    }
    EXPECT_EQ(m.meanGammaZ(), 1.0);
    EXPECT_EQ(m.stdKineticEnergy(), 0.0);
}

TEST_F(BeamMomentsTest, SingleParticleAndZeroExtentHaveNoSpread) {
    std::vector<Particle> one(1);
    one[0].r            = Vec3(1e-3, -2e-3, 0.5);
    one[0].p            = Vec3(0.0, 0.01, 3.0);
    one[0].massGeV      = 1.0;
    const BeamMoments m = compute(*makeContainer(one), opalx::moments::allMoments);
    for (int d = 0; d < 3; ++d) {
        EXPECT_EQ(m.minR()[d], m.maxR()[d]);
        EXPECT_EQ(m.rmsR()[d], 0.0);
        EXPECT_EQ(m.correlationRP()[d], 0.0);
        EXPECT_EQ(m.normalizedEmittance()[d], 0.0);
    }

    // A flat beam (all particles at y = 0 with py = 0) has no y correlation, not NaN.
    auto flat = correlatedSet(30, 17, 1.0);
    for (Particle& q : flat) {
        q.r[1] = 0.0;
        q.p[1] = 0.0;
    }
    const BeamMoments f = compute(*makeContainer(flat), Moment::Spread);
    EXPECT_EQ(f.correlationRP()[1], 0.0);
    EXPECT_EQ(f.normalizedEmittance()[1], 0.0);
}

TEST_F(BeamMomentsTest, UnionEqualsConcatenatedContainer) {
    // Two beams with different species masses (reference masses) and one empty container.
    auto first  = correlatedSet(80, 19, Physics::m_p);
    auto second = correlatedSet(50, 23, Physics::m_e);
    for (Particle& q : second) {
        q.r = Vec3(q.r[0] + 0.5, q.r[1] - 0.3, q.r[2] + 0.2);
        q.p = Vec3(0.5 * q.p[0], 0.5 * q.p[1], 0.5 * q.p[2]);
    }
    auto a     = makeContainer(first);
    auto b     = makeContainer(second);
    auto empty = makeContainer({});
    const PartData proton(1.0, Physics::m_p * Units::GeV2eV, 1.0);
    const PartData electron(-1.0, Physics::m_e * Units::GeV2eV, 1.0);
    a->setReference(&proton);
    b->setReference(&electron);

    const std::vector<const PC_t*> containers = {a.get(), empty.get(), b.get()};
    const BeamMoments m                       = compute(containers, opalx::moments::allMoments);

    std::vector<Particle> all = first;
    all.insert(all.end(), second.begin(), second.end());
    EXPECT_EQ(m.count(), all.size());
    expectMatchesReference(m, reference(all), 1e-12);

    // The empty container must not pull the box to the origin.
    const std::vector<const PC_t*> withEmpty = {empty.get(), b.get()};
    const BeamMoments boxed                  = compute(withEmpty, Moment::Bounds);
    const BeamMoments alone                  = compute(*b, Moment::Bounds);
    for (int d = 0; d < 3; ++d) {
        EXPECT_EQ(boxed.minR()[d], alone.minR()[d]);
        EXPECT_EQ(boxed.maxR()[d], alone.maxR()[d]);
    }
    const std::vector<const PC_t*> onlyEmpty = {empty.get()};
    const BeamMoments none                   = compute(onlyEmpty, Moment::Bounds);
    for (int d = 0; d < 3; ++d) {
        EXPECT_EQ(none.minR()[d], 0.0);
        EXPECT_EQ(none.maxR()[d], 0.0);
    }
}

TEST_F(BeamMomentsTest, PolarizationOnlyWhenEveryContainerTracksSpin) {
    const auto set = correlatedSet(60, 29, 1.0);
    auto spin      = makeContainer(set, true);
    auto plain     = makeContainer(set, false);

    const BeamMoments m = compute(*spin, Moment::Polarization);
    const Reference ref = reference(set);
    for (int d = 0; d < 3; ++d) {
        expectRelative(m.meanPolarization()[d], ref.meanPol[d], 1e-13, "meanPol");
        // One-pass variance as on master: the spread is 2-10 % of the mean here.
        expectRelative(m.rmsPolarization()[d], ref.rmsPol[d], 1e-11, "rmsPol");
    }
    const long double magnitude = std::sqrt(
            ref.meanPol[0] * ref.meanPol[0] + ref.meanPol[1] * ref.meanPol[1]
            + ref.meanPol[2] * ref.meanPol[2]);
    expectRelative(m.meanPolarizationMagnitude(), magnitude, 1e-13, "polMagnitude");

    const std::vector<const PC_t*> mixed = {spin.get(), plain.get()};
    EXPECT_FALSE(compute(mixed, Moment::Polarization).available().contains(Moment::Polarization));
}

namespace {
    double dot3(const Vec3& a, const Vec3& b) { return a(0) * b(0) + a(1) * b(1) + a(2) * b(2); }

    /// Boost a rest-frame momentum u' (beta*gamma) to the frame in which the rest frame moves
    /// with velocity beta.
    Vec3 boostToLab(const Vec3& uRest, const Vec3& beta) {
        const double gamma     = 1.0 / std::sqrt(1.0 - dot3(beta, beta));
        const double gammaRest = std::sqrt(1.0 + dot3(uRest, uRest));
        const double shift = gamma * gamma / (gamma + 1.0) * dot3(beta, uRest) + gamma * gammaRest;
        Vec3 u;
        for (int d = 0; d < 3; ++d) {
            u(d) = uRest(d) + shift * beta(d);
        }
        return u;
    }
}  // namespace

// An anisotropic set with zero total momentum is boosted to the lab; compute() must recover the
// temperatures k_B T_i = m c^2 <u'_i^2 / gamma'> of the rest-frame set.
TEST_F(BeamMomentsTest, RestFrameTemperatureRecoversBoostedAnisotropicSet) {
    std::vector<Vec3> uRest = {
            Vec3(1.0e-2, 0.0, 0.0), Vec3(0.0, 2.0e-2, 0.0), Vec3(0.0, 0.0, 5.0e-3),
            Vec3(5.0e-3, -7.0e-3, 2.0e-3)};
    const std::size_t nHalf = uRest.size();
    for (std::size_t i = 0; i < nHalf; ++i) {
        uRest.push_back(Vec3(-uRest[i](0), -uRest[i](1), -uRest[i](2)));
    }
    const std::size_t n = uRest.size();

    const double massGeV = Physics::m_p;
    const PartData proton(1.0, massGeV * Units::GeV2eV, 1.0);
    Vec3 expected(0.0);
    for (const Vec3& u : uRest) {
        const double gammaRest = std::sqrt(1.0 + dot3(u, u));
        for (int d = 0; d < 3; ++d) {
            expected(d) += massGeV * Units::GeV2eV * u(d) * u(d) / gammaRest / Physics::kB / n;
        }
    }

    const Vec3 zAxis(0.0, 0.0, 1.0);
    const Vec3 oblique(1.0 / std::sqrt(14.0), -2.0 / std::sqrt(14.0), 3.0 / std::sqrt(14.0));
    for (const Vec3& direction : {zAxis, oblique}) {
        for (const double gammaBoost : {1.0, 10.0, 100.0}) {
            const double betaBoost = std::sqrt(1.0 - 1.0 / (gammaBoost * gammaBoost));
            const Vec3 beta(
                    betaBoost * direction(0), betaBoost * direction(1), betaBoost * direction(2));
            std::vector<Particle> set(n);
            for (std::size_t i = 0; i < n; ++i) {
                set[i].r = Vec3(1.0e-3 * i, -2.0e-3 * i * i, 3.0e-3 * (n - i));
                set[i].p = boostToLab(uRest[i], beta);
            }
            // The macro-particle mass differs from the species mass, which must be used.
            auto pc = makeContainer(set, false, 1.0e3 * massGeV);
            pc->setReference(&proton);

            const Vec3 temperature = compute(*pc, Moment::Spread).restFrameTemperature();
            for (int d = 0; d < 3; ++d) {
                EXPECT_NEAR(temperature(d), expected(d), 1e-8 * expected(d))
                        << "gamma = " << gammaBoost << ", direction = " << direction
                        << ", component " << d;
            }
        }
    }
}
