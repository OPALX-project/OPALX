/**
 * \file TestMarkOutsideAperture.cpp
 * \brief Unit tests for ElementBase::markOutsideAperture.
 *
 * ---------------------------------------------------------------------------
 * Coverage
 * ---------------------------------------------------------------------------
 *
 * 1. Base implementation (DriftRep on a minimal 8^3 mesh)
 *    - marks exactly the particles inside the z-window and outside the
 *      aperture in InvalidMask, and returns the newly-marked count
 *    - honors setFlagDeleteOnTransverseExit(false)
 *    - never clears marks set by other producers and does not recount them
 *    - marked particles are removed by deleteInvalidParticles()
 *    - empty container / nullptr return 0
 *
 * 2. z-window source (MonitorRep)
 *    - the window comes from getFieldExtent(), which is centered on the
 *      monitor plane ([-halfLength, +halfLength]), not [0, L)
 *
 * 3. Arc-coordinate override (SBendRep)
 *    - the z-window and the aperture are measured in arc coordinates
 *      relative to the design orbit, not in the straight entrance frame
 *
 * Particles are created directly in element-local coordinates (z measured
 * from the entrance face), matching the frame ParallelTracker's element loop
 * provides via transformBunch().
 */

#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <memory>
#include <vector>

#include "Ippl.h"

#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/MonitorRep.h"
#include "BeamlineCore/SBendRep.h"
#include "PartBunch/ParticleContainer.hpp"

namespace {

    using PC_t = ParticleContainer<double, 3>;

    class MarkOutsideApertureTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }

        /// Drift body of length L with a circular aperture of radius r.
        static DriftRep makeDrift(double length, double radius) {
            DriftRep drift("test_drift");
            drift.getGeometry().setElementLength(length);
            drift.setAperture(ApertureType::ELLIPTICAL, {radius, radius});
            return drift;
        }

        /// Build a minimal ParticleContainer on an 8^3 periodic mesh over [-4, 4]^3.
        std::shared_ptr<PC_t> makeContainer() {
            ippl::Vector<int, 3> nr        = 8;
            ippl::Vector<double, 3> rmin   = -4.0;
            ippl::Vector<double, 3> rmax   = 4.0;
            ippl::Vector<double, 3> origin = rmin;
            ippl::Vector<double, 3> hr     = (rmax - rmin) / ippl::Vector<double, 3>(nr);
            std::array<bool, 3> decomp     = {true, true, true};

            ippl::NDIndex<3> domain;
            for (unsigned i = 0; i < 3; i++) {
                domain[i] = ippl::Index(nr[i]);
            }

            Mesh_t<3> mesh(domain, hr, origin);
            FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, true);

            std::shared_ptr<PC_t> pc = std::make_shared<PC_t>(mesh, fl);
            pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
            return pc;
        }

        /// Create local particles at the given element-local positions.
        void createParticlesAt(
                std::shared_ptr<PC_t>& pc, const std::vector<std::array<double, 3>>& positions) {
            const size_t n = positions.size();
            if (n == 0) return;

            pc->createParticles(n);

            auto R_host = pc->R.getHostMirror();
            for (size_t i = 0; i < n; ++i) {
                R_host(i)[0] = positions[i][0];
                R_host(i)[1] = positions[i][1];
                R_host(i)[2] = positions[i][2];
            }
            Kokkos::deep_copy(pc->R.getView(), R_host);
            Kokkos::fence();
        }

        /// Host copy of the InvalidMask for assertions.
        static std::vector<bool> invalidMaskOnHost(const std::shared_ptr<PC_t>& pc) {
            auto mask = pc->InvalidMask.getHostMirror();
            Kokkos::deep_copy(mask, pc->InvalidMask.getView());
            std::vector<bool> result(pc->getLocalNum());
            for (size_t i = 0; i < result.size(); ++i) {
                result[i] = mask(i);
            }
            return result;
        }
    };

    // ================================================================
    // Base implementation (DriftRep)
    // ================================================================

    TEST_F(MarkOutsideApertureTest, MarksOnlyParticlesInsideZAndOutsideAperture) {
        DriftRep drift = makeDrift(0.1, 0.01);
        auto pc        = makeContainer();

        createParticlesAt(
                pc, {
                            {0.000, 0.000, 0.05},   // 0: in z, in aperture -> kept
                            {0.020, 0.000, 0.05},   // 1: in z, out of aperture -> marked
                            {0.000, -0.02, 0.05},   // 2: in z, out of aperture -> marked
                            {0.020, 0.000, -0.05},  // 3: upstream of body -> kept
                            {0.020, 0.000, 0.15},   // 4: downstream of body -> kept
                            {0.009, 0.000, 0.00},   // 5: entrance face, in aperture -> kept
                    });

        EXPECT_EQ(drift.markOutsideAperture(pc), 2u);

        const std::vector<bool> mask     = invalidMaskOnHost(pc);
        const std::vector<bool> expected = {false, true, true, false, false, false};
        EXPECT_EQ(mask, expected);
    }

    TEST_F(MarkOutsideApertureTest, RespectsDeleteOnTransverseExitFlag) {
        DriftRep drift = makeDrift(0.1, 0.01);
        drift.setFlagDeleteOnTransverseExit(false);
        auto pc = makeContainer();

        createParticlesAt(pc, {{0.02, 0.0, 0.05}});

        EXPECT_EQ(drift.markOutsideAperture(pc), 0u);
        EXPECT_EQ(invalidMaskOnHost(pc), std::vector<bool>({false}));
    }

    TEST_F(MarkOutsideApertureTest, NeverClearsExistingMarksAndDoesNotRecountThem) {
        DriftRep drift = makeDrift(0.1, 0.01);
        auto pc        = makeContainer();

        // The first particle passes the aperture but is pre-marked by another
        // producer (e.g. a decay process); the second is out of the aperture
        // and pre-marked -> already counted by that producer, not recounted.
        createParticlesAt(pc, {{0.0, 0.0, 0.05}, {0.02, 0.0, 0.05}});

        auto mask_host = pc->InvalidMask.getHostMirror();
        mask_host(0)   = true;
        mask_host(1)   = true;
        Kokkos::deep_copy(pc->InvalidMask.getView(), mask_host);
        Kokkos::fence();

        EXPECT_EQ(drift.markOutsideAperture(pc), 0u);
        EXPECT_EQ(invalidMaskOnHost(pc), std::vector<bool>({true, true}));
    }

    TEST_F(MarkOutsideApertureTest, MarkedParticlesAreDeletedByDeleteInvalidParticles) {
        DriftRep drift = makeDrift(0.1, 0.01);
        auto pc        = makeContainer();

        createParticlesAt(
                pc, {
                            {0.000, 0.00, 0.05},  // kept
                            {0.020, 0.00, 0.05},  // marked
                            {0.000, 0.02, 0.05},  // marked
                            {0.005, 0.00, 0.05},  // kept
                    });

        EXPECT_EQ(drift.markOutsideAperture(pc), 2u);
        EXPECT_EQ(pc->deleteInvalidParticles(), 2u);
        EXPECT_EQ(pc->getTotalNum(), 2u);

        // survivors are the in-aperture particles
        auto R_host = pc->R.getHostMirror();
        Kokkos::deep_copy(R_host, pc->R.getView());
        for (size_t i = 0; i < pc->getLocalNum(); ++i) {
            const double x = R_host(i)[0];
            const double y = R_host(i)[1];
            EXPECT_LT(x * x + y * y, 0.01 * 0.01);
        }
    }

    TEST_F(MarkOutsideApertureTest, EmptyContainerAndNullptrReturnZero) {
        DriftRep drift = makeDrift(0.1, 0.01);
        auto pc        = makeContainer();

        EXPECT_EQ(drift.markOutsideAperture(pc), 0u);
        EXPECT_EQ(drift.markOutsideAperture(nullptr), 0u);
    }

    // ================================================================
    // z-window source (MonitorRep)
    // ================================================================

    TEST_F(MarkOutsideApertureTest, UsesBodyWindowNotFieldExtent) {
        // The marking window is the geometric body [0, L), NOT getFieldExtent().
        // A monitor is the sharpest case: its field extent is centered on the
        // recording plane ([-halfLength, +halfLength]) and so reaches upstream
        // of the body, while the body starts at z = 0.
        const double length = 0.01;  // what OpalMonitor::update installs
        MonitorRep monitor("test_monitor");
        monitor.getGeometry().setElementLength(length);
        monitor.setAperture(ApertureType::ELLIPTICAL, {0.01, 0.01});

        double fieldBegin = 0.0, fieldEnd = 0.0;
        monitor.getFieldExtent(fieldBegin, fieldEnd);
        ASSERT_LT(fieldBegin, 0.0) << "field extent must reach upstream for this test to bite";

        auto pc = makeContainer();
        createParticlesAt(
                pc, {
                            {0.020, 0.0, 0.5 * length},       // in body, out of aperture -> marked
                            {0.020, 0.0, 0.5 * fieldBegin},   // in field extent but upstream of
                                                              // the body -> kept
                            {0.020, 0.0, 2.0 * length},       // downstream of the body -> kept
                            {0.005, 0.0, 0.5 * length},       // in body, in aperture -> kept
                    });

        EXPECT_EQ(monitor.markOutsideAperture(pc), 1u);
        EXPECT_EQ(invalidMaskOnHost(pc), std::vector<bool>({true, false, false, false}));
    }

    // ================================================================
    // Arc-coordinate override (SBendRep)
    // ================================================================

    TEST_F(MarkOutsideApertureTest, SBendMeasuresApertureInArcCoordinates) {
        // 1 m body with curvature 1/m (1 rad bend), 5 cm circular aperture.
        const double curvature = 1.0;
        const double length    = 1.0;
        const double radius    = 1.0 / curvature;

        SBendRep sbend("test_sbend");
        sbend.getGeometry().setElementLength(length);
        sbend.getGeometry().setCurvature(curvature);
        sbend.setAperture(ApertureType::ELLIPTICAL, {0.05, 0.05});

        // Entry-frame position of a particle at bend angle phi with radial
        // offset d from the design orbit (circle of radius 1/curvature
        // centered at x = -radius).
        const auto entryPos = [&](double phi, double d) -> std::array<double, 3> {
            return {-radius + (radius + d) * std::cos(phi), 0.0, (radius + d) * std::sin(phi)};
        };

        const std::array<double, 3> onOrbit  = entryPos(0.5, 0.0);   // arc (0, 0, 0.5)
        const std::array<double, 3> offOrbit = entryPos(0.5, 0.1);   // arc (0.1, 0, 0.5)

        // Verify the frame conversion assumptions of this test.
        const Vector_t<double, 3> arcOn = GeometryHelper::toBendArcCoords(
                Vector_t<double, 3>(onOrbit[0], onOrbit[1], onOrbit[2]), curvature, length);
        ASSERT_NEAR(arcOn(0), 0.0, 1.0e-12);
        ASSERT_NEAR(arcOn(2), 0.5, 1.0e-12);

        // In the straight entrance frame the on-orbit particle is > 0.05 off
        // axis and the off-orbit one is not — a straight-frame check would
        // invert the expected marking.
        ASSERT_GT(std::abs(onOrbit[0]), 0.05);
        ASSERT_LT(std::abs(offOrbit[0]), 0.05);

        auto pc = makeContainer();
        createParticlesAt(pc, {onOrbit, offOrbit});

        EXPECT_EQ(sbend.markOutsideAperture(pc), 1u);
        EXPECT_EQ(invalidMaskOnHost(pc), std::vector<bool>({false, true}));
    }

    TEST_F(MarkOutsideApertureTest, SBendMarksOnlyInsideTheBodyNotTheFringe) {
        // With a non-zero gap the Enge field extends one fringe half width past
        // each face, but the aperture belongs to the body: a particle still
        // upstream of the entrance face must not be scraped.
        const double curvature = 1.0;
        const double length    = 1.0;
        const double radius    = 1.0 / curvature;

        SBendRep sbend("test_sbend_fringe");
        sbend.getGeometry().setElementLength(length);
        sbend.getGeometry().setCurvature(curvature);
        sbend.setFullGap(0.02);
        sbend.setAperture(ApertureType::ELLIPTICAL, {0.05, 0.05});

        double fieldBegin = 0.0, fieldEnd = 0.0;
        sbend.getFieldExtent(fieldBegin, fieldEnd);
        ASSERT_LT(fieldBegin, 0.0) << "gap must produce a fringe reaching upstream";
        ASSERT_GT(fieldEnd, length) << "gap must produce a fringe reaching downstream";

        // Straight entrance frame: for z < 0 the arc transform is the identity,
        // so these sit at arc s < 0 -- inside the field extent, outside the body.
        const double zUp = 0.5 * fieldBegin;

        // Downstream of the exit face, on the straight exit tangent.
        const auto exitPos = [&](double extra, double d) -> std::array<double, 3> {
            const double phi = curvature * length;
            const double cx  = -radius + (radius + d) * std::cos(phi);
            const double cz  = (radius + d) * std::sin(phi);
            return {cx - extra * std::sin(phi), 0.0, cz + extra * std::cos(phi)};
        };
        const std::array<double, 3> downstream = exitPos(0.5 * (fieldEnd - length), 0.1);

        auto pc = makeContainer();
        createParticlesAt(
                pc, {
                            {0.10, 0.0, zUp},  // upstream fringe, outside aperture -> kept
                            downstream,        // downstream fringe, outside aperture -> kept
                    });

        EXPECT_EQ(sbend.markOutsideAperture(pc), 0u);
        EXPECT_EQ(invalidMaskOnHost(pc), std::vector<bool>({false, false}));
    }

}  // namespace
