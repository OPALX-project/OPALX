/**
 * \file TestCollimator.cpp
 * \brief Unit tests for the COLLIMATOR element (CollimatorRep).
 *
 * ---------------------------------------------------------------------------
 * Coverage
 * ---------------------------------------------------------------------------
 *
 * 1. Element identity
 *    - getType(), getFieldExtent() (the body span), one time step per element,
 *      clone() copies name, geometry and aperture
 *
 * 2. Absorption
 *    - the inherited ElementBase::markOutsideAperture marks particles inside
 *      the body window [0, L) but outside an elliptical / rectangular aperture
 *    - setFlagDeleteOnTransverseExit(false) turns the collimator off
 *    - marked particles are removed by deleteInvalidParticles()
 *
 * A collimator adds no aperture logic of its own; these tests exist to pin
 * down that a Collimator really inherits the generic behaviour (the shape and
 * z-window semantics themselves are covered by TestMarkOutsideAperture).
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

#include "BeamlineCore/CollimatorRep.h"
#include "PartBunch/ParticleContainer.hpp"

namespace {

    using PC_t = ParticleContainer<double, 3>;

    class CollimatorTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }

        /// Collimator body of length L with the given aperture.
        static CollimatorRep makeCollimator(
                double length, ApertureType type, double xLimit, double yLimit) {
            CollimatorRep coll("test_collimator");
            coll.getGeometry().setElementLength(length);
            coll.setAperture(type, {xLimit, yLimit});
            return coll;
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
    // Element identity
    // ================================================================

    TEST_F(CollimatorTest, GetType) {
        CollimatorRep coll = makeCollimator(0.1, ApertureType::ELLIPTICAL, 0.01, 0.01);
        EXPECT_EQ(coll.getType(), ElementType::COLLIMATOR);
        EXPECT_EQ(coll.getTypeString(), "Collimator");
    }

    TEST_F(CollimatorTest, RequiredNumberOfTimeStepsIsOne) {
        // The ElementBase default of 10 makes OrbitThreader::checkElementLengths
        // throw for short elements, which a collimator typically is.
        CollimatorRep coll = makeCollimator(0.1, ApertureType::ELLIPTICAL, 0.01, 0.01);
        EXPECT_EQ(coll.getRequiredNumberOfTimeSteps(), 1);
    }

    TEST_F(CollimatorTest, FieldExtentIsTheBodySpan) {
        CollimatorRep coll = makeCollimator(0.1, ApertureType::ELLIPTICAL, 0.01, 0.01);

        double zBegin = -1.0, zEnd = -1.0;
        coll.getFieldExtent(zBegin, zEnd);
        EXPECT_DOUBLE_EQ(zBegin, 0.0);
        EXPECT_DOUBLE_EQ(zEnd, 0.1);
    }

    TEST_F(CollimatorTest, CloneCopiesGeometryAndAperture) {
        CollimatorRep coll = makeCollimator(0.1, ApertureType::RECTANGULAR, 0.01, 0.02);

        std::unique_ptr<ElementBase> copy(coll.clone());
        ASSERT_NE(copy, nullptr);
        EXPECT_EQ(copy->getName(), "test_collimator");
        EXPECT_EQ(copy->getType(), ElementType::COLLIMATOR);
        EXPECT_DOUBLE_EQ(copy->getGeometry().getElementLength(), 0.1);

        const auto aperture = copy->getAperture();
        EXPECT_EQ(aperture.first, ApertureType::RECTANGULAR);
        ASSERT_EQ(aperture.second.size(), 2u);
        EXPECT_DOUBLE_EQ(aperture.second[0], 0.01);
        EXPECT_DOUBLE_EQ(aperture.second[1], 0.02);
    }

    // ================================================================
    // Absorption
    // ================================================================

    TEST_F(CollimatorTest, MarksParticlesOutsideAnEllipticalAperture) {
        CollimatorRep coll = makeCollimator(0.1, ApertureType::ELLIPTICAL, 0.01, 0.005);
        auto pc            = makeContainer();

        createParticlesAt(
                pc, {
                            {0.000, 0.000, 0.05},   // 0: on axis -> kept
                            {0.020, 0.000, 0.05},   // 1: outside in x -> marked
                            {0.000, -0.01, 0.05},   // 2: outside in y -> marked
                            {0.009, 0.000, 0.05},   // 3: inside -> kept
                            {0.009, 0.004, 0.05},   // 4: outside (x/a)^2+(y/b)^2 > 1 -> marked
                            {0.020, 0.000, -0.05},  // 5: upstream of the body -> kept
                            {0.020, 0.000, 0.15},   // 6: downstream of the body -> kept
                    });

        EXPECT_EQ(coll.markOutsideAperture(pc), 3u);
        EXPECT_EQ(
                invalidMaskOnHost(pc),
                std::vector<bool>({false, true, true, false, true, false, false}));
    }

    TEST_F(CollimatorTest, MarksParticlesOutsideARectangularAperture) {
        CollimatorRep coll = makeCollimator(0.1, ApertureType::RECTANGULAR, 0.01, 0.005);
        auto pc            = makeContainer();

        createParticlesAt(
                pc, {
                            {0.009, 0.004, 0.05},  // 0: inside the rectangle -> kept
                            {0.011, 0.000, 0.05},  // 1: outside in x -> marked
                            {0.000, 0.006, 0.05},  // 2: outside in y -> marked
                    });

        // The corner particle 0 would be outside an ellipse with the same
        // half-widths; the shapes must not be confused.
        EXPECT_EQ(coll.markOutsideAperture(pc), 2u);
        EXPECT_EQ(invalidMaskOnHost(pc), std::vector<bool>({false, true, true}));
    }

    TEST_F(CollimatorTest, RespectsDeleteOnTransverseExitFlag) {
        CollimatorRep coll = makeCollimator(0.1, ApertureType::ELLIPTICAL, 0.01, 0.01);
        coll.setFlagDeleteOnTransverseExit(false);
        auto pc = makeContainer();

        createParticlesAt(pc, {{0.02, 0.0, 0.05}});

        EXPECT_EQ(coll.markOutsideAperture(pc), 0u);
        EXPECT_EQ(invalidMaskOnHost(pc), std::vector<bool>({false}));
    }

    TEST_F(CollimatorTest, AbsorbedParticlesAreDeleted) {
        CollimatorRep coll = makeCollimator(0.1, ApertureType::ELLIPTICAL, 0.01, 0.01);
        auto pc            = makeContainer();

        createParticlesAt(
                pc, {
                            {0.000, 0.00, 0.05},  // kept
                            {0.020, 0.00, 0.05},  // absorbed
                            {0.000, 0.02, 0.05},  // absorbed
                            {0.005, 0.00, 0.05},  // kept
                    });

        EXPECT_EQ(coll.markOutsideAperture(pc), 2u);
        EXPECT_EQ(pc->deleteInvalidParticles(), 2u);
        EXPECT_EQ(pc->getTotalNum(), 2u);
    }

}  // namespace
