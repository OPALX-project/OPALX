/**
 * \file TestBox.cpp
 * \brief Unit tests for the BOX element (BoxRep).
 *
 * ---------------------------------------------------------------------------
 * Coverage
 * ---------------------------------------------------------------------------
 *
 * 1. Element identity
 *    - getType(), getFieldExtent() (the body span), one time step per element,
 *      clone() copies name, geometry and block size
 *
 * 2. Absorption
 *    - Box::markOutsideAperture marks particles inside the block
 *      |x| < halfWidth, |y| < halfHeight, 0 <= z < L, and nothing else
 *    - the aperture stays wide open, so the element is still registered when
 *      the reference path passes beside the block
 *    - setFlagDeleteOnTransverseExit(false) turns the box off
 *    - marked particles are removed by deleteInvalidParticles()
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

#include "BeamlineCore/BoxRep.h"
#include "PartBunch/ParticleContainer.hpp"

namespace {

    using PC_t = ParticleContainer<double, 3>;

    class BoxTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }

        /// Block of length L with the given half sizes.
        static BoxRep makeBox(double length, double halfWidth, double halfHeight) {
            BoxRep box("test_box");
            box.getGeometry().setElementLength(length);
            box.setHalfSize(halfWidth, halfHeight);
            return box;
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

    TEST_F(BoxTest, GetType) {
        BoxRep box = makeBox(0.1, 0.01, 0.005);
        EXPECT_EQ(box.getType(), ElementType::BOX);
        EXPECT_EQ(box.getTypeString(), "Box");
    }

    TEST_F(BoxTest, RequiredNumberOfTimeStepsIsOne) {
        BoxRep box = makeBox(0.1, 0.01, 0.005);
        EXPECT_EQ(box.getRequiredNumberOfTimeSteps(), 1);
    }

    TEST_F(BoxTest, FieldExtentIsTheBodySpan) {
        BoxRep box = makeBox(0.1, 0.01, 0.005);

        double zBegin = -1.0, zEnd = -1.0;
        box.getFieldExtent(zBegin, zEnd);
        EXPECT_DOUBLE_EQ(zBegin, 0.0);
        EXPECT_DOUBLE_EQ(zEnd, 0.1);
    }

    TEST_F(BoxTest, CloneCopiesGeometryAndSize) {
        BoxRep box = makeBox(0.1, 0.01, 0.02);

        std::unique_ptr<ElementBase> copy(box.clone());
        ASSERT_NE(copy, nullptr);
        EXPECT_EQ(copy->getName(), "test_box");
        EXPECT_EQ(copy->getType(), ElementType::BOX);
        EXPECT_DOUBLE_EQ(copy->getGeometry().getElementLength(), 0.1);

        const Box* copyBox = dynamic_cast<const Box*>(copy.get());
        ASSERT_NE(copyBox, nullptr);
        EXPECT_DOUBLE_EQ(copyBox->getHalfWidth(), 0.01);
        EXPECT_DOUBLE_EQ(copyBox->getHalfHeight(), 0.02);
    }

    TEST_F(BoxTest, ApertureStaysWideOpen) {
        // A block beside the reference path must still be registered by
        // ElementBase::isInside, which tests the reference particle against
        // the aperture.
        BoxRep box = makeBox(0.1, 0.01, 0.005);
        EXPECT_TRUE(box.isInside(Vector_t<double, 3>({-0.5, 0.3, 0.05})));
    }

    // ================================================================
    // Absorption
    // ================================================================

    TEST_F(BoxTest, MarksParticlesInsideTheBlock) {
        BoxRep box = makeBox(0.1, 0.01, 0.005);
        auto pc    = makeContainer();

        createParticlesAt(
                pc, {
                            {0.000, 0.000, 0.05},    // 0: centre -> marked
                            {0.009, -0.004, 0.05},   // 1: inside near a corner -> marked
                            {0.011, 0.000, 0.05},    // 2: outside in x -> kept
                            {0.000, -0.006, 0.05},   // 3: outside in y -> kept
                            {0.010, 0.000, 0.05},    // 4: on the x face -> kept
                            {0.000, 0.000, -0.001},  // 5: upstream of the body -> kept
                            {0.000, 0.000, 0.0},     // 6: on the entrance face -> marked
                            {0.000, 0.000, 0.1},     // 7: on the exit face -> kept
                    });

        EXPECT_EQ(box.markOutsideAperture(pc), 3u);
        EXPECT_EQ(
                invalidMaskOnHost(pc),
                std::vector<bool>({true, true, false, false, false, false, true, false}));
    }

    TEST_F(BoxTest, RespectsDeleteOnTransverseExitFlag) {
        BoxRep box = makeBox(0.1, 0.01, 0.005);
        box.setFlagDeleteOnTransverseExit(false);
        auto pc = makeContainer();

        createParticlesAt(pc, {{0.0, 0.0, 0.05}});

        EXPECT_EQ(box.markOutsideAperture(pc), 0u);
        EXPECT_EQ(invalidMaskOnHost(pc), std::vector<bool>({false}));
    }

    TEST_F(BoxTest, AbsorbedParticlesAreDeleted) {
        BoxRep box = makeBox(0.1, 0.01, 0.005);
        auto pc    = makeContainer();

        createParticlesAt(
                pc, {
                            {0.000, 0.000, 0.05},  // absorbed
                            {0.020, 0.000, 0.05},  // kept
                            {0.000, 0.020, 0.05},  // kept
                            {0.005, 0.001, 0.05},  // absorbed
                    });

        EXPECT_EQ(box.markOutsideAperture(pc), 2u);
        EXPECT_EQ(pc->deleteInvalidParticles(), 2u);
        EXPECT_EQ(pc->getTotalNum(), 2u);
    }

}  // namespace
