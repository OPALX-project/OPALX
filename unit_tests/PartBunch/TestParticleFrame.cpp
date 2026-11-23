/**
 * \file TestParticleFrame.cpp
 * \brief Unit tests for ParticleFrame and ParticleContainer::enterFrame / adoptParticles.
 *
 * Covers:
 * - A request for the frame the data is already in does not touch it.
 * - A change from or to the reference frame applies the stored transform itself.
 * - A direct change between two element frames equals the path through the reference frame.
 * - `clearFields` zeroes E and B instead of rotating them.
 * - A FrameScope returns the container to its reference frame when an exception unwinds it.
 * - `adoptParticles` moves only the appended particles.
 * - The reference update refuses data outside the reference frame.
 */

#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <memory>
#include <stdexcept>
#include <vector>

#include "Ippl.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Utilities/OpalException.h"

namespace {

    using PC_t   = ParticleContainer<double, 3>;
    using Vec    = Vector_t<double, 3>;
    using Kind   = ParticleFrame::Kind;
    using Points = std::vector<Vec>;

    class ParticleFrameTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }

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
            auto pc = std::make_shared<PC_t>(mesh, fl);
            pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
            return pc;
        }

        /// Append particles with R = r[i], P = r[i] scaled, E = 2 r[i], B = 3 r[i].
        static void append(PC_t& pc, const Points& r) {
            const size_t first = pc.getLocalNum();
            pc.createParticles(r.size());
            auto R = pc.R.getHostMirror(), P = pc.P.getHostMirror();
            auto E = pc.E.getHostMirror(), B = pc.B.getHostMirror();
            Kokkos::deep_copy(R, pc.R.getView());
            Kokkos::deep_copy(P, pc.P.getView());
            Kokkos::deep_copy(E, pc.E.getView());
            Kokkos::deep_copy(B, pc.B.getView());
            for (size_t i = 0; i < r.size(); ++i) {
                R(first + i) = r[i];
                P(first + i) = Vec(0.1 * r[i][0], -0.2 * r[i][1], 5.0 + r[i][2]);
                E(first + i) = 2.0 * r[i];
                B(first + i) = 3.0 * r[i];
            }
            Kokkos::deep_copy(pc.R.getView(), R);
            Kokkos::deep_copy(pc.P.getView(), P);
            Kokkos::deep_copy(pc.E.getView(), E);
            Kokkos::deep_copy(pc.B.getView(), B);
        }

        template <class Attribute>
        static Points host(const Attribute& a, size_t n) {
            auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a.getView());
            Points out(n);
            for (size_t i = 0; i < n; ++i)
                out[i] = h(i);
            return out;
        }

        static ParticleFrame elementFrame(Vec origin, Vec axis, double tilt) {
            const CoordinateSystemTrafo onAxis(origin, getQuaternion(axis, Vec(0, 0, 1)));
            const CoordinateSystemTrafo tilted(
                    Vec(0.0),
                    getQuaternion(Vec(std::sin(tilt), 0.0, std::cos(tilt)), Vec(0, 0, 1)));
            return ParticleFrame(Kind::Element, tilted * onAxis);
        }

        const Points r_ = {Vec(1e-3, -2e-3, 3e-3), Vec(-4e-3, 5e-3, -6e-3), Vec(0.0, 0.0, 0.0)};
    };

    TEST_F(ParticleFrameTest, RequestForCurrentFrameDoesNotTouchData) {
        auto pc = makeContainer();
        append(*pc, r_);
        const auto frame = elementFrame(Vec(0.2, -0.1, 1.5), Vec(0.1, 0.05, 1.0), 0.3);
        pc->enterFrame(frame);
        const auto before = host(pc->R, r_.size());

        pc->enterFrame(elementFrame(Vec(0.2, -0.1, 1.5), Vec(0.1, 0.05, 1.0), 0.3));
        EXPECT_EQ(pc->frame().kind(), Kind::Element);
        const auto after = host(pc->R, r_.size());
        for (size_t i = 0; i < r_.size(); ++i)
            for (unsigned d = 0; d < 3; ++d)
                EXPECT_EQ(after[i][d], before[i][d]);

        pc->enterReferenceFrame();
        EXPECT_TRUE(pc->frame().isReference());
    }

    /// Same origin and quaternion, bit for bit.
    static void expectSameTransform(
            const CoordinateSystemTrafo& a, const CoordinateSystemTrafo& b) {
        for (unsigned d = 0; d < 3; ++d)
            EXPECT_EQ(a.getOrigin()[d], b.getOrigin()[d]);
        for (unsigned d = 0; d < 4; ++d)
            EXPECT_EQ(a.getRotation()(d), b.getRotation()(d));
    }

    TEST_F(ParticleFrameTest, ChangeFromOrToReferenceUsesStoredTransform) {
        const auto frame = elementFrame(Vec(0.2, -0.1, 1.5), Vec(0.1, 0.05, 1.0), 0.3);
        expectSameTransform(ParticleFrame().to(frame), frame.fromReference());
        expectSameTransform(frame.to(ParticleFrame()), frame.toReference());

        // The kernel may contract to FMA on a device; compare to the host within rounding.
        auto pc = makeContainer();
        append(*pc, r_);
        pc->enterFrame(frame);
        const auto R = host(pc->R, r_.size()), P = host(pc->P, r_.size());
        for (size_t i = 0; i < r_.size(); ++i) {
            const Vec r = frame.fromReference().transformTo(r_[i]);
            const Vec p = frame.fromReference().rotateTo(
                    Vec(0.1 * r_[i][0], -0.2 * r_[i][1], 5.0 + r_[i][2]));
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_NEAR(R[i][d], r[d], 1e-15);
                EXPECT_NEAR(P[i][d], p[d], 1e-14);
            }
        }
    }

    TEST_F(ParticleFrameTest, DirectChangeEqualsPathThroughReference) {
        const auto first  = elementFrame(Vec(0.2, -0.1, 1.5), Vec(0.1, 0.05, 1.0), 0.3);
        const auto second = elementFrame(Vec(-0.3, 0.4, -2.0), Vec(-0.2, 0.1, 1.0), -0.7);

        auto direct = makeContainer(), viaReference = makeContainer();
        append(*direct, r_);
        append(*viaReference, r_);
        direct->enterFrame(first);
        direct->enterFrame(second);
        viaReference->enterFrame(first);
        viaReference->enterReferenceFrame();
        viaReference->enterFrame(second);

        const auto a = host(direct->R, r_.size()), b = host(viaReference->R, r_.size());
        const auto ea = host(direct->E, r_.size()), eb = host(viaReference->E, r_.size());
        for (size_t i = 0; i < r_.size(); ++i)
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_NEAR(a[i][d], b[i][d], 1e-14);
                EXPECT_NEAR(ea[i][d], eb[i][d], 1e-14);
            }
        EXPECT_TRUE(direct->frame() == second);
    }

    TEST_F(ParticleFrameTest, ClearFieldsZeroesInsteadOfRotating) {
        auto pc = makeContainer();
        append(*pc, r_);
        const auto expectZeroFields = [&] {
            const auto E = host(pc->E, r_.size()), B = host(pc->B, r_.size());
            for (size_t i = 0; i < r_.size(); ++i)
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_EQ(E[i][d], 0.0);
                    EXPECT_EQ(B[i][d], 0.0);
                }
        };
        pc->enterFrame(elementFrame(Vec(0.2, -0.1, 1.5), Vec(0.1, 0.05, 1.0), 0.3), true);
        expectZeroFields();

        pc->E = Vec(1.0, 2.0, 3.0);
        pc->B = Vec(-1.0, 0.5, 2.0);
        pc->enterFrame(pc->frame(), true);
        expectZeroFields();
    }

    TEST_F(ParticleFrameTest, FrameScopeReturnsToReferenceFrameOnException) {
        auto pc = makeContainer();
        append(*pc, r_);
        try {
            const PC_t::FrameScope scope(
                    *pc, elementFrame(Vec(0.2, -0.1, 1.5), Vec(0.1, 0.05, 1.0), 0.3));
            EXPECT_FALSE(pc->frame().isReference());
            throw std::runtime_error("unwind");
        } catch (const std::runtime_error&) {
        }
        EXPECT_TRUE(pc->frame().isReference());
        const auto R = host(pc->R, r_.size());
        for (size_t i = 0; i < r_.size(); ++i)
            for (unsigned d = 0; d < 3; ++d)
                EXPECT_NEAR(R[i][d], r_[i][d], 1e-14);
    }

    TEST_F(ParticleFrameTest, AdoptMovesOnlyAppendedParticles) {
        auto pc = makeContainer();
        append(*pc, r_);
        const auto old = host(pc->R, r_.size());

        const auto source    = elementFrame(Vec(0.0, 0.0, -0.5), Vec(0.0, 0.1, 1.0), 0.0);
        const Points created = {Vec(1e-3, 1e-3, 0.0), Vec(-1e-3, 2e-3, 1e-4)};
        append(*pc, created);
        pc->adoptParticles(source, r_.size());

        EXPECT_TRUE(pc->frame().isReference());
        const auto R = host(pc->R, r_.size() + created.size());
        for (size_t i = 0; i < r_.size(); ++i)
            for (unsigned d = 0; d < 3; ++d)
                EXPECT_EQ(R[i][d], old[i][d]);
        for (size_t i = 0; i < created.size(); ++i) {
            const Vec expected = source.toReference().transformTo(created[i]);
            for (unsigned d = 0; d < 3; ++d)
                EXPECT_NEAR(R[r_.size() + i][d], expected[d], 1e-15);
        }
    }

    TEST_F(ParticleFrameTest, ReferenceUpdateRequiresReferenceFrame) {
        auto pc = makeContainer();
        append(*pc, r_);
        pc->getRefPartP() = Vec(0.0, 0.0, 5.0);
        pc->enterFrame(elementFrame(Vec(0.2, -0.1, 1.5), Vec(0.1, 0.05, 1.0), 0.3));
        EXPECT_THROW(pc->updateRefToLabCSTrafo(1e-12), OpalException);
        EXPECT_THROW(pc->setToLabTrafo(CoordinateSystemTrafo()), OpalException);
        pc->enterReferenceFrame();
        EXPECT_NO_THROW(pc->updateRefToLabCSTrafo(1e-12));
    }

}  // namespace
