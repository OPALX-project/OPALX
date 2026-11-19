/**
 * \file TestElementFrame.cpp
 * \brief Unit tests for the frame elements request (ElementBase::fieldFrame) and for aperture
 * checks that read local positions through a RigidTransform.
 *
 * Covers:
 * - Elements that declare ElementFrame::None leave R, P, E and B unchanged in apply(pc), so the
 *   tracker may skip the call and its transform.
 * - Elements with fields keep the default ElementFrame::Local.
 * - markOutsideAperture with a transform marks the same particles as on a bunch moved into the
 *   element frame, and does not move the particles.
 */

#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <memory>
#include <vector>

#include "Ippl.h"

#include "BeamlineCore/BoxRep.h"
#include "BeamlineCore/CollimatorRep.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/LaserRep.h"
#include "BeamlineCore/MarkerRep.h"
#include "BeamlineCore/MonitorRep.h"
#include "BeamlineCore/MultipoleRep.h"
#include "BeamlineCore/ProbeRep.h"
#include "PartBunch/ParticleContainer.hpp"

namespace {

    using PC_t = ParticleContainer<double, 3>;
    using Vec  = Vector_t<double, 3>;

    class ElementFrameTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }

        std::shared_ptr<PC_t> makeContainer(const std::vector<Vec>& r) {
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

            pc->createParticles(r.size());
            auto R = pc->R.getHostMirror(), P = pc->P.getHostMirror();
            auto E = pc->E.getHostMirror(), B = pc->B.getHostMirror();
            for (size_t i = 0; i < r.size(); ++i) {
                R(i) = r[i];
                P(i) = Vec(0.01, -0.02, 2.0);
                E(i) = Vec(1.0, 2.0, 3.0);
                B(i) = Vec(-0.1, 0.2, -0.3);
            }
            Kokkos::deep_copy(pc->R.getView(), R);
            Kokkos::deep_copy(pc->P.getView(), P);
            Kokkos::deep_copy(pc->E.getView(), E);
            Kokkos::deep_copy(pc->B.getView(), B);
            return pc;
        }

        template <class Attribute>
        static std::vector<Vec> host(const Attribute& a, size_t n) {
            auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a.getView());
            std::vector<Vec> out(n);
            for (size_t i = 0; i < n; ++i)
                out[i] = h(i);
            return out;
        }

        static std::vector<int> invalid(const PC_t& pc) {
            auto h = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), pc.InvalidMask.getView());
            std::vector<int> out(pc.getLocalNum());
            for (size_t i = 0; i < out.size(); ++i)
                out[i] = h(i) ? 1 : 0;
            return out;
        }

        /// apply(pc) of a frameless element must leave every particle attribute unchanged.
        static void expectNoOp(ElementBase& element, const std::shared_ptr<PC_t>& pc) {
            SCOPED_TRACE(element.getName());
            const size_t n = pc->getLocalNum();
            const auto R0  = host(pc->R, n);
            const auto P0  = host(pc->P, n);
            const auto E0  = host(pc->E, n);
            const auto B0  = host(pc->B, n);
            element.getGeometry().setElementLength(1.0);
            EXPECT_EQ(element.fieldFrame(), ElementFrame::None);
            element.apply(pc);
            const auto R = host(pc->R, n), P = host(pc->P, n);
            const auto E = host(pc->E, n), B = host(pc->B, n);
            for (size_t i = 0; i < n; ++i)
                for (unsigned d = 0; d < 3; ++d) {
                    EXPECT_EQ(R[i][d], R0[i][d]);
                    EXPECT_EQ(P[i][d], P0[i][d]);
                    EXPECT_EQ(E[i][d], E0[i][d]);
                    EXPECT_EQ(B[i][d], B0[i][d]);
                }
        }

        /// Positions inside a 1 m body, inside and outside a 1 cm aperture.
        const std::vector<Vec> r_ = {
                Vec(0.0, 0.0, 0.5), Vec(0.02, 0.0, 0.5), Vec(0.0, -0.03, 0.2),
                Vec(0.004, 0.003, 0.9), Vec(0.5, 0.5, 1.5)};
    };

    TEST_F(ElementFrameTest, FramelessElementsDoNotTouchParticles) {
        auto pc = makeContainer(r_);
        DriftRep drift("drift");
        MarkerRep marker("marker");
        CollimatorRep collimator("collimator");
        BoxRep box("box");
        LaserRep laser("laser");
        ProbeRep probe("probe");
        for (ElementBase* element :
             std::vector<ElementBase*>{&drift, &marker, &collimator, &box, &laser, &probe})
            expectNoOp(*element, pc);
    }

    TEST_F(ElementFrameTest, FieldElementsKeepLocalFrame) {
        MultipoleRep multipole("multipole");
        MonitorRep monitor("monitor");
        EXPECT_EQ(multipole.fieldFrame(), ElementFrame::Local);
        EXPECT_EQ(monitor.fieldFrame(), ElementFrame::Local);
    }

    TEST_F(ElementFrameTest, ApertureThroughTransformMatchesMovedBunch) {
        // Element frame: shifted by more than the aperture and tilted relative to the particles'
        // frame, so a transform applied the wrong way round marks different particles.
        const CoordinateSystemTrafo toLocal(
                Vec(0.008, -0.006, -0.1), getQuaternion(Vec(0.01, -0.02, 1.0), Vec(0, 0, 1)));
        DriftRep drift("drift");
        drift.getGeometry().setElementLength(1.0);
        drift.setAperture(ApertureType::ELLIPTICAL, {0.01, 0.01});

        // r_ holds element-local positions; the particles are stored in the other frame.
        std::vector<Vec> stored;
        for (const auto& local : r_)
            stored.push_back(toLocal.transformFrom(local));
        const std::vector<int> expected = {0, 1, 1, 0, 0};

        auto lazy = makeContainer(stored), moved = makeContainer(stored);
        const auto R0 = host(lazy->R, stored.size());
        EXPECT_EQ(drift.markOutsideAperture(lazy, toLocal.rigid()), 2u);
        EXPECT_EQ(invalid(*lazy), expected);

        moved->enterFrame(ParticleFrame(ParticleFrame::Kind::Element, toLocal, &drift));
        EXPECT_EQ(drift.markOutsideAperture(moved, RigidTransform{}), 2u);
        EXPECT_EQ(invalid(*moved), expected);

        const auto R = host(lazy->R, stored.size());
        for (size_t i = 0; i < stored.size(); ++i)
            for (unsigned d = 0; d < 3; ++d)
                EXPECT_EQ(R[i][d], R0[i][d]);
    }

}  // namespace
