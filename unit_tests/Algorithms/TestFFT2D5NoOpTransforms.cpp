#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include "BeamlineCore/MonitorRep.h"
#include "BeamlineCore/SBendRep.h"
#include "Ippl.h"
#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/FFT2D5/MomentUpdatePolicy.h"
#include "SpaceCharge/FFT2D5/NoOpTransformPolicy.h"

namespace {
    using namespace opalx::spacecharge::fft2d5_detail;
    using PC  = ParticleContainer<double, 3>;
    using Vec = Vector_t<double, 3>;

    class FieldDrift : public DriftRep {
    public:
        FieldDrift() : DriftRep("custom-field-drift") {}
        void apply(const std::shared_ptr<PC>& pc) override {
            const auto r = pc->R.getView(), b = pc->B.getView();
            Kokkos::parallel_for(
                    "NoOpTest::field", pc->getLocalNum(), KOKKOS_LAMBDA(size_t i) {
                        b(i) += Vec(.1 * r(i)[0], .2 * r(i)[1], .3 * r(i)[2]);
                    });
        }
    };

    class FFT2D5NoOpTransforms : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }
        static void TearDownTestSuite() { ippl::finalize(); }

        std::shared_ptr<PC> particles(const std::vector<Vec>& localPositions) {
            const ippl::Vector<int, 3> nr(8);
            const Vec origin(-4.), spacing(1.);
            ippl::NDIndex<3> domain;
            for (unsigned d = 0; d < 3; ++d)
                domain[d] = ippl::Index(nr[d]);
            Mesh_t<3> mesh(domain, spacing, origin);
            FieldLayout_t<3> layout(MPI_COMM_WORLD, domain, {true, true, true}, true);
            auto pc = std::make_shared<PC>(mesh, layout);
            pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
            pc->createParticles(localPositions.size());
            auto r = pc->R.getHostMirror();
            for (size_t i = 0; i < localPositions.size(); ++i)
                r(i) = localPositions[i];
            Kokkos::deep_copy(pc->R.getView(), r);
            Kokkos::deep_copy(pc->P.getView(), Vec(.03, -.02, .4));
            Kokkos::deep_copy(pc->E.getView(), Vec(.2, .3, -.1));
            Kokkos::deep_copy(pc->B.getView(), Vec(-.4, .1, .2));
            return pc;
        }

        std::vector<Vec> positions() {
            // Rank zero owns all particles, exercising a real empty MPI rank.
            return ippl::Comm->rank() == 0 ? std::vector<Vec>{Vec(.1, .02, .3), Vec(-.2, .01, .7)}
                                           : std::vector<Vec>{};
        }

        CoordinateSystemTrafo misalignedFrame() {
            const double a = .19;
            return CoordinateSystemTrafo(
                           Vec(.4, -.7, .2),
                           Quaternion(
                                   std::cos(a), std::sin(a) / std::sqrt(3.),
                                   std::sin(a) / std::sqrt(3.), std::sin(a) / std::sqrt(3.)))
                   * CoordinateSystemTrafo(
                           Vec(-.03, .02, -.01), Quaternion(std::cos(.07), 0., std::sin(.07), 0.));
        }

        void sameState(const PC& a, const PC& b) {
            for (const auto& pair :
                 {std::make_pair(a.R.getView(), b.R.getView()),
                  std::make_pair(a.P.getView(), b.P.getView()),
                  std::make_pair(a.E.getView(), b.E.getView()),
                  std::make_pair(a.B.getView(), b.B.getView())}) {
                auto ah = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pair.first);
                auto bh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pair.second);
                for (size_t i = 0; i < a.getLocalNum(); ++i)
                    for (unsigned d = 0; d < 3; ++d)
                        EXPECT_EQ(ah(i)[d], bh(i)[d]);
            }
        }
    };

    TEST_F(FFT2D5NoOpTransforms, FilteringIsConservativeAndFallbackKeepsEveryCallback) {
        DriftRep drift("drift");
        MarkerRep marker("marker");
        MonitorRep monitor("monitor");
        SBendRep bend("bend");
        FieldDrift custom;
        for (auto* element : std::vector<ElementBase*>{&drift, &marker, &monitor, &bend, &custom}) {
            for (auto operation :
                 {FrameOperation::All, FrameOperation::Fields, FrameOperation::Monitors,
                  FrameOperation::Apertures}) {
                EXPECT_TRUE(callbackNeedsFrame(*element, operation, true, false));
            }
            EXPECT_TRUE(callbackNeedsFrame(*element, FrameOperation::All, true, true));
            EXPECT_EQ(
                    callbackNeedsFrame(*element, FrameOperation::Monitors, true, true),
                    element == &monitor);
            element->setFlagDeleteOnTransverseExit(false);
            EXPECT_FALSE(callbackNeedsFrame(*element, FrameOperation::Apertures, true, true));
            element->setFlagDeleteOnTransverseExit(true);
            EXPECT_TRUE(callbackNeedsFrame(*element, FrameOperation::Apertures, true, true));
        }
        EXPECT_FALSE(callbackNeedsFrame(drift, FrameOperation::Fields, true, true));
        EXPECT_FALSE(callbackNeedsFrame(marker, FrameOperation::Fields, true, true));
        EXPECT_TRUE(callbackNeedsFrame(custom, FrameOperation::Fields, true, true));
        EXPECT_TRUE(callbackNeedsFrame(bend, FrameOperation::Fields, true, true));
        EXPECT_FALSE(callbackNeedsFrame(monitor, FrameOperation::Fields, true, true));
        EXPECT_TRUE(callbackNeedsFrame(monitor, FrameOperation::Fields, false, true));
    }

    TEST_F(FFT2D5NoOpTransforms, SkippedCallbacksLeaveArraysAndMomentValidityUntouched) {
        auto pc = particles(positions()), reference = particles(positions());
        DriftRep drift("empty");
        drift.setFlagDeleteOnTransverseExit(false);
        unsigned callbacks = 0;
        pc->markMomentsClean();
        for (auto operation :
             {FrameOperation::Fields, FrameOperation::Monitors, FrameOperation::Apertures}) {
            EXPECT_FALSE(inElementFrame(drift, *pc, misalignedFrame(), operation, true, true, [&] {
                ++callbacks;
            }));
        }
        EXPECT_EQ(callbacks, 0u);
        EXPECT_FALSE(pc->isMomentsDirty());
        sameState(*pc, *reference);
    }

    TEST_F(FFT2D5NoOpTransforms, UnknownFieldCallbackMatchesLegacyMisalignedFrameSequence) {
        auto pc = particles(positions()), reference = particles(positions());
        FieldDrift element;
        const auto frame = misalignedFrame();
        reference->transformBunch(frame);
        element.apply(reference);
        reference->transformBunch(frame.inverted());
        unsigned calls = 0;
        EXPECT_TRUE(inElementFrame(element, *pc, frame, FrameOperation::Fields, true, true, [&] {
            ++calls;
            element.apply(pc);
        }));
        EXPECT_EQ(calls, 1u);
        EXPECT_TRUE(pc->isMomentsDirty());
        sameState(*pc, *reference);
    }

    TEST_F(FFT2D5NoOpTransforms, MonitorCallbackRetainedAndSeesElementLocalCoordinates) {
        auto pc = particles(positions());
        MonitorRep monitor("observation");
        const auto frame = misalignedFrame();
        // A host mirror view may alias live CPU storage; this snapshot must own its data.
        auto initial = Kokkos::create_mirror(pc->R.getView());
        Kokkos::deep_copy(initial, pc->R.getView());
        unsigned calls = 0;
        EXPECT_TRUE(inElementFrame(monitor, *pc, frame, FrameOperation::Monitors, true, true, [&] {
            ++calls;
            auto local = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
            for (size_t i = 0; i < pc->getLocalNum(); ++i) {
                const auto expected = frame.transformTo(initial(i));
                for (unsigned d = 0; d < 3; ++d)
                    EXPECT_NEAR(local(i)[d], expected[d], 1.e-15);
            }
        }));
        EXPECT_EQ(calls, 1u);
    }

    TEST_F(FFT2D5NoOpTransforms, EnabledApertureRetainsEntranceExitAndMisalignmentBehavior) {
        DriftRep aperture("aperture");
        aperture.getGeometry().setElementLength(1.);
        aperture.setAperture(ApertureType::RECTANGULAR, {1., .2});
        aperture.setFlagDeleteOnTransverseExit(true);
        for (const auto& frame : {CoordinateSystemTrafo(), misalignedFrame()}) {
            std::vector<Vec> local =
                    ippl::Comm->rank() == 0
                            ? std::vector<
                                      Vec>{Vec(0., 0., .5), Vec(1.1, 0., .5), Vec(1.1, 0., 0.), Vec(1.1, 0., 1.), Vec(0., 0., .999)}
                            : std::vector<Vec>{};
            std::vector<Vec> tracker;
            for (const auto& r : local)
                tracker.push_back(frame.transformFrom(r));
            auto reference = particles(tracker), pc = particles(tracker);
            reference->transformBunch(frame);
            const auto expectedCount = aperture.markOutsideAperture(reference);
            reference->transformBunch(frame.inverted());
            size_t actualCount = 0;
            EXPECT_TRUE(inElementFrame(
                    aperture, *pc, frame, FrameOperation::Apertures, true, true, [&] {
                        actualCount = aperture.markOutsideAperture(pc);
                    }));
            EXPECT_EQ(actualCount, expectedCount);
            auto a = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), pc->InvalidMask.getView());
            auto b = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), reference->InvalidMask.getView());
            for (size_t i = 0; i < pc->getLocalNum(); ++i)
                EXPECT_EQ(a(i), b(i));
            sameState(*pc, *reference);
        }
    }

    TEST_F(FFT2D5NoOpTransforms, EligibilityRejectsActiveEmissionAndOtherSolverModes) {
        MomentUpdateEligibility eligible{true, true, false, false, false, false, true, 1};
        EXPECT_TRUE(eligible.mayDefer());
        for (int changed = 0; changed < 7; ++changed) {
            auto state = eligible;
            if (changed == 0) state.fft2d5 = false;
            if (changed == 1) state.emissionActive = true;
            if (changed == 2) state.midpoint = false;
            if (changed == 3) state.boundaryControlled = true;
            if (changed == 4) state.spinEnabled = true;
            if (changed == 5) state.globalProcesses = true;
            if (changed == 6) state.containerCount = 2;
            EXPECT_FALSE(state.mayDefer());
        }
    }
}  // namespace
