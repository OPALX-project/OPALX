#include "gtest/gtest.h"

#include "Ippl.h"

#include "../KeywordObjects.h"
#include "Algorithms/BeamStop.h"
#include "Algorithms/StepSizeConfig.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"
#include "Structure/Beam.h"

#include <cmath>
#include <memory>
#include <vector>

extern Inform* gmsg;

namespace {
    using Vector = Vector_t<double, 3>;

    class BeamStopTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
            gmsg = new Inform(nullptr, -1);
        }
        static void TearDownTestSuite() {
            delete gmsg;
            gmsg = nullptr;
            ippl::finalize();
        }

        /// A bunch of @p beams containers whose frames are the lab frame.
        static std::unique_ptr<PartBunch_t> makeBunch(const size_t beams) {
            const auto beam = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
            return std::make_unique<PartBunch_t>(
                    std::vector<double>(beams, 1.0), std::vector<double>(beams, 1.0),
                    std::vector(beams, opalx::run::describeBeamForBunch(*beam)),
                    std::vector<size_t>(beams, 8), 1.0,
                    opalx::spacecharge::CartesianDomainConfig3D{});
        }

        /// Particles on the axis at @p z [m].
        static void placeOnAxis(ParticleContainer_t& pc, const std::vector<double>& z) {
            if (pc.getLocalNum() < z.size()) pc.createParticles(z.size() - pc.getLocalNum());
            auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc.R.getView());
            for (size_t i = 0; i < z.size(); ++i)
                host(i) = Vector({0.0, 0.0, z[i]});
            Kokkos::deep_copy(pc.R.getView(), host);
        }

        /// A tracked bunch of one beam at path length @p s with beta gamma = 1 along z.
        static std::unique_ptr<PartBunch_t> movingBunch(const std::vector<double>& z, double s) {
            auto bunch = makeBunch(1);
            auto pc    = bunch->getParticleContainer(0);
            placeOnAxis(*pc, z);
            bunch->resetPcActive();
            bunch->setdT(dt);
            pc->setRefPartR(Vector({0.0, 0.0, s}));
            pc->setRefPartP(Vector({0.0, 0.0, 1.0}));
            pc->set_sPos(s);
            return bunch;
        }

        static StepSizeConfig schedule(const std::vector<double>& sStops) {
            StepSizeConfig stepSizes;
            for (double sStop : sStops)
                stepSizes.push_back(dt, sStop, 100);
            stepSizes.resetIterator();
            return stepSizes;
        }

        static constexpr double dt = 1.0e-10;
        /// Path length of one step at beta gamma = 1.
        const double drift = dt * Physics::c / std::sqrt(2.0);
    };

    TEST_F(BeamStopTest, ReachedWithinHalfAStep) {
        EXPECT_FALSE(beam_stop::reached(0.90, 1.0, 0.1));
        EXPECT_TRUE(beam_stop::reached(0.96, 1.0, 0.1));
        EXPECT_TRUE(beam_stop::reached(1.20, 1.0, 0.1));
    }

    TEST_F(BeamStopTest, PlaneIsMovedToZStopAlongTheMomentum) {
        BeamStop stop;
        stop.zStop = 2.0;
        // Reference 3 cm before ZSTOP, moving along (3, 0, 4) / 5.
        beam_stop::setPlane(stop, Vector({1.0, 0.0, 1.0}), Vector({0.3, 0.0, 0.4}), 1.97);
        ASSERT_TRUE(stop.planeSet);
        EXPECT_NEAR(stop.planeNormal(0), 0.6, 1e-15);
        EXPECT_NEAR(stop.planeNormal(2), 0.8, 1e-15);
        EXPECT_NEAR(stop.planeOrigin(0), 1.0 + 0.03 * 0.6, 1e-15);
        EXPECT_NEAR(stop.planeOrigin(2), 1.0 + 0.03 * 0.8, 1e-15);
    }

    TEST_F(BeamStopTest, LastParticleDecidesWhetherThePlaneIsPassed) {
        Kokkos::View<Vector*> positions("positions", 3);
        auto host = Kokkos::create_mirror_view(positions);
        host(0)   = Vector({0.0, 0.0, 2.1});
        host(1)   = Vector({0.0, 0.0, 1.95});  // still behind the plane at z = 2
        host(2)   = Vector({0.0, 0.0, 2.3});
        Kokkos::deep_copy(positions, host);

        const Vector origin({0.0, 0.0, 2.0}), normal({0.0, 0.0, 1.0});
        EXPECT_NEAR(beam_stop::minDistanceToPlane(positions, 3, origin, normal), -0.05, 1e-15);
        EXPECT_NEAR(beam_stop::minDistanceToPlane(positions, 1, origin, normal), 0.1, 1e-15);
        EXPECT_GT(beam_stop::minDistanceToPlane(positions, 0, origin, normal), 0.0);
    }

    TEST_F(BeamStopTest, RunEndsWhenEveryZStopBeamIsDoneAndNoneRunsToMaxSteps) {
        // Without any ZSTOP the run ends by MAXSTEPS.
        EXPECT_FALSE(beam_stop::runFinished({BeamStop{}, BeamStop{}}));

        BeamStop primary;
        primary.zStop = 10.0;
        BeamStop other;
        EXPECT_FALSE(beam_stop::runFinished({primary, other}));
        primary.done = true;
        EXPECT_TRUE(beam_stop::runFinished({primary, other}));

        // A second beam with its own ZSTOP keeps the run going until it is done too.
        BeamStop second;
        second.zStop = 20.0;
        EXPECT_FALSE(beam_stop::runFinished({primary, second}));
        second.done = true;
        EXPECT_TRUE(beam_stop::runFinished({primary, second}));

        // UNTILMAXSTEPS keeps the run going regardless.
        other.untilMaxSteps = true;
        EXPECT_FALSE(beam_stop::runFinished({primary, second, other}));
    }

    TEST_F(BeamStopTest, AnnounceKeepsOneStopPerContainer) {
        BeamStop primary;
        primary.zStop = 1.0;
        BeamStop untilMaxSteps;
        untilMaxSteps.untilMaxSteps = true;
        BeamStops stops({primary, untilMaxSteps, untilMaxSteps});

        stops.announce(2);
        ASSERT_EQ(stops.stops().size(), 2u);
        EXPECT_EQ(stops.stops()[0].zStop, 1.0);
        EXPECT_TRUE(stops.stops()[1].untilMaxSteps);

        // A container without a stop of its own is tracked while the run lasts.
        stops.announce(3);
        ASSERT_EQ(stops.stops().size(), 3u);
        EXPECT_FALSE(stops.stops()[2].zStop);
        EXPECT_FALSE(stops.stops()[2].untilMaxSteps);
    }

    TEST_F(BeamStopTest, DoneBeamsStayStopped) {
        auto bunch = makeBunch(2);
        placeOnAxis(*bunch->getParticleContainer(0), {0.0});
        placeOnAxis(*bunch->getParticleContainer(1), {0.0});
        bunch->resetPcActive();
        BeamStop done;
        done.zStop = 1.0;
        done.done  = true;
        BeamStops stops({BeamStop{}, done});

        stops.keepDoneBeamsStopped(*bunch);
        EXPECT_TRUE(bunch->isPcActive(0));
        EXPECT_FALSE(bunch->pcAtSStop(0));
        EXPECT_FALSE(bunch->isPcActive(1));
        EXPECT_TRUE(bunch->pcAtSStop(1));
    }

    TEST_F(BeamStopTest, PrimaryReferenceAtSegmentStopSwitchesSegment) {
        auto bunch         = movingBunch({0.0}, 1.0 - 0.6 * drift);
        StepSizeConfig dts = schedule({1.0, 2.0});
        BeamStops stops({BeamStop{}});
        EXPECT_FALSE(stops.applyStopRules(*bunch, dts, false).switchSegment);

        bunch->getParticleContainer(0)->set_sPos(1.0 - 0.4 * drift);
        const auto rules = stops.applyStopRules(*bunch, dts, false);
        EXPECT_TRUE(rules.switchSegment);
        EXPECT_FALSE(rules.runFinished);

        // The last segment has no next one.
        ++dts;
        bunch->getParticleContainer(0)->set_sPos(2.0);
        EXPECT_FALSE(stops.applyStopRules(*bunch, dts, false).switchSegment);
    }

    TEST_F(BeamStopTest, ZStopPlaneIsSetByTheReferenceAndPassedByTheLastParticle) {
        auto bunch               = movingBunch({1.05, 0.99, 1.02}, 1.0 - 0.6 * drift);
        auto pc                  = bunch->getParticleContainer(0);
        const StepSizeConfig dts = schedule({10.0});
        BeamStop stop;
        stop.zStop = 1.0;
        BeamStops stops({stop});

        EXPECT_FALSE(stops.applyStopRules(*bunch, dts, false).runFinished);
        EXPECT_FALSE(stops.stops()[0].planeSet);

        // Within half a step of ZSTOP the plane is set at ZSTOP; a particle is still behind it.
        pc->setRefPartR(Vector({0.0, 0.0, 1.0 - 0.4 * drift}));
        pc->set_sPos(1.0 - 0.4 * drift);
        EXPECT_FALSE(stops.applyStopRules(*bunch, dts, false).runFinished);
        ASSERT_TRUE(stops.stops()[0].planeSet);
        EXPECT_NEAR(stops.stops()[0].planeOrigin(2), 1.0, 1e-15);
        EXPECT_FALSE(stops.stops()[0].done);
        EXPECT_TRUE(bunch->isPcActive(0));

        placeOnAxis(*pc, {1.05, 1.001, 1.02});
        EXPECT_TRUE(stops.applyStopRules(*bunch, dts, false).runFinished);
        EXPECT_TRUE(stops.stops()[0].done);
        EXPECT_TRUE(bunch->pcAtSStop(0));
    }

    TEST_F(BeamStopTest, ReferenceHitEndsTheRun) {
        auto bunch               = movingBunch({0.0}, 0.0);
        const StepSizeConfig dts = schedule({10.0});
        BeamStops stops({BeamStop{}});
        EXPECT_FALSE(stops.applyStopRules(*bunch, dts, false).runFinished);
        EXPECT_TRUE(stops.applyStopRules(*bunch, dts, true).runFinished);
    }
}  // namespace
