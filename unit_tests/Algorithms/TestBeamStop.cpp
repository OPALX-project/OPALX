#include "gtest/gtest.h"

#include "Ippl.h"

#include "Algorithms/BeamStop.h"

#include <vector>

namespace {
    using Vector = Vector_t<double, 3>;

    class BeamStopTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }
        static void TearDownTestSuite() { ippl::finalize(); }
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
}  // namespace
