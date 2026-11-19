#include <gtest/gtest.h>

#include "SpaceCharge/SpaceChargeSolver.h"
#include "Utilities/OpalException.h"

#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

namespace opalx::spacecharge {
    namespace {

        class RecordingAlgorithm final : public SpaceChargeAlgorithm {
        public:
            explicit RecordingAlgorithm(SpaceChargeSolveResult result) : result_m(result) {}

            SpaceChargeSolveResult solve(const SpaceChargeSolveContext& context) override {
                ++calls;
                activity.assign(context.trackingActive().begin(), context.trackingActive().end());
                return result_m;
            }

            SolveFrame solveFrame() const override { return SolveFrame::Reference; }

            SpaceChargeSolveResult result_m;
            std::size_t calls = 0;
            std::vector<std::uint8_t> activity;
        };

        class ThrowingAlgorithm final : public SpaceChargeAlgorithm {
        public:
            SpaceChargeSolveResult solve(const SpaceChargeSolveContext&) override {
                throw std::runtime_error("deliberate algorithm failure");
            }

            SolveFrame solveFrame() const override { return SolveFrame::Reference; }
        };

        TEST(SpaceChargeSolverTest, ValidatesContainerCountAndAggregatesResults) {
            auto algorithm = std::make_unique<RecordingAlgorithm>(SpaceChargeSolveResult{
                    .backendSolves = 3, .redistributions = 1, .reportedBins = 7});
            RecordingAlgorithm* observer = algorithm.get();
            SpaceChargeSolver solver(std::move(algorithm), 2);

            const std::vector<std::uint8_t> activity{1, 0};
            SpaceChargeSolveContext context(activity, {});
            solver.solve(context);
            solver.solve(context);

            EXPECT_EQ(observer->calls, 2u);
            EXPECT_EQ(observer->activity, activity);
            EXPECT_EQ(solver.backendSolveCount(), 6u);
            EXPECT_EQ(solver.redistributionCount(), 2u);
            EXPECT_EQ(solver.reportedBinCount(), 7);

            const std::vector<std::uint8_t> wrongActivity{1};
            SpaceChargeSolveContext wrongContext(wrongActivity, {});
            EXPECT_THROW(solver.solve(wrongContext), OpalException);
            EXPECT_EQ(observer->calls, 2u);
        }

        TEST(SpaceChargeSolverTest, RejectsInvalidConstruction) {
            EXPECT_THROW(SpaceChargeSolver(nullptr, 1), OpalException);
            EXPECT_THROW(
                    SpaceChargeSolver(
                            std::make_unique<RecordingAlgorithm>(SpaceChargeSolveResult{}), 0),
                    OpalException);
        }

        TEST(SpaceChargeSolverTest, PropagatesAlgorithmFailures) {
            SpaceChargeSolver solver(std::make_unique<ThrowingAlgorithm>(), 1);
            const std::vector<std::uint8_t> activity{1};
            SpaceChargeSolveContext context(activity, {});
            EXPECT_THROW(solver.solve(context), std::runtime_error);
            EXPECT_EQ(solver.backendSolveCount(), 0u);
        }

        TEST(SolveFrameTest, RequestsBuildTheirFrames) {
            using Vec = ippl::Vector<double, 3>;
            const Vec meanR(1e-3, -2e-3, 5e-4), meanP(0.02, -0.01, 3.0), fallback(0, 0, 1);
            const double norm =
                    std::sqrt(meanP[0] * meanP[0] + meanP[1] * meanP[1] + meanP[2] * meanP[2]);
            const auto expectNear = [](const Vec& a, const Vec& b) {
                for (unsigned d = 0; d < 3; ++d) EXPECT_NEAR(a[d], b[d], 1e-14) << d;
            };

            EXPECT_TRUE(makeSolveFrame(SolveFrame::Reference, meanR, meanP, fallback, 7.0)
                                .isReference());

            // Centroid: the centroid is the origin and the mean momentum points along +z.
            const auto centroid =
                    makeSolveFrame(SolveFrame::Centroid, meanR, meanP, fallback, 7.0);
            expectNear(centroid.fromReference().transformTo(meanR), Vec(0.0));
            expectNear(
                    centroid.fromReference().rotateTo(meanP), Vec(0, 0, norm));

            // PathLength: the reference particle sits at z = path length on the aligned axis.
            const auto path = makeSolveFrame(SolveFrame::PathLength, meanR, meanP, fallback, 7.0);
            expectNear(path.fromReference().transformTo(Vec(0.0)), Vec(0, 0, 7.0));
            expectNear(path.fromReference().rotateTo(meanP), Vec(0, 0, norm));

            // A vanishing mean momentum takes the fallback direction.
            const auto still =
                    makeSolveFrame(SolveFrame::Centroid, meanR, Vec(0.0), Vec(1, 0, 0), 0.0);
            expectNear(still.fromReference().rotateTo(Vec(1, 0, 0)), Vec(0, 0, 1));
        }

    }  // namespace
}  // namespace opalx::spacecharge
