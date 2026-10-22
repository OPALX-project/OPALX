#include <gtest/gtest.h>

#include "SpaceCharge/SpaceChargeSolver.h"
#include "Utilities/OpalException.h"

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

            SpaceChargeSolveResult result_m;
            std::size_t calls = 0;
            std::vector<std::uint8_t> activity;
        };

        class ThrowingAlgorithm final : public SpaceChargeAlgorithm {
        public:
            SpaceChargeSolveResult solve(const SpaceChargeSolveContext&) override {
                throw std::runtime_error("deliberate algorithm failure");
            }
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

    }  // namespace
}  // namespace opalx::spacecharge
