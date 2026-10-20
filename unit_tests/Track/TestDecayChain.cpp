#include <gtest/gtest.h>

#include "Track/DecayChain.h"
#include "Utilities/OpalException.h"

#include <string>
#include <vector>

namespace {

    using Names = std::vector<std::string>;
    using Roots = std::vector<std::size_t>;

    TEST(DecayChainTest, ChainFollowsFirstBeam) {
        // PI -> MU -> E: both daughters follow the pion beam.
        EXPECT_EQ(
                resolveDecayRoots({"PI", "MU", "E"}, {"MU", "E", ""}, {true, true, false}),
                (Roots{0, 0, 0}));
    }

    TEST(DecayChainTest, ChainInReverseOrder) {
        EXPECT_EQ(
                resolveDecayRoots({"E", "MU", "PI"}, {"", "E", "MU"}, {false, true, true}),
                (Roots{2, 2, 2}));
    }

    TEST(DecayChainTest, IndependentBeamsFollowThemselves) {
        EXPECT_EQ(resolveDecayRoots({"A", "B"}, {"", ""}, {true, false}), (Roots{0, 1}));
        EXPECT_EQ(
                resolveDecayRoots({"MU", "E", "P"}, {"E", "", ""}, {true, false, false}),
                (Roots{0, 0, 2}));
    }

    TEST(DecayChainTest, RejectsBadSettings) {
        // DAUGHTERBEAM without DECAY
        EXPECT_THROW(resolveDecayRoots({"MU", "E"}, {"E", ""}, {false, false}), OpalException);
        // DAUGHTERBEAM not in the TRACK
        EXPECT_THROW(resolveDecayRoots({"MU"}, {"E"}, {true}), OpalException);
        // a beam names itself
        EXPECT_THROW(resolveDecayRoots({"MU"}, {"MU"}, {true}), OpalException);
        // two beams name the same daughter beam
        EXPECT_THROW(
                resolveDecayRoots({"PI", "MU", "E"}, {"E", "E", ""}, {true, true, false}),
                OpalException);
        // loop
        EXPECT_THROW(resolveDecayRoots({"A", "B"}, {"B", "A"}, {true, true}), OpalException);
    }

}  // namespace
