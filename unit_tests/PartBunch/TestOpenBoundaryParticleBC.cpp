/**
 * @file TestOpenBoundaryParticleBC.cpp
 * @brief Verify that open field solver sets particle boundary condition to NO.
 */

#include <gtest/gtest.h>

#include "PartBunch/PartBunch.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Structure/Beam.h"
#include "Structure/FieldSolverCmd.h"
#include "Structure/DataSink.h"
#include "Utilities/Options.h"
#include "Ippl.h"

// Reuse the fixture from TestBinnedFieldSolver for convenience.
#include "TestBinnedFieldSolver.cpp"

// The fixture class BinnedFieldSolverSmokeTest is defined in TestBinnedFieldSolver.cpp.

TEST_F(BinnedFieldSolverSmokeTest, OpenSolver_ParticleBCIsNo) {
    // Build a bunch with OPEN solver (non-P3M).
    rebuildOpenBunchWithGreensFunction("STANDARD");
    // The container should use spatial layout, not P3M layout.
    EXPECT_FALSE(pc->hasP3MLayout());
    // The particle BC should be NO (open).
    auto bcArray = pc->getPL().getParticleBC();
    for (const auto& bc : bcArray) {
        EXPECT_EQ(bc, ippl::BC::NO);
    }
}
