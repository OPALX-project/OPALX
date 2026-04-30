//
// Unit tests for LossDataSink collective-save edge cases.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//

#include "Structure/LossDataSink.h"

#include "gtest/gtest.h"

class TestLossDataSink : public testing::Test {
public:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() { ippl::finalize(); }
};

TEST_F(TestLossDataSink, SaveWithoutParticlesReturnsWithoutCollectiveAbort) {
    /**
     * Temporal monitors may trigger save() while the reference-particle path is
     * being traced, before any particle sample has been stored locally. This
     * regression verifies that the empty-save fast path participates in the
     * collective emptiness check without relying on an in-place MPI_Reduce.
     */
    LossDataSink sink("unit_test_empty_loss_sink", false, CollectionType::TEMPORAL);
    EXPECT_NO_THROW(sink.save());
}
