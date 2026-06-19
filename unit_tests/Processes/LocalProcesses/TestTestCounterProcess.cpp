//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <memory>
#include <vector>

#include "Ippl.h"

#include "BeamlineCore/DriftRep.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Processes/LocalProcesses/TestCounterProcess.h"

namespace {

    using PC_t = ParticleContainer<double, 3>;

    class TestCounterProcessTest : public ::testing::Test {
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
            for (unsigned i = 0; i < 3; ++i) {
                domain[i] = ippl::Index(nr[i]);
            }

            Mesh_t<3> mesh(domain, hr, origin);
            FieldLayout_t<3> fl(MPI_COMM_WORLD, domain, decomp, true);
            auto pc = std::make_shared<PC_t>(mesh, fl);
            pc->setBunchStateHandler(std::make_shared<BunchStateHandler>());
            return pc;
        }

        /// Create particles on the beam axis with the given longitudinal positions.
        void createParticlesAtZ(std::shared_ptr<PC_t>& pc, const std::vector<double>& zPositions) {
            if (zPositions.empty()) {
                return;
            }

            pc->createParticles(zPositions.size());
            auto R_host = pc->R.getHostMirror();
            auto P_host = pc->P.getHostMirror();

            for (size_t i = 0; i < zPositions.size(); ++i) {
                R_host(i)[0] = 0.0;
                R_host(i)[1] = 0.0;
                R_host(i)[2] = zPositions[i];
                P_host(i)[0] = 0.0;
                P_host(i)[1] = 0.0;
                P_host(i)[2] = 1.0;
            }

            Kokkos::deep_copy(pc->R.getView(), R_host);
            Kokkos::deep_copy(pc->P.getView(), P_host);
            Kokkos::fence();
        }

        /// Drift of length 1 m as the host element.
        DriftRep makeElement() {
            DriftRep drift("D1");
            drift.getGeometry().setElementLength(1.0);
            return drift;
        }

        LocalProcessContext makeContext(long long step) {
            return LocalProcessContext{1.0e-12, step, 0, CoordinateSystemTrafo()};
        }
    };

}  // namespace

TEST_F(TestCounterProcessTest, CountsParticlesInsideElementBody) {
    auto pc       = makeContainer();
    DriftRep elem = makeElement();
    // 3 inside [0, 1), 3 outside.
    createParticlesAtZ(pc, {0.0, 0.5, 0.999, -0.1, 1.0, 2.5});

    TestCounterProcess process("TP1");
    const size_t marked = process.apply(elem, *pc, makeContext(7));

    EXPECT_EQ(marked, 0u);
    EXPECT_EQ(process.getInvocationCount(), 1ll);
    EXPECT_EQ(process.getLastInsideCount(), 3u);
    EXPECT_EQ(process.getTotalInsideCount(), 3u);
    EXPECT_EQ(process.getLastGlobalTrackStep(), 7ll);
}

TEST_F(TestCounterProcessTest, RepeatedApplyAccumulates) {
    auto pc       = makeContainer();
    DriftRep elem = makeElement();
    createParticlesAtZ(pc, {0.25, 0.75, 5.0});

    TestCounterProcess process("TP1");
    process.apply(elem, *pc, makeContext(1));
    process.apply(elem, *pc, makeContext(2));

    EXPECT_EQ(process.getInvocationCount(), 2ll);
    EXPECT_EQ(process.getLastInsideCount(), 2u);
    EXPECT_EQ(process.getTotalInsideCount(), 4u);
    EXPECT_EQ(process.getLastGlobalTrackStep(), 2ll);
}

TEST_F(TestCounterProcessTest, EmptyContainerCountsZeroButInvokes) {
    auto pc       = makeContainer();
    DriftRep elem = makeElement();

    TestCounterProcess process("TP1");
    const size_t marked = process.apply(elem, *pc, makeContext(3));

    EXPECT_EQ(marked, 0u);
    EXPECT_EQ(process.getInvocationCount(), 1ll);
    EXPECT_EQ(process.getLastInsideCount(), 0u);
    EXPECT_EQ(process.getTotalInsideCount(), 0u);
}

TEST_F(TestCounterProcessTest, DoesNotModifyParticleState) {
    auto pc       = makeContainer();
    DriftRep elem = makeElement();
    const std::vector<double> zPositions = {0.1, 0.6, 1.4};
    createParticlesAtZ(pc, zPositions);

    TestCounterProcess process("TP1");
    process.apply(elem, *pc, makeContext(1));

    auto R_host = pc->R.getHostMirror();
    auto P_host = pc->P.getHostMirror();
    Kokkos::deep_copy(R_host, pc->R.getView());
    Kokkos::deep_copy(P_host, pc->P.getView());

    for (size_t i = 0; i < zPositions.size(); ++i) {
        EXPECT_DOUBLE_EQ(R_host(i)[0], 0.0);
        EXPECT_DOUBLE_EQ(R_host(i)[1], 0.0);
        EXPECT_DOUBLE_EQ(R_host(i)[2], zPositions[i]);
        EXPECT_DOUBLE_EQ(P_host(i)[0], 0.0);
        EXPECT_DOUBLE_EQ(P_host(i)[1], 0.0);
        EXPECT_DOUBLE_EQ(P_host(i)[2], 1.0);
    }
    EXPECT_EQ(pc->getLocalNum(), zPositions.size());
}
