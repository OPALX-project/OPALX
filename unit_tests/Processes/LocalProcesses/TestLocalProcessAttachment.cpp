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
#include "BeamlineCore/DriftRep.h"
#include "Processes/LocalProcesses/LocalProcess.h"
#include "gtest/gtest.h"

#include <memory>
#include <string>
#include <vector>

namespace {

    /// Minimal concrete process: does nothing, used only to test attachment plumbing.
    class NoOpProcess : public LocalProcess {
    public:
        explicit NoOpProcess(const std::string& name) : LocalProcess(name) {}

        std::size_t apply(Component&, ParticleContainer<double, 3>&,
                          const LocalProcessContext&) override {
            return 0;
        }
    };

}  // namespace

TEST(TestLocalProcessAttachment, DefaultElementHasNoProcesses) {
    DriftRep drift("D1");

    EXPECT_FALSE(drift.hasLocalProcesses());
    EXPECT_TRUE(drift.getLocalProcesses().empty());
}

TEST(TestLocalProcessAttachment, SetLocalProcessesStoresInOrder) {
    DriftRep drift("D1");

    auto p1 = std::make_shared<NoOpProcess>("P1");
    auto p2 = std::make_shared<NoOpProcess>("P2");
    drift.setLocalProcesses({p1, p2});

    EXPECT_TRUE(drift.hasLocalProcesses());
    const auto& processes = drift.getLocalProcesses();
    ASSERT_EQ(processes.size(), 2u);
    EXPECT_EQ(processes[0].get(), p1.get());
    EXPECT_EQ(processes[1].get(), p2.get());
    EXPECT_EQ(processes[0]->getName(), "P1");
    EXPECT_EQ(processes[1]->getName(), "P2");
}

TEST(TestLocalProcessAttachment, AddLocalProcessAppends) {
    DriftRep drift("D1");

    auto p1 = std::make_shared<NoOpProcess>("P1");
    auto p2 = std::make_shared<NoOpProcess>("P2");
    drift.addLocalProcess(p1);
    drift.addLocalProcess(p2);

    const auto& processes = drift.getLocalProcesses();
    ASSERT_EQ(processes.size(), 2u);
    EXPECT_EQ(processes[0].get(), p1.get());
    EXPECT_EQ(processes[1].get(), p2.get());
}

TEST(TestLocalProcessAttachment, ProcessesSurviveClone) {
    DriftRep drift("D1");

    auto p1 = std::make_shared<NoOpProcess>("P1");
    drift.addLocalProcess(p1);

    std::unique_ptr<ElementBase> clone(drift.clone());
    ASSERT_NE(clone, nullptr);

    EXPECT_TRUE(clone->hasLocalProcesses());
    const auto& processes = clone->getLocalProcesses();
    ASSERT_EQ(processes.size(), 1u);
    EXPECT_EQ(processes[0].get(), p1.get());
}
