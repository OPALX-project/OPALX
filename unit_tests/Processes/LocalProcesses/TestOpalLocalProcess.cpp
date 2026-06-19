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
#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalDrift.h"
#include "Processes/LocalProcesses/OpalTestProcess.h"
#include "Processes/LocalProcesses/TestCounterProcess.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

#include <memory>
#include <string>
#include <vector>

namespace {

    /// Register a TESTPROCESS definition clone under the given name in OpalData.
    OpalTestProcess* defineTestProcess(const std::string& name) {
        OpalTestProcess exemplar;
        OpalTestProcess* definition = exemplar.clone(name);
        OpalData::getInstance()->define(definition);
        return definition;
    }

    std::unique_ptr<OpalDrift> makeDrift(const std::vector<std::string>& processNames) {
        auto drift = std::make_unique<OpalDrift>();
        Attributes::setReal(*drift->findAttribute("L"), 2.0);
        Attributes::setReal(*drift->findAttribute("ELEMEDGE"), 0.0);
        if (!processNames.empty()) {
            Attributes::setUpperCaseStringArray(*drift->findAttribute("PROCESSES"), processNames);
        }
        drift->update();
        return drift;
    }

}  // namespace

TEST(TestOpalLocalProcess, CloneProvidesCachedProcessInstance) {
    OpalTestProcess exemplar;
    std::unique_ptr<OpalTestProcess> definition(exemplar.clone("TPCACHED"));

    std::shared_ptr<LocalProcess> process = definition->getProcess();
    ASSERT_NE(process, nullptr);
    EXPECT_EQ(process->getName(), "TPCACHED");
    EXPECT_NE(dynamic_cast<TestCounterProcess*>(process.get()), nullptr);

    // Repeated calls must hand out the same runtime instance.
    EXPECT_EQ(definition->getProcess().get(), process.get());
}

TEST(TestOpalLocalProcess, FindResolvesRegisteredDefinition) {
    OpalTestProcess* definition = defineTestProcess("TPFIND");

    EXPECT_EQ(OpalLocalProcess::find("TPFIND"), definition);
    EXPECT_THROW(OpalLocalProcess::find("TPDOESNOTEXIST"), OpalException);
}

TEST(TestOpalLocalProcess, ElementUpdateAttachesProcesses) {
    OpalTestProcess* definition = defineTestProcess("TPATTACH");

    auto drift = makeDrift({"TPATTACH"});

    ElementBase* element = drift->getElement();
    ASSERT_NE(element, nullptr);
    ASSERT_TRUE(element->hasLocalProcesses());

    const auto& processes = element->getLocalProcesses();
    ASSERT_EQ(processes.size(), 1u);
    EXPECT_EQ(processes[0].get(), definition->getProcess().get());
}

TEST(TestOpalLocalProcess, UnknownProcessNameThrowsOnUpdate) {
    EXPECT_THROW(makeDrift({"TPUNKNOWN"}), OpalException);
}

TEST(TestOpalLocalProcess, RedefinitionIsRejected) {
    defineTestProcess("TPREDEF");

    OpalTestProcess exemplar;
    OpalTestProcess* replacement = exemplar.clone("TPREDEF");

    bool threw = false;
    try {
        OpalData::getInstance()->define(replacement);
    } catch (const OpalException&) {
        // define() does not take ownership when it rejects the object.
        threw = true;
        delete replacement;
    }
    EXPECT_TRUE(threw);
}
