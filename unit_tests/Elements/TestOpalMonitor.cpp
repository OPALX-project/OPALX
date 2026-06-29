//
// Unit tests for MONITOR element configuration and geometric support.
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

#include "Attributes/Attributes.h"
#include "BeamlineCore/MonitorRep.h"
#include "Elements/OpalMonitor.h"
#include "gtest/gtest.h"

extern Inform* gmsg;

class TestOpalMonitor : public testing::Test {
public:
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
};

TEST_F(TestOpalMonitor, UpdateAllowsZeroBodyLengthAndPreservesSamplingSupport) {
    OpalMonitor ui;
    std::unique_ptr<OpalMonitor> instance(ui.clone("MON0"));
    Attributes::setReal(instance->itsAttr[OpalMonitor::LENGTH], 0.0);
    Attributes::setPredefinedString(instance->itsAttr[OpalMonitor::TYPE], "TEMPORAL");
    Attributes::setString(instance->itsAttr[OpalMonitor::OUTFN], "mon_zero");

    EXPECT_NO_THROW(instance->update());

    const auto* monitor = dynamic_cast<const MonitorRep*>(instance->getElement());
    ASSERT_NE(monitor, nullptr);
    EXPECT_NEAR(monitor->getElementLength(), 0.0, 1.0e-12);

    EXPECT_TRUE(monitor->isInside(ippl::Vector<double, 3>(0.0, 0.0, 0.0049)));
    EXPECT_FALSE(monitor->isInside(ippl::Vector<double, 3>(0.0, 0.0, 0.0051)));
}

TEST_F(TestOpalMonitor, ExplicitBodyLengthExpandsMonitorSupportInterval) {
    OpalMonitor ui;
    std::unique_ptr<OpalMonitor> instance(ui.clone("MON1"));
    Attributes::setReal(instance->itsAttr[OpalMonitor::LENGTH], 0.03);
    Attributes::setPredefinedString(instance->itsAttr[OpalMonitor::TYPE], "SPATIAL");
    Attributes::setString(instance->itsAttr[OpalMonitor::OUTFN], "mon_long");

    EXPECT_NO_THROW(instance->update());

    const auto* monitor = dynamic_cast<const MonitorRep*>(instance->getElement());
    ASSERT_NE(monitor, nullptr);
    EXPECT_NEAR(monitor->getElementLength(), 0.03, 1.0e-12);

    EXPECT_TRUE(monitor->isInside(ippl::Vector<double, 3>(0.0, 0.0, 0.014)));
    EXPECT_FALSE(monitor->isInside(ippl::Vector<double, 3>(0.0, 0.0, 0.016)));
}
