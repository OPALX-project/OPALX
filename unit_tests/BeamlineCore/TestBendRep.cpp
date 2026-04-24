/**
 * @file TestBendRep.cpp
 * @brief Regression tests for the OPALX-native SBEND/RBEND representation layer.
 *
 * The tests cover the first port of the analytic bend representation:
 * geometry ownership, dipole-field channels, and the curved placement hooks
 * used by the element-placement bridge.
 */

#include <gtest/gtest.h>

#include "BeamlineCore/RBendRep.h"
#include "BeamlineCore/SBendRep.h"
#include "Channels/Channel.h"

#include <cmath>
#include <memory>

namespace {

    TEST(SBendRep, GeometryAndDesignPathFollowPlanarArc) {
        SBendRep bend("SBEND");
        bend.getGeometry() = PlanarArcGeometry(2.0, M_PI / 12.0);
        bend.setElementLength(2.0);
        bend.setBendAngle(M_PI / 6.0);
        bend.setExitAngle(M_PI / 18.0);

        EXPECT_EQ(bend.getType(), ElementType::SBEND);
        EXPECT_NEAR(bend.getGeometry().getElementLength(), 2.0, 1e-12);
        EXPECT_NEAR(bend.getBendAngle(), M_PI / 6.0, 1e-12);

        const double radius = 2.0 / (M_PI / 6.0);
        const double chord  = 2.0 * radius * std::sin(M_PI / 12.0);
        EXPECT_NEAR(bend.getChordLength(), chord, 1e-12);

        const auto path = bend.getDesignPath(7);
        ASSERT_GE(path.size(), 7u);

        const auto entry = bend.getEdgeToBegin().getOrigin();
        const auto exit  = bend.getEdgeToEnd().getOrigin();
        EXPECT_NEAR(path.front()(0), entry(0), 1e-12);
        EXPECT_NEAR(path.front()(1), entry(1), 1e-12);
        EXPECT_NEAR(path.front()(2), entry(2), 1e-12);
        EXPECT_NEAR(path.back()(0), exit(0), 1e-12);
        EXPECT_NEAR(path.back()(1), exit(1), 1e-12);
        EXPECT_NEAR(path.back()(2), exit(2), 1e-12);
    }

    TEST(SBendRep, DipoleChannelUsesCurrentZeroBasedMultipoleConvention) {
        SBendRep bend("SBEND");

        Channel* by = bend.getChannel("BY", false);
        ASSERT_NE(by, nullptr);
        EXPECT_TRUE(by->set(1.75));
        EXPECT_DOUBLE_EQ(bend.getB(), 1.75);
        EXPECT_DOUBLE_EQ(bend.getField().getNormalComponent(0), 1.75);
        delete by;
    }

    TEST(RBendRep, ExitAngleTracksRectangularBendConvention) {
        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(1.2);
        bend.getGeometry().setBendAngle(M_PI / 9.0);
        bend.setElementLength(1.2);
        bend.setBendAngle(M_PI / 9.0);
        bend.setEntranceAngle(M_PI / 18.0);

        EXPECT_EQ(bend.getType(), ElementType::RBEND);
        EXPECT_NEAR(bend.getExitAngle(), M_PI / 18.0, 1e-12);
        EXPECT_NEAR(bend.getChordLength(), 1.2, 1e-12);

        const auto path = bend.getDesignPath(5);
        ASSERT_GE(path.size(), 5u);
        const auto entry = bend.getEdgeToBegin().getOrigin();
        const auto exit  = bend.getEdgeToEnd().getOrigin();
        EXPECT_NEAR(path.front()(0), entry(0), 1e-12);
        EXPECT_NEAR(path.back()(0), exit(0), 1e-12);
    }

    TEST(RBendRep, ClonePreservesGeometryAndField) {
        RBendRep bend("RBEND");
        bend.getGeometry().setElementLength(0.8);
        bend.getGeometry().setBendAngle(M_PI / 8.0);
        bend.setElementLength(0.8);
        bend.setBendAngle(M_PI / 8.0);
        bend.setB(0.9);

        std::unique_ptr<ElementBase> clone(bend.clone());
        auto* copy = dynamic_cast<RBendRep*>(clone.get());
        ASSERT_NE(copy, nullptr);
        EXPECT_NEAR(copy->getElementLength(), 0.8, 1e-12);
        EXPECT_NEAR(copy->getBendAngle(), M_PI / 8.0, 1e-12);
        EXPECT_NEAR(copy->getB(), 0.9, 1e-12);
    }

}  // namespace
