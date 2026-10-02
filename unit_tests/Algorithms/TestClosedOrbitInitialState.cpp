#include "Algorithms/ClosedOrbitInitialState.h"
#include "Attributes/Attributes.h"
#include "Track/CofCmd.h"
#include "Utilities/Util.h"
#include "gtest/gtest.h"

using V = ippl::Vector<double, 3>;

TEST(ClosedOrbitInitialState, FrameMapsOriginAndFullMomentum) {
    ClosedOrbitInitialState s;
    s.position = V(.12, -.03, 2.4);
    s.section  = CoordinateSystemTrafo(V(1, 2, 3), Quaternion(std::cos(.3), 0, std::sin(.3), 0));
    s.momentum = s.section.rotateFrom(V(.12, -.08, 1.7));
    const auto frame = s.frame();
    const V origin   = frame.transformFrom(V(0.0));
    const V p        = frame.rotateFrom(V(0, 0, euclidean_norm(s.momentum)));
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_NEAR(origin[d], s.position[d], 1e-15);
        EXPECT_NEAR(p[d], s.momentum[d], 1e-15);
    }
    const V offset(.01, -.02, .03), localP(.02, .01, 1.8);
    const V placed = frame.transformFrom(offset), rotated = frame.rotateFrom(localP);
    EXPECT_NEAR(euclidean_norm(V(placed - s.position)), euclidean_norm(offset), 1e-15);
    EXPECT_NEAR(euclidean_norm(rotated), euclidean_norm(localP), 1e-15);
    const auto back = frame.transformTo(placed);
    for (unsigned d = 0; d < 3; ++d)
        EXPECT_NEAR(back[d], offset[d], 1e-15);
}

TEST(ClosedOrbitInitialState, RejectInvalidState) {
    ClosedOrbitInitialState s;
    EXPECT_THROW(s.frame(), OpalException);
    s.momentum = V(0, 0, -1);
    EXPECT_THROW(s.frame(), OpalException);
    s.momentum    = V(0, 0, 1);
    s.position[0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(s.frame(), OpalException);
    s.position[0] = 0;
    s.time        = std::numeric_limits<double>::infinity();
    EXPECT_THROW(s.frame(), OpalException);
}

TEST(ClosedOrbitInitialState, CompatibilityAndIndependentSnapshots) {
    PartData ref;
    ref.setM(1e9);
    ref.setQ(1);
    ref.setP(2e9);
    ClosedOrbitInitialState s;
    s.massEV   = ref.getM();
    s.chargeE  = ref.getQ();
    s.momentum = V(0, 0, 2);
    s.line     = "RING1";
    s.species  = "PROTON";
    EXPECT_NO_THROW(s.validate("RING1", "PROTON", ref));
    EXPECT_THROW(s.validate("RING2", "PROTON", ref), OpalException);
    EXPECT_THROW(s.validate("RING1", "ELECTRON", ref), OpalException);
    auto other        = s;
    other.position[0] = .1;
    other.momentum[2] = 3;
    EXPECT_DOUBLE_EQ(s.position[0], 0);
    EXPECT_DOUBLE_EQ(s.momentum[2], 2);
    EXPECT_THROW(other.validate("RING1", "PROTON", ref), OpalException);
    ref.setQ(-1);
    EXPECT_THROW(s.validate("RING1", "PROTON", ref), OpalException);
}

TEST(ClosedOrbitInitialState, CofAttributesAreOnSingleCommand) {
    CofCmd command;
    for (const auto* name :
         {"LINE", "BEAM", "DT", "X", "PX", "Y", "PY", "XTOL", "PTOL", "OUTPUT", "FDSTEP", "SCALES"})
        EXPECT_NE(command.findAttribute(name), nullptr) << name;
    EXPECT_TRUE(command.findAttribute("OUTPUT")->defaultUsed());
    EXPECT_TRUE(Attributes::getString(*command.findAttribute("OUTPUT")).empty());
    std::unique_ptr<CofCmd> copy(command.clone("IC"));
    EXPECT_EQ(copy->getOpalName(), "IC");
    EXPECT_TRUE(copy->findAttribute("OUTPUT")->defaultUsed());
}
