//
// Tests for buildElement() and buildElements(): runtime elements built from a run description.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include <gtest/gtest.h>

#include <memory>
#include <sstream>
#include <string>

#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "AbsBeamline/Multipole.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "OpalConfigure/Configure.h"
#include "Track/BuildElements.h"
#include "Utilities/OpalException.h"
#include "Utility/Inform.h"

extern Inform* gmsg;

using namespace opalx::run;

namespace {

    AttributeEntry set(const std::string& name, AttributeValue value) {
        AttributeEntry entry;
        entry.name        = name;
        entry.hasValue    = true;
        entry.defaultUsed = false;
        entry.value       = std::move(value);
        return entry;
    }

    ElementDescription quadrupole(const std::string& name) {
        ElementDescription element;
        element.object.keyword = "QUADRUPOLE";
        element.object.name    = name;
        element.object.attributes.push_back(set("L", 0.1));
        element.object.attributes.push_back(set("K1", 2.5));
        element.object.attributes.push_back(set("ELEMEDGE", 0.3));
        element.elementType = "Multipole";
        return element;
    }

    class BuildElementsTest : public ::testing::Test {
    protected:
        void SetUp() override {
            gmsg = new Inform("", output_m, INFORM_ALL_NODES);
            Configure::configure();
        }

        void TearDown() override {
            OpalData::deleteInstance();
            delete gmsg;
            gmsg = nullptr;
        }

        std::ostringstream output_m;
    };

}  // namespace

TEST_F(BuildElementsTest, ElementMatchesTheOneOfTheInputObject) {
    // The same QUADRUPOLE as an input object, the way the parser makes it.
    std::unique_ptr<Object> input(OpalData::getInstance()->makeObject("QUADRUPOLE", "Q1"));
    auto* quad = dynamic_cast<Element*>(input.get());
    ASSERT_NE(quad, nullptr);
    Attributes::setReal(*quad->findAttribute("L"), 0.1);
    Attributes::setReal(*quad->findAttribute("K1"), 2.5);
    Attributes::setReal(*quad->findAttribute("ELEMEDGE"), 0.3);
    quad->update();
    const auto* expected = dynamic_cast<const Multipole*>(quad->getElement());
    ASSERT_NE(expected, nullptr);

    const auto built = buildElement(quadrupole("Q1"));
    const auto* got  = dynamic_cast<const Multipole*>(built.get());
    ASSERT_NE(got, nullptr);
    EXPECT_EQ(got->getName(), "Q1");
    EXPECT_EQ(got->getGeometry().getElementLength(), expected->getGeometry().getElementLength());
    EXPECT_EQ(got->getGeometry().getElementLength(), 0.1);
    EXPECT_EQ(got->getElementPosition(), expected->getElementPosition());
    EXPECT_EQ(got->getNormalComponent(1), expected->getNormalComponent(1));
    EXPECT_NE(got->getNormalComponent(1), 0.0);
}

TEST_F(BuildElementsTest, AttributesTheInputDidNotSetKeepTheKeywordDefault) {
    ElementDescription element               = quadrupole("Q2");
    element.object.attributes[1].defaultUsed = true;  // K1 = 2.5 was a default: not copied
    const auto built                         = buildElement(element);
    const auto* got                          = dynamic_cast<const Multipole*>(built.get());
    ASSERT_NE(got, nullptr);
    EXPECT_EQ(got->getNormalComponent(1), 0.0);
}

TEST_F(BuildElementsTest, UnknownKeywordOrFailedAttributeThrows) {
    ElementDescription unknown = quadrupole("Q3");
    unknown.object.keyword     = "NOT_A_KEYWORD";
    EXPECT_THROW(buildElement(unknown), OpalException);

    ElementDescription failed         = quadrupole("Q4");
    failed.object.attributes[1].error = "Can't find variable 'X'";
    EXPECT_THROW(buildElement(failed), OpalException);
}

TEST_F(BuildElementsTest, RingElementsGetTheRingMembership) {
    ElementDescription marker;
    marker.object.keyword = "MARKER";
    marker.object.name    = "#S";
    marker.elementType    = "Marker";

    LineDescription line;
    line.name        = "R1";
    line.topology    = "RING";
    line.elements    = {marker, quadrupole("Q5")};
    const auto built = buildElements(line);
    ASSERT_EQ(built.size(), 2u);
    EXPECT_EQ(built[0]->getBeamlineMembership().topology, BeamlineTopology::LINE);
    EXPECT_EQ(built[1]->getBeamlineMembership().topology, BeamlineTopology::RING);
    EXPECT_EQ(built[1]->getBeamlineMembership().ownerName, "R1");
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ippl::initialize(argc, argv);
    const int result = RUN_ALL_TESTS();
    endfieldmodel::EndFieldModelManager::clearEFMManager();
    ippl::finalize();
    return result;
}
