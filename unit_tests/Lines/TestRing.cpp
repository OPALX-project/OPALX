#include <gtest/gtest.h>

#include <memory>
#include <sstream>
#include <string>

#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/MarkerRep.h"
#include "Beamlines/FlaggedBeamline.h"
#include "Lines/Ring.h"
#include "OpalParser/SimpleStatement.h"
#include "OpalParser/Token.h"
#include "Utilities/ParseError.h"

namespace {
    constexpr char testFile[] = "TestRing.in";

    Token delimiter(char value) { return Token(testFile, 1, Token::IS_DELIMITER, value); }

    Token word(const char* value) { return Token(testFile, 1, Token::IS_WORD, value); }

    SimpleStatement statement(Statement::TokenList& tokens) {
        return SimpleStatement(testFile, tokens);
    }

    class TestableRing : public Ring {
    public:
        FlaggedBeamline* line() { return fetchLine(); }
    };
}  // namespace

class RingTest : public ::testing::Test {
protected:
    Ring* makeRing(const std::string& name) { return exemplar_m.clone(name); }

private:
    Ring exemplar_m;
};

TEST_F(RingTest, ClonesAndPrintsRingDefinition) {
    std::unique_ptr<Ring> ring(makeRing("TEST_RING"));

    std::ostringstream output;
    ring->print(output);
    EXPECT_EQ(output.str(), "TEST_RING:RING=();\n");
    EXPECT_EQ(ring->getElement()->getBeamlineTopology(), BeamlineTopology::RING);
    EXPECT_EQ(ring->getElement()->getBeamlineOwnerName(), "TEST_RING");
}

TEST(RingMembershipTest, DefaultsValidatesAndSurvivesClone) {
    DriftRep drift("D");
    EXPECT_EQ(drift.getBeamlineTopology(), BeamlineTopology::LINE);
    EXPECT_TRUE(drift.getBeamlineOwnerName().empty());

    EXPECT_THROW(drift.setBeamlineMembership(BeamlineTopology::RING), GeneralOpalException);
    EXPECT_THROW(
            drift.setBeamlineMembership(BeamlineTopology::LINE, "NOT_ALLOWED"),
            GeneralOpalException);

    drift.setBeamlineMembership(BeamlineTopology::RING, "R1");
    std::unique_ptr<ElementBase> clone(drift.clone());
    EXPECT_EQ(clone->getBeamlineTopology(), BeamlineTopology::RING);
    EXPECT_EQ(clone->getBeamlineOwnerName(), "R1");

    drift.clearBeamlineMembership();
    EXPECT_EQ(drift.getBeamlineTopology(), BeamlineTopology::LINE);
    EXPECT_TRUE(drift.getBeamlineOwnerName().empty());
}

TEST(RingMembershipTest, ClonesAndTagsOccurrencesRecursively) {
    auto prototype = std::make_shared<DriftRep>("D");
    auto marker    = std::make_shared<MarkerRep>("#S");
    auto nested    = std::make_shared<FlaggedBeamline>("SUBLINE");
    nested->push_back(FlaggedElmPtr(ElmPtr(prototype)));

    TestableRing ring;
    ring.line()->push_back(FlaggedElmPtr(ElmPtr(marker)));
    ring.line()->push_back(FlaggedElmPtr(ElmPtr(prototype)));
    ring.line()->push_back(FlaggedElmPtr(ElmPtr(prototype)));
    ring.line()->push_back(FlaggedElmPtr(ElmPtr(nested)));
    ring.prepareForTracking();

    auto member = ring.line()->begin();
    EXPECT_EQ(member->getElement(), marker.get());
    EXPECT_EQ(member->getElement()->getBeamlineTopology(), BeamlineTopology::LINE);

    ElementBase* first  = (++member)->getElement();
    ElementBase* second = (++member)->getElement();
    EXPECT_NE(first, prototype.get());
    EXPECT_NE(second, prototype.get());
    EXPECT_NE(first, second);
    EXPECT_EQ(first->getBeamlineOwnerName(), "RING");
    EXPECT_EQ(second->getBeamlineOwnerName(), "RING");
    EXPECT_EQ(prototype->getBeamlineTopology(), BeamlineTopology::LINE);

    auto* nestedOccurrence = dynamic_cast<FlaggedBeamline*>((++member)->getElement());
    ASSERT_NE(nestedOccurrence, nullptr);
    EXPECT_NE(nestedOccurrence, nested.get());
    EXPECT_EQ(nestedOccurrence->getBeamlineTopology(), BeamlineTopology::RING);
    EXPECT_EQ(nestedOccurrence->getBeamlineOwnerName(), "RING");
    ASSERT_EQ(nestedOccurrence->size(), 1);
    EXPECT_NE(nestedOccurrence->front().getElement(), prototype.get());
    EXPECT_EQ(
            nestedOccurrence->front().getElement()->getBeamlineTopology(), BeamlineTopology::RING);
    EXPECT_EQ(nestedOccurrence->front().getElement()->getBeamlineOwnerName(), "RING");
}

TEST(RingMembershipTest, KeepsSharedPrototypeIndependentAcrossRings) {
    auto prototype = std::make_shared<DriftRep>("SHARED");
    TestableRing first;
    TestableRing second;
    first.setOpalName("R1");
    second.setOpalName("R2");
    first.line()->push_back(FlaggedElmPtr(ElmPtr(prototype)));
    second.line()->push_back(FlaggedElmPtr(ElmPtr(prototype)));

    first.prepareForTracking();
    second.prepareForTracking();

    ElementBase* firstOccurrence  = first.line()->front().getElement();
    ElementBase* secondOccurrence = second.line()->front().getElement();
    EXPECT_NE(firstOccurrence, secondOccurrence);
    EXPECT_NE(firstOccurrence, prototype.get());
    EXPECT_NE(secondOccurrence, prototype.get());
    EXPECT_EQ(firstOccurrence->getBeamlineOwnerName(), "R1");
    EXPECT_EQ(secondOccurrence->getBeamlineOwnerName(), "R2");
    EXPECT_EQ(prototype->getBeamlineTopology(), BeamlineTopology::LINE);
    EXPECT_TRUE(prototype->getBeamlineOwnerName().empty());
}

TEST_F(RingTest, RejectsReflection) {
    Statement::TokenList tokens{
            delimiter('='), delimiter('('), delimiter('-'), word("DRIFT"), delimiter(')')};
    SimpleStatement input = statement(tokens);
    std::unique_ptr<Ring> ring(makeRing("REFLECTED_RING"));

    EXPECT_THROW(ring->parse(input), ParseError);
}

TEST_F(RingTest, RejectsRepetition) {
    Statement::TokenList tokens{delimiter('='), delimiter('('), Token(testFile, 1, "2", 2),
                                delimiter('*'), word("DRIFT"),  delimiter(')')};
    SimpleStatement input = statement(tokens);
    std::unique_ptr<Ring> ring(makeRing("REPEATED_RING"));

    EXPECT_THROW(ring->parse(input), ParseError);
}
