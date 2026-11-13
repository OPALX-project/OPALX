#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <sstream>
#include <string>

#include "../KeywordObjects.h"
#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/MarkerRep.h"
#include "BeamlineCore/SBendRep.h"
#include "Ippl.h"
#include "Lines/Ring.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/OpalParser.h"
#include "OpalParser/SimpleStatement.h"
#include "OpalParser/StringStream.h"
#include "OpalParser/Token.h"
#include "Track/BuildElements.h"
#include "Track/FreezeRun.h"
#include "Utilities/OpalException.h"
#include "Utilities/ParseError.h"
#include "Utility/Inform.h"

extern Inform* gmsg;

namespace {
    constexpr char testFile[] = "TestRing.in";

    Token delimiter(char value) { return Token(testFile, 1, Token::IS_DELIMITER, value); }

    Token word(const char* value) { return Token(testFile, 1, Token::IS_WORD, value); }

    SimpleStatement statement(Statement::TokenList& tokens) {
        return SimpleStatement(testFile, tokens);
    }
}  // namespace

class RingTest : public ::testing::Test {
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

    Ring* makeRing(const std::string& name) {
        return opalx::test::makeKeywordObject<Ring>("RING", name).release();
    }

    void run(const std::string& text) {
        std::string line = text;
        for (char& c : line) {
            if (c == '\n') c = ' ';
        }
        OpalParser parser;
        parser.run(new StringStream(line));
    }

    /// The elements built for the occurrences of the parsed LINE or RING @p name.
    std::vector<std::shared_ptr<ElementBase>> build(const std::string& name) {
        return opalx::run::buildElements(opalx::run::describeLine(*BeamSequence::find(name)));
    }

    std::ostringstream output_m;
};

TEST_F(RingTest, PrintsRingDefinitionAndIsDescribedAsRing) {
    std::unique_ptr<Ring> ring(makeRing("TEST_RING"));

    std::ostringstream output;
    ring->print(output);
    EXPECT_EQ(output.str(), "TEST_RING:RING=();\n");
    const opalx::run::LineDescription line = opalx::run::describeLine(*ring);
    EXPECT_EQ(line.topology, "RING");
    EXPECT_EQ(line.name, "TEST_RING");
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

TEST_F(RingTest, BuildsAndTagsOneElementPerOccurrenceRecursively) {
    run(R"(
        D: DRIFT, L = 0.5, ELEMEDGE = 0;
        SUBLINE: LINE = (D);
        R: RING = (D, D, SUBLINE);
    )");
    // #S, D, D, then SUBLINE with its own #S and #E, then #E.
    const auto elements = build("R");
    ASSERT_EQ(elements.size(), 7u);
    // D built on its own, outside any RING: not tagged.
    const auto definition    = opalx::test::buildFrom(*Element::find("D"));
    const ElementBase* input = definition.get();

    EXPECT_EQ(elements[0]->getName(), "#S");
    EXPECT_EQ(elements[0]->getBeamlineTopology(), BeamlineTopology::LINE);

    const ElementBase* first  = elements[1].get();
    const ElementBase* second = elements[2].get();
    EXPECT_NE(first, input);
    EXPECT_NE(second, input);
    EXPECT_NE(first, second);
    EXPECT_EQ(first->getBeamlineOwnerName(), "R");
    EXPECT_EQ(second->getBeamlineOwnerName(), "R");
    EXPECT_EQ(input->getBeamlineTopology(), BeamlineTopology::LINE);

    // The member of the nested LINE carries the RING name, not the LINE name.
    const ElementBase* nested = elements[4].get();
    EXPECT_EQ(nested->getName(), "D");
    EXPECT_NE(nested, input);
    EXPECT_EQ(nested->getBeamlineTopology(), BeamlineTopology::RING);
    EXPECT_EQ(nested->getBeamlineOwnerName(), "R");
}

TEST_F(RingTest, KeepsSharedElementIndependentAcrossRings) {
    run(R"(
        DS: DRIFT, L = 0.5, ELEMEDGE = 0;
        R1: RING = (DS);
        R2: RING = (DS);
    )");
    const auto first  = build("R1");
    const auto second = build("R2");
    ASSERT_EQ(first.size(), 3u);
    ASSERT_EQ(second.size(), 3u);
    // DS built on its own, outside any RING: not tagged.
    const auto definition    = opalx::test::buildFrom(*Element::find("DS"));
    const ElementBase* input = definition.get();

    EXPECT_NE(first[1].get(), second[1].get());
    EXPECT_NE(first[1].get(), input);
    EXPECT_NE(second[1].get(), input);
    EXPECT_EQ(first[1]->getBeamlineOwnerName(), "R1");
    EXPECT_EQ(second[1]->getBeamlineOwnerName(), "R2");
    EXPECT_EQ(input->getBeamlineTopology(), BeamlineTopology::LINE);
    EXPECT_TRUE(input->getBeamlineOwnerName().empty());
}

TEST(RingMembershipTest, CircumferenceCountsNominalArcsAndDriftsNotFringeSupports) {
    auto bend = std::make_shared<SBendRep>("B");
    bend->getGeometry().setElementLength(1.0);
    bend->getGeometry().setCurvature(0.5);
    auto drift = std::make_shared<DriftRep>("D");
    drift->getGeometry().setElementLength(0.4);
    auto marker = std::make_shared<MarkerRep>("#S");
    // Two occurrences of a cell (bend, drift); markers do not count.
    const std::vector<std::shared_ptr<ElementBase>> ring{marker, bend, drift, bend, drift};
    for (const double gap : {0.0, 0.1, 0.2}) {
        bend->setFullGap(gap);
        bend->setFringeIntegral(0.1);
        EXPECT_DOUBLE_EQ(opalx::run::ringLength(ring), 2.8);
        double begin = 0.0, end = 0.0;
        bend->getFieldExtent(begin, end);
        EXPECT_DOUBLE_EQ(end - begin, 1.0 + 10.0 * gap);
        EXPECT_EQ(bend->isInside(Vector_t<double, 3>(0, 0, -0.01)), gap > 0.0);
    }
    // Rectangular-bend geometry stores a chord, not a design arc.
    drift->getGeometry() = Geometry::makeRBend(0.4, 0.6);
    EXPECT_NEAR(opalx::run::ringLength(ring), 2.0 * (1.0 + 0.4 * 0.3 / std::sin(0.3)), 1.e-14);
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

TEST_F(RingTest, LineMacroIsCopiedInWithItsMarkersAndReflection) {
    run(R"(
        D: DRIFT, L = 0.5, ELEMEDGE = 0;
        Q: QUADRUPOLE, L = 0.1, ELEMEDGE = 0;
        M(X): LINE = (X, -D);
        L: LINE = (-M(Q));
    )");
    const opalx::run::LineDescription line = opalx::run::describeLine(*BeamSequence::find("L"));
    // #S, then the expansion of M with its own #S and #E, in order, every member reflected
    // unless it was reflected in M, then #E.
    const std::vector<std::pair<std::string, bool>> expected{
            {"#S", false}, {"#S", true}, {"Q", true}, {"D", false}, {"#E", true}, {"#E", false}};
    ASSERT_EQ(line.elements.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        EXPECT_EQ(line.elements[i].object.name, expected[i].first) << i;
        EXPECT_EQ(line.elements[i].reflected, expected[i].second) << i;
    }
}

TEST_F(RingTest, LineThatContainsItselfIsRefused) {
    run(R"(
        D: DRIFT, L = 0.5, ELEMEDGE = 0;
        L: LINE = (D);
        L: LINE = (D, L);
    )");
    EXPECT_THROW(opalx::run::describeLine(*BeamSequence::find("L")), OpalException);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ippl::initialize(argc, argv);
    const int result = RUN_ALL_TESTS();
    endfieldmodel::EndFieldModelManager::clearEFMManager();
    ippl::finalize();
    return result;
}
