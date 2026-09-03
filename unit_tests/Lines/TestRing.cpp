#include <gtest/gtest.h>

#include <memory>
#include <sstream>
#include <string>

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
