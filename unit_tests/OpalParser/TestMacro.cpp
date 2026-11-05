#include <gtest/gtest.h>

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AbstractObjects/OpalData.h"
#include "BasicActions/Value.h"
#include "Ippl.h"
#include "OpalParser/CompoundStatement.h"
#include "OpalParser/MacroCmd.h"
#include "OpalParser/MacroStream.h"
#include "OpalParser/OpalParser.h"
#include "OpalParser/SimpleStatement.h"
#include "OpalParser/StringStream.h"
#include "Utilities/ParseError.h"
#include "ValueDefinitions/RealConstant.h"
#include "ValueDefinitions/RealVariable.h"
#include "ValueDefinitions/RealVector.h"
#include "ValueDefinitions/StringConstant.h"

extern Inform* gmsg;

namespace {
    // Tokenize a single statement without the parser's input stack. In factory
    // tests the first token is the first formal/actual, just after the opening '('.
    std::unique_ptr<SimpleStatement> statement(const std::string& text) {
        auto result = std::make_unique<SimpleStatement>("TestMacro", 1);
        StringStream source(text);
        for (Token token = source.readToken(); !token.isEOF(); token = source.readToken()) {
            result->append(token);
        }
        result->start();
        return result;
    }

    class RecordingParser : public Parser {
    public:
        mutable std::vector<std::vector<std::string>> calls;

        void parse(Statement&) const override {}
        Statement* readStatement(TokenStream*) const override { return nullptr; }
        void run(TokenStream* input) const override {
            std::unique_ptr<TokenStream> owned(input);
            std::vector<std::string> tokens;
            for (Token token = input->readToken(); !token.isEOF(); token = input->readToken()) {
                tokens.push_back(token.getLex());
            }
            calls.push_back(tokens);
        }
    };

    std::unique_ptr<Object> makeTemplate(MacroCmd& model, const std::string& text) {
        auto stat = statement(text);
        StringStream unused("");
        return std::unique_ptr<Object>(model.makeTemplate("TEST", unused, *stat));
    }

    std::unique_ptr<Object> makeInstance(
            Object& model, const std::string& actuals, const Parser& parser) {
        auto stat = statement(actuals);
        return std::unique_ptr<Object>(model.makeInstance("TEST", *stat, &parser));
    }

    class CountingMacro : public MacroCmd {
    public:
        CountingMacro(const std::string& name, MacroCmd* parent, int& destroyed, bool fail)
            : MacroCmd(name, parent), destroyed_m(destroyed), fail_m(fail) {}
        ~CountingMacro() override { ++destroyed_m; }
        void execute() override {
            if (fail_m) throw std::runtime_error("test execution failure");
        }

    private:
        int& destroyed_m;
        bool fail_m;
    };

    class CountingModel : public MacroCmd {
    public:
        int destroyed{0};
        bool fail{false};
        Object* makeInstance(const std::string& name, Statement&, const Parser*) override {
            return new CountingMacro(name, this, destroyed, fail);
        }
    };

    class CountingStream : public StringStream {
    public:
        explicit CountingStream(int& destroyed) : StringStream(""), destroyed_m(destroyed) {}
        ~CountingStream() override { ++destroyed_m; }

    private:
        int& destroyed_m;
    };

    class ThrowingParser : public OpalParser {
    public:
        using OpalParser::run;
        void run() const override { throw std::runtime_error("test parser failure"); }
    };

    class MacroTest : public ::testing::Test {
    protected:
        void SetUp() override {
            gmsg = new Inform("", output_m, INFORM_ALL_NODES);
            gmsg->setOutputLevel(5);
            auto* data = OpalData::getInstance();
            data->create(new RealConstant());
            data->create(new RealVariable());
            data->create(new RealVector());
            data->create(new StringConstant());
            data->create(new MacroCmd());
            data->create(new Value());
        }

        void TearDown() override {
            OpalData::deleteInstance();
            delete gmsg;
            gmsg = nullptr;
        }

        void run(const std::string& text) {
            // StringStream tokenizes a single line; fold manuscript-style input
            // lines into spaces without changing any token boundaries.
            std::string line = text;
            for (char& c : line) {
                if (c == '\n') c = ' ';
            }
            parser_m.run(new StringStream(line));
        }

        double real(const std::string& name) {
            auto* value = dynamic_cast<RealVariable*>(OpalData::getInstance()->find(name));
            if (!value) throw std::runtime_error("Missing real variable " + name);
            return value->getReal();
        }

        OpalParser parser_m;
        std::ostringstream output_m;
    };
}  // namespace

TEST_F(MacroTest, Issue583MinimalReproducer) {
    run("TEST_MACRO(VALUE): MACRO { REAL X = 1; } TEST_MACRO(0);");
    EXPECT_DOUBLE_EQ(real("X"), 1.0);
}

TEST_F(MacroTest, DocumentedShowitAndKsExamples) {
    run(R"(
        SHOWIT(X): MACRO { VALUE, VALUE={X}; }
        KS(X,Y): MACRO { Y = 3e-3*X + 0.5e-3; }
        SHOWIT(PI);
        REAL X = 2;
        REAL Y;
        KS(X,Y);
    )");
    EXPECT_DOUBLE_EQ(real("Y"), 0.0065);
    EXPECT_NE(output_m.str().find("3.14159"), std::string::npos);
}

TEST_F(MacroTest, ZeroArgumentsEmptyBodyAndRepeatedCalls) {
    run(R"(
        REAL COUNT = 0;
        TICK(): MACRO { COUNT = COUNT + 1; }
        EMPTY(): MACRO {}
        TICK(); EMPTY(); TICK(); TICK();
    )");
    EXPECT_DOUBLE_EQ(real("COUNT"), 3.0);
}

TEST_F(MacroTest, NestedCallsAndCompoundBlocksAreReusable) {
    run(R"(
        REAL TOTAL = 0;
        ADD(X): MACRO { TOTAL = TOTAL + (X); }
        TWICE(X): MACRO { { ADD(X); ADD(X); } }
        TWICE(2); TWICE(3);
        REAL I = 0;
        WHILE (I < 3) { ADD(1); I = I + 1; }
        IF (TOTAL == 13) { ADD(2); } ELSE { ADD(100); }
    )");
    EXPECT_DOUBLE_EQ(real("TOTAL"), 15.0);
    EXPECT_DOUBLE_EQ(real("I"), 3.0);
}

TEST_F(MacroTest, VectorAndNestedFunctionArguments) {
    run(R"(
        SAVE(V): MACRO { VECTOR SAVED = V; }
        SAVE({1, MAX(2,3), 4});
        REAL RESULT;
        SET(X,Y): MACRO { Y = X; }
        SET(MAX(2, MIN(5,3)), RESULT);
    )");
    auto* saved = dynamic_cast<RealVector*>(OpalData::getInstance()->find("SAVED"));
    ASSERT_NE(saved, nullptr);
    EXPECT_DOUBLE_EQ(saved->getRealComponent(1), 1.0);
    EXPECT_DOUBLE_EQ(saved->getRealComponent(2), 3.0);
    EXPECT_DOUBLE_EQ(saved->getRealComponent(3), 4.0);
    EXPECT_DOUBLE_EQ(real("RESULT"), 3.0);
}

TEST_F(MacroTest, SubstitutionPreservesExpressionPrecedence) {
    run(R"(
        RAW(X,Y): MACRO { Y = 2*X; }
        GROUPED(X,Y): MACRO { Y = 2*(X); }
        REAL A; REAL B;
        RAW(1+2,A); GROUPED(1+2,B);
    )");
    EXPECT_DOUBLE_EQ(real("A"), 4.0);
    EXPECT_DOUBLE_EQ(real("B"), 6.0);
}

TEST_F(MacroTest, SuccessfulAndThrowingInstancesAreDestroyed) {
    auto* model = new CountingModel();
    model->setOpalName("COUNTED");
    OpalData::getInstance()->create(model);
    run("COUNTED(); COUNTED();");
    EXPECT_EQ(model->destroyed, 2);
    model->fail = true;
    auto stat   = statement("COUNTED()");
    EXPECT_THROW(parser_m.parse(*stat), std::runtime_error);
    EXPECT_EQ(model->destroyed, 3);
    EXPECT_THROW(run("COUNTED();"), std::runtime_error);
    EXPECT_EQ(model->destroyed, 4);
    EXPECT_TRUE(OpalParser::readToken().isEOF());
    model->fail = false;
    run("COUNTED();");
    EXPECT_EQ(model->destroyed, 5);
}

TEST_F(MacroTest, DirectCallsRejectUnknownMacrosAndUnclosedArgumentLists) {
    for (const std::string text : {"UNKNOWN()", "UNKNOWN(1", "UNKNOWN((1)"}) {
        auto stat = statement(text);
        EXPECT_THROW(parser_m.parse(*stat), ParseError);
    }
}

TEST(MacroFactoryTest, ReplacesWholeWordsOnlyAndPreservesQuotedStrings) {
    MacroCmd model;
    RecordingParser parser;
    auto macro    = makeTemplate(model, R"(X): MACRO { X XX "X" 'X'; })");
    auto instance = makeInstance(*macro, "42)", parser);
    instance->execute();
    instance->execute();
    ASSERT_EQ(parser.calls.size(), 2u);
    EXPECT_EQ(parser.calls[0], (std::vector<std::string>{"42", "XX", "X", "X", ";"}));
    EXPECT_EQ(parser.calls[0], parser.calls[1]);
}

TEST(MacroFactoryTest, DelimitersInsideBracesBracketsAndStringsStayWithinArguments) {
    MacroCmd model;
    RecordingParser parser;
    auto macro    = makeTemplate(model, "X,Y): MACRO { X; Y; }");
    auto instance = makeInstance(*macro, R"(A[1,2], {"a,b", F(3,4)}))", parser);
    instance->execute();
    ASSERT_EQ(parser.calls.size(), 1u);
    EXPECT_EQ(
            parser.calls[0], (std::vector<std::string>{
                                     "A", "[", "1", ",", "2", "]", ";", "{", "a,b", ",", "F", "(",
                                     "3", ",", "4", ")", "}", ";"}));
}

TEST(MacroFactoryTest, RepeatedInstantiationsDoNotMutateTemplateOrResubstituteActuals) {
    MacroCmd model;
    RecordingParser parser;
    auto macro  = makeTemplate(model, "X,Y): MACRO { X; Y; }");
    auto first  = makeInstance(*macro, "Y,X)", parser);
    auto second = makeInstance(*macro, "1,2)", parser);
    first->execute();
    second->execute();
    first->execute();
    ASSERT_EQ(parser.calls.size(), 3u);
    EXPECT_EQ(parser.calls[0], (std::vector<std::string>{"Y", ";", "X", ";"}));
    EXPECT_EQ(parser.calls[1], (std::vector<std::string>{"1", ";", "2", ";"}));
    EXPECT_EQ(parser.calls[0], parser.calls[2]);
}

TEST(MacroFactoryTest, RejectsWrongArityEmptyArgumentsAndBadDelimiters) {
    MacroCmd model;
    RecordingParser parser;
    auto macro = makeTemplate(model, "X,Y): MACRO { X; Y; }");
    for (const std::string text : {")", "1)", "1,2,3)", ",2)", "1,)", "1,,2)", "1,2", "{1,2],3)"}) {
        SCOPED_TRACE(text);
        EXPECT_THROW(makeInstance(*macro, text, parser), ParseError);
    }
    auto valid = makeInstance(*macro, "1,2)", parser);
    EXPECT_NO_THROW(valid->execute());
}

TEST(MacroFactoryTest, RejectsDuplicateAndMalformedFormalsAndBodies) {
    MacroCmd model;
    for (const std::string text :
         {"X,X): MACRO {}", "X,): MACRO {}", "X): MACRO", "X): MACRO { X;"}) {
        SCOPED_TRACE(text);
        EXPECT_THROW(makeTemplate(model, text), ParseError);
    }
    EXPECT_THROW(model.execute(), ParseError);
    EXPECT_THROW(model.clone("INVALID"), ParseError);
    auto macro = makeTemplate(model, "): MACRO {}");
    EXPECT_THROW(macro->execute(), ParseError);
    EXPECT_EQ(macro->getCategory(), "MACRO");
    EXPECT_FALSE(macro->shouldTrace());
    EXPECT_FALSE(macro->shouldUpdate());
}

TEST(MacroStreamTest, EmptyRewindAndPutBack) {
    MacroStream stream("test");
    stream.start();
    EXPECT_TRUE(stream.readToken().isEOF());
    Token a("test", 1, Token::IS_WORD, "A");
    Token b("test", 1, Token::IS_WORD, "B");
    stream.append(a);
    stream.append(b);
    stream.start();
    EXPECT_EQ(stream.readToken().getWord(), "A");
    stream.putBack(a);
    EXPECT_EQ(stream.readToken().getWord(), "A");
    EXPECT_EQ(stream.readToken().getWord(), "B");
    EXPECT_TRUE(stream.readToken().isEOF());
    stream.putBack(b);
    stream.start();
    EXPECT_EQ(stream.readToken().getWord(), "A");
}

TEST(ParserStreamTest, OwnsInputAndRestoresStackOnException) {
    int destroyed = 0;
    ThrowingParser parser;
    EXPECT_THROW(parser.run(new CountingStream(destroyed)), std::runtime_error);
    EXPECT_EQ(destroyed, 1);
    EXPECT_TRUE(OpalParser::readToken().isEOF());
    OpalParser normal;
    normal.run(new CountingStream(destroyed));
    EXPECT_EQ(destroyed, 2);
    EXPECT_TRUE(OpalParser::readToken().isEOF());
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ippl::initialize(argc, argv);
    const int result = RUN_ALL_TESTS();
    OpalData::deleteInstance();
    ippl::finalize();
    return result;
}
