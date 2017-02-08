#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

#include "Test1.h"
#include "Algorithms/StatisticalErrors.h"
#include "Structure/OpalInputInterpreter.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>


TEST(StatisticalErrorTests, CPPCommentTest) {
    std::istringstream in(inputCPPCommentTest);
    OpalInputInterpreter interpreter;

    interpreter.parse(in, interpreter.sourceParts_m.end());

    ast::variant_t front = interpreter.sourceParts_m.front();
    EXPECT_EQ(ast::TEXT, boost::get<int>(front));

    interpreter.sourceParts_m.pop_front();

    front = interpreter.sourceParts_m.front();
    EXPECT_EQ(outputCPPCommentTest, boost::get<std::string>(front));
}


TEST(StatisticalErrorTests, CCommentTest) {
    std::istringstream in(inputCCommentTest);
    OpalInputInterpreter interpreter;

    interpreter.parse(in, interpreter.sourceParts_m.end());

    ast::variant_t front = interpreter.sourceParts_m.front();
    EXPECT_EQ(ast::TEXT, boost::get<int>(front));

    interpreter.sourceParts_m.pop_front();

    front = interpreter.sourceParts_m.front();
    EXPECT_EQ(outputCCommentTest, boost::get<std::string>(front));
}

TEST(StatisticalErrorTests, TGaussTest) {
    std::istringstream in(inputTGaussTest);
    OpalInputInterpreter interpreter;

    interpreter.parse(in, interpreter.sourceParts_m.end());

    ast::variant_t front = interpreter.sourceParts_m.front();
    EXPECT_EQ(ast::TEXT, boost::get<int>(front));

    interpreter.sourceParts_m.pop_front();

    front = interpreter.sourceParts_m.front();
    EXPECT_EQ(outputTGaussTest1, boost::get<std::string>(front));

    interpreter.sourceParts_m.pop_front();

    front = interpreter.sourceParts_m.front();
    EXPECT_EQ(ast::TGAUSS, boost::get<int>(front));

    interpreter.sourceParts_m.pop_front();

    front = interpreter.sourceParts_m.front();
    EXPECT_EQ(-3.141e-3, boost::get<double>(front));

    ast::variant_t back = interpreter.sourceParts_m.back();
    EXPECT_EQ(outputTGaussTest2, boost::get<std::string>(back));
}

TEST(StatisticalErrorTests, GaussTest) {
    std::istringstream in(inputGaussTest);
    OpalInputInterpreter interpreter;

    interpreter.parse(in, interpreter.sourceParts_m.end());

    ast::variant_t front = interpreter.sourceParts_m.front();
    EXPECT_EQ(ast::TEXT, boost::get<int>(front));

    interpreter.sourceParts_m.pop_front();

    front = interpreter.sourceParts_m.front();
    EXPECT_EQ(outputGaussTest1, boost::get<std::string>(front));

    interpreter.sourceParts_m.pop_front();

    front = interpreter.sourceParts_m.front();
    EXPECT_EQ(ast::GAUSS, boost::get<int>(front));

    interpreter.sourceParts_m.pop_front();

    front = interpreter.sourceParts_m.front();
    EXPECT_EQ(ast::TEXT, boost::get<int>(front));

    interpreter.sourceParts_m.pop_front();

    front = interpreter.sourceParts_m.front();
    EXPECT_EQ(outputGaussTest2, boost::get<std::string>(front));

    interpreter.sourceParts_m.pop_front();

    front = interpreter.sourceParts_m.front();
    EXPECT_EQ(ast::TGAUSS, boost::get<int>(front));

    interpreter.sourceParts_m.pop_front(); // -3.141e-3
    interpreter.sourceParts_m.pop_front(); // ast::TEXT
    interpreter.sourceParts_m.pop_front();

    front = interpreter.sourceParts_m.front();
    EXPECT_EQ(outputTGaussTest2, boost::get<std::string>(front));
}

TEST(StatisticalErrorTests, CallTest) {
    std::istringstream in(inputCallTest);
    OpalInputInterpreter interpreter;

    interpreter.parse(in, interpreter.sourceParts_m.end());

    ast::variant_t front = interpreter.sourceParts_m.front();
    EXPECT_EQ(ast::TEXT, boost::get<int>(front));

    interpreter.sourceParts_m.pop_front();

    front = interpreter.sourceParts_m.front();
    EXPECT_EQ(outputCallTest1, boost::get<std::string>(front));

    interpreter.sourceParts_m.pop_front();

    front = interpreter.sourceParts_m.front();
    EXPECT_EQ(ast::CALL, boost::get<int>(front));

    interpreter.sourceParts_m.pop_front();

    front = interpreter.sourceParts_m.front();
    EXPECT_EQ("path/to/some/OPAL_input-file.ver3.141.txt", boost::get<std::string>(front));
}

TEST(StatisticalErrorTests, FileTest) {
    std::ofstream ofh("FileTest.in");
    ofh << inputGaussTest;
    ofh.close();

    OpalInputInterpreter interpreter("FileTest.in");

    auto it = interpreter.sourceParts_m.begin();
    std::advance(it, 4);

    EXPECT_EQ(outputGaussTest2, boost::get<std::string>(*it));

    boost::filesystem::remove("FileTest.in");
}

TEST(StatisticalErrorTests, SubFileTest) {
    std::ofstream ofh("FileTest.in");
    ofh << inputSubFileTest1;
    ofh.close();
    ofh.open("SubFileTest.in");
    ofh << inputSubFileTest2;
    ofh.close();

    OpalInputInterpreter interpreter("FileTest.in");

    auto it = interpreter.sourceParts_m.begin();
    std::advance(it, 4);

    EXPECT_EQ(ast::GAUSS, boost::get<int>(*it));

    std::advance(it, 4);
    EXPECT_EQ(3.1e-4, boost::get<double>(*it));

    std::advance(it, 4);
    EXPECT_EQ(outputSubFileTest1, boost::get<std::string>(*it));

    boost::filesystem::remove("SubFileTest.in");
    boost::filesystem::remove("FileTest.in");
}

TEST(StatisticalErrorTests, ASTTest) {
    std::ofstream ofh("FileTest.in");
    ofh << inputSubFileTest1;
    ofh.close();
    ofh.open("SubFileTest.in");
    ofh << inputSubFileTest2;
    ofh.close();

    OpalInputInterpreter interpreter("FileTest.in");
    std::string output = interpreter.processASTReference();

    EXPECT_EQ(outputASTTest + originalTrack, output);

    boost::filesystem::remove("SubFileTest.in");
    boost::filesystem::remove("FileTest.in");
}

TEST(StatisticalErrorTests, ReplaceTest) {
    std::ofstream ofh("FileTest.in");
    ofh << inputSubFileTest1;
    ofh.close();
    ofh.open("SubFileTest.in");
    ofh << inputSubFileTest2;
    ofh.close();

    OpalInputInterpreter interpreter("FileTest.in");
    interpreter.replaceString("STATISTICAL-ERRORS\\((.*?),\\s*\\d*\\s*,\\s*\\d*\\s*\\)",
                              "${1}");
    std::string output = interpreter.processASTReference();

    EXPECT_EQ(outputASTTest + replacedTrack, output);

    boost::filesystem::remove("SubFileTest.in");
    boost::filesystem::remove("FileTest.in");
}