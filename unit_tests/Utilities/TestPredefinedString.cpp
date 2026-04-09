#include <gtest/gtest.h>

#include "Attributes/Attributes.h"
#include "Attributes/PredefinedString.h"
#include "Distribution/Distribution.h"
#include "OpalParser/SimpleStatement.h"
#include "Utilities/ParseError.h"

using namespace Attributes;

TEST(PredefinedStringTest, TestDistributionType) {
    Distribution dist;

    auto* typeAttr = dist.findAttribute("TYPE");
    ASSERT_NE(typeAttr, nullptr);

    auto* typeAttribute = dynamic_cast<PredefinedString*>(&typeAttr->getHandler());
    ASSERT_NE(typeAttribute, nullptr);

    {
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING, "GAUSS");
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_NO_THROW(typeAttribute->parse(*typeAttr, statement, true));
        EXPECT_EQ(Attributes::getString(*typeAttr), "GAUSS");
    }

    {
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING,
                     "MULTIVARIATEGAUSS");
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_NO_THROW(typeAttribute->parse(*typeAttr, statement, true));
        EXPECT_EQ(Attributes::getString(*typeAttr), "MULTIVARIATEGAUSS");
    }

    {
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING, "gauss");
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_NO_THROW(typeAttribute->parse(*typeAttr, statement, true));
        EXPECT_EQ(Attributes::getString(*typeAttr), "GAUSS");
    }

    {
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING, "Gauss");
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_NO_THROW(typeAttribute->parse(*typeAttr, statement, true));
        EXPECT_EQ(Attributes::getString(*typeAttr), "GAUSS");
    }

    {
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING, "flattop");
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_NO_THROW(typeAttribute->parse(*typeAttr, statement, true));
        EXPECT_EQ(Attributes::getString(*typeAttr), "FLATTOP");
    }

    {
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING, "FromFile");
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_NO_THROW(typeAttribute->parse(*typeAttr, statement, true));
        EXPECT_EQ(Attributes::getString(*typeAttr), "FROMFILE");
    }

    {
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING, "guass");
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_THROW(typeAttribute->parse(*typeAttr, statement, true),
                     ParseError);
    }

    {
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING, "");
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_THROW(typeAttribute->parse(*typeAttr, statement, true), ParseError);
    }
}

