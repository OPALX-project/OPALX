#include "Utilities/DomainError.h"
#include "Utilities/GeneralOpalException.h"
#include "Utilities/OpalException.h"
#include "Utilities/ParseError.h"
#include "Utilities/SDDSParser/SDDSParserException.h"

#include <gtest/gtest.h>

TEST(OpalExceptionTest, StoresLocationAndMessage) {
    OpalException ex("method", "message");

    EXPECT_EQ(ex.where(), "method");
    EXPECT_EQ(ex.what(), "message");
}

TEST(OpalExceptionTest, CopyPreservesLocationAndMessage) {
    OpalException ex("method", "message");
    OpalException copy(ex);

    EXPECT_EQ(copy.where(), "method");
    EXPECT_EQ(copy.what(), "message");
}

TEST(OpalExceptionTest, GenericExceptionIsCaughtAsOpalException) {
    try {
        throw GeneralOpalException("generic", "failure");
    } catch (const OpalException& ex) {
        EXPECT_EQ(ex.where(), "generic");
        EXPECT_EQ(ex.what(), "failure");
        return;
    }

    FAIL() << "GeneralOpalException was not caught as OpalException";
}

TEST(OpalExceptionTest, ParserExceptionsAreCaughtAsOpalException) {
    EXPECT_THROW(throw ParseError("parse", "bad token"), OpalException);
    EXPECT_THROW(throw SDDSParserException("sdds", "bad column"), OpalException);
}

TEST(OpalExceptionTest, ArithmeticDerivedExceptionsAreCaughtAsOpalException) {
    EXPECT_THROW(throw DomainError("domain"), OpalException);
}
