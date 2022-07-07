#include "gtest/gtest.h"
#include "Classic/AbsBeamline/EndFieldModel/AsymmetricEnge.h"

TEST(AsymmetricEngeTest, FunctionTest) {
    endfieldmodel::AsymmetricEnge enge({1.0}, 10.0, 3.0, 
                                       {4.0}, 11.0, 6.0);
}
