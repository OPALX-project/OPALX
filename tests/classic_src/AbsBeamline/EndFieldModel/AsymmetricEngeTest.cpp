#include <math.h>

#include "gtest/gtest.h"
#include "Classic/AbsBeamline/EndFieldModel/AsymmetricEnge.h"

double engeTest(double delta_s, double c, double lambda_end) {
    double value = 1/(1+exp(c*delta_s/lambda_end));
    return value;
}

double tryEnge(double s, double c1, double c2, double l1, double l2, double x1, double x2) {
    double out = -1.0;
    out += engeTest(s-x1, c1, l1);
    out += engeTest(-x2-s, c2, l2);
    return out;
}

TEST(AsymmetricEngeTest, FunctionTest) {
    endfieldmodel::AsymmetricEnge enge({0.0, 1.0}, 10.0, 3.0, 
                                       {0.0, 4.0}, 11.0, 6.0);
    for (double s = -20.0; s < 20.0; s += 1.0) {
        double test = enge.function(s, 0);
        double ref = tryEnge(s, 1, 4, 3, 6, 10, 11);
        EXPECT_NEAR(test, ref, 1e-12) << s << " " << test << " " << ref;
    }

}
