#include "gtest/gtest.h"

#include "Index/SOffset.h"

constexpr unsigned Dim = 2;
constexpr double   roundOffError = 1e-10;

TEST(Index, SOffset)
{
    // testing SOffset creation
    SOffset<Dim> A(1, 2);
    SOffset<Dim> B(1, 1);
    EXPECT_NEAR(A[0],1, roundOffError);
    EXPECT_NEAR(A[1],2, roundOffError);
    EXPECT_NEAR(B[0],1, roundOffError);
    EXPECT_NEAR(B[1],1, roundOffError);

    // testing SOffset [] operator
    SOffset<Dim> soLeft;
    SOffset<Dim> soRight;
    for (unsigned int d=0; d < Dim; d++) {
        soLeft[d] = 0;
        soRight[d] = (d == 1 ? 1 : 0);
    }
    EXPECT_NEAR(soLeft[0], 0, roundOffError);
    EXPECT_NEAR(soLeft[1], 0, roundOffError);
    EXPECT_NEAR(soRight[0], 0, roundOffError);
    EXPECT_NEAR(soRight[1], 1, roundOffError);

    // testing SOffset addition, subtraction, and copy constructor
    SOffset<Dim> C(A + B);
    SOffset<Dim> D(A - B);
    EXPECT_NEAR(C[0], 2, roundOffError);
    EXPECT_NEAR(C[1], 3, roundOffError);
    EXPECT_NEAR(D[0], 0, roundOffError);
    EXPECT_NEAR(D[1], 1, roundOffError);

    // testing SOffset +=, -= operators
    D += C;
    EXPECT_NEAR(D[0], 2, roundOffError);
    EXPECT_NEAR(D[1], 4, roundOffError);
    C -= D;
    EXPECT_NEAR(C[0],  0, roundOffError);
    EXPECT_NEAR(C[1], -1, roundOffError);

    // testing SOffset comparisons
    EXPECT_TRUE(B < A);
    EXPECT_TRUE(B <= A);
    EXPECT_TRUE(D > C);
    EXPECT_TRUE(D >= C);
    EXPECT_TRUE(A != B);
    EXPECT_TRUE(B == B);

    // testing containment check
    NDIndex<Dim> N1(Index(0,1), Index(0,1));
    EXPECT_TRUE(! A.inside(N1));
    EXPECT_TRUE(  B.inside(N1));
}