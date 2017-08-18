/*
 *  Copyright (c) 2017, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include <cmath>
#include <fstream>
#include <sstream>

#include "gtest/gtest.h"

#include "Classic/AbsBeamline/EndFieldModel/Tanh.h"
#include "Classic/AbsBeamline/ScalingFFAGMagnet.h"
#include "Classic/AbsBeamline/Offset.h"

class ScalingFFAGMagnetTest : public ::testing::Test { 
public: 
    ScalingFFAGMagnetTest() : sector_m(NULL), fout_m() { 
    }
 
    void SetUp( ) { 
        sector_m = new ScalingFFAGMagnet("test");
        // characteristic length is R*dphi => 0.6545 m
        endfieldmodel::Tanh* tanh = new endfieldmodel::Tanh(psi0_m, psi0_m/5., 20);
        sector_m->setEndField(tanh);
        sector_m->setTanDelta(tan(M_PI/4.));
        sector_m->setR0(r0_m);
        sector_m->setRMin(0.);
        sector_m->setRMax(r0_m*2.);
        sector_m->setDipoleConstant(1.);
        sector_m->setFieldIndex(15);
        sector_m->setMaxOrder(10);
        sector_m->setPhiStart(psi0_m);
        sector_m->setPhiEnd(psi0_m*4.);
        sector_m->setAzimuthalExtent(M_PI);
        sector_m->setVerticalExtent(1.); // 1 m
        sector_m->initialise();
    }
 
    void TearDown( ) { 
        delete sector_m;
        sector_m = NULL;
    }

    ~ScalingFFAGMagnetTest() {
    }

    ScalingFFAGMagnet* sector_m;
    std::ofstream fout_m;
    double r0_m = 24; // m
    double magnetLength_m = 0.63; // m
    double psi0_m = magnetLength_m*2*M_PI/r0_m; // radians

    double getDBDu(int i, int j, Vector_t pos, double delta, bool isCartesian) {
        Vector_t mom(0., 0., 0.);
        Vector_t E(0., 0., 0.);
        Vector_t BUpper(0., 0., 0.);
        double t = 0;
        double derivative = 0;
        Vector_t posUpper = pos;
        posUpper[j] += delta;
        if (isCartesian) {
            sector_m->apply(posUpper, mom, t, E, BUpper);
        } else {
            sector_m->getFieldValueCylindrical(posUpper, BUpper);
        }
        Vector_t BLower(0., 0., 0.);
        Vector_t posLower = pos;
        posLower[j] -= delta;
        if (isCartesian) {
            sector_m->apply(posLower, mom, t, E, BLower);
        } else {
            sector_m->getFieldValueCylindrical(posLower, BLower);
        }
        derivative = (BUpper[i]-BLower[i])/(posUpper[j]-posLower[j]);
        return derivative;
    }

    double getDivBCart(Vector_t pos, Vector_t delta) {
        Vector_t div_elements(0., 0., 0.);
        // dB_x/dx ~ (B_x(x+dx) - B_x(x-dx))/2/dx
        for (size_t i = 0; i < 3; ++i) {
            div_elements[i] = getDBDu(i, i, pos, delta[i], true);
        }
        return div_elements[0] + div_elements[1] + div_elements[2];
    }

    Vector_t getCurlBCart(Vector_t posCart, Vector_t delta) {
        Vector_t curlElements(0., 0., 0.);
        // dB_x/dx ~ (B_x(x+dx) - B_x(x-dx))/2/dx
        curlElements[0] += getDBDu(2, 1, posCart, delta[1], true); // dBz/dy
        curlElements[0] -= getDBDu(1, 2, posCart, delta[2], true); // dBy/dz
        curlElements[1] += getDBDu(0, 2, posCart, delta[2], true); // dBx/dz
        curlElements[1] -= getDBDu(2, 0, posCart, delta[0], true); // dBz/dx
        curlElements[2] += getDBDu(1, 0, posCart, delta[0], true); // dBy/dx
        curlElements[2] -= getDBDu(0, 1, posCart, delta[1], true); // dBx/dy
        return curlElements;
    }

    Vector_t getCurlBCyl(Vector_t posCyl, Vector_t delta) {
        // posCyl goes like r, z, phi
        Vector_t curlElements(0., 0., 0.);
        Vector_t B(0., 0., 0.);
        sector_m->getFieldValueCylindrical(posCyl, B);

        curlElements[0] += getDBDu(2, 1, posCyl, delta[1], false); // dBphi/dz
        curlElements[0] -= getDBDu(1, 2, posCyl, delta[2], false)/posCyl[0]; // 1/r*dBz/dphi

        curlElements[1] += getDBDu(2, 0, posCyl, delta[0], false); // dBphi/dr
        curlElements[1] -= getDBDu(0, 2, posCyl, delta[2], false)/posCyl[0]; // dBr/dphi/r
        curlElements[1] += B[2]/posCyl[0]; // Bphi/r

        curlElements[2] += getDBDu(0, 1, posCyl, delta[1], false); // dBr/dz
        curlElements[2] -= getDBDu(1, 0, posCyl, delta[0], false); // dBz/dr

        return curlElements;
    }

    double getDivBCylindrical(Vector_t posCyl, Vector_t delta) {
        Vector_t B(0., 0., 0.);
        Vector_t div_elements(0., 0., 0.);
        sector_m->getFieldValueCylindrical(posCyl, B);
        div_elements[0] = B[0]/posCyl[0] + getDBDu(0, 0, posCyl, delta[0], false);
        div_elements[1] = getDBDu(1, 1, posCyl, delta[1], false);
        div_elements[2] = getDBDu(2, 2, posCyl, delta[2], false)/posCyl[0];
        // dB_x/dx ~ (B_x(x+dx) - B_x(x-dx))/2/dx
        return div_elements[0] + div_elements[1] + div_elements[2];
    }

    void printCoefficients() {
        std::vector<std::vector<double> > coeff = sector_m->getDfCoefficients();
        std::cerr << "Coefficients" << std::endl;
        for (size_t n = 0; n < coeff.size(); ++n) {
            for (size_t i = 0; i < coeff[n].size(); ++i) {
                std::cerr << coeff[n][i] << " ";
            }
            std::cerr << std::endl;
        }
    }

    bool printLine(Vector_t posCyl, double aux, std::ofstream& fout, double maxwell_tolerance) {
        double r = posCyl[0];
        double y = posCyl[1];
        double phi = posCyl[2];
        double x = r*sin(phi);
        double z = r*cos(phi);
        Vector_t posCart(x, y, z);
        Vector_t mom(0., 0., 0.);
        Vector_t E(0., 0., 0.);
        Vector_t B(0., 0., 0.);
        Vector_t BCart(0., 0., 0.);
        double t = 0;
        
        sector_m->getFieldValueCylindrical(posCyl, B);
        bool outsideFieldMap = sector_m->apply(posCart, mom, t, E, BCart);
        double delta = y/1000.;
        if (delta > 1e-3 or delta < 1e-9) {
            delta = 1e-3;
        }
        double divBCyl =  getDivBCylindrical(posCyl, Vector_t(delta, delta, delta/r0_m));
        Vector_t curlBCyl = getCurlBCyl(posCyl, Vector_t(delta, delta, delta/r0_m)); //curlBCyl[0] should be cancelled by 2nd term
        double divBCart = getDivBCart(posCart, Vector_t(delta, delta, delta));
        Vector_t curlBCart = getCurlBCart(posCart, Vector_t(delta, delta, delta)); //curlBCyl[0] should be cancelled by 2nd term
        std::stringstream ssout;
        ssout << "order " << sector_m->getMaxOrder()
             << " (x,y,z)  " << x << "    " << y << "    " << z
             << " (r,phi) " << r << "    " << phi
             << " B(r,z,phi) " << B[0] << "    " << B[1] << "    " << B[2]
             << " DivB_Cyl: " << divBCyl
             << " CurlB_Cyl: " << curlBCyl
             << " B(x,y,z) " << BCart[0] << "    " << BCart[1] << "    " << BCart[2]
             << " DivB_Cart: " << divBCart
             << " CurlB_Cart: " << curlBCart
             << " Aux: " << aux;
        fout << ssout.str() << std::endl;
        bool passtest = !outsideFieldMap;
        passtest &= divBCart < maxwell_tolerance;
        // for (size_t i = 0; i < 3; ++i) {
        //     passtest &= curlBCart[i] < maxwell_tolerance;
        // }
        if (!passtest) {
            std::cout << ssout.str() << std::endl;
        }
        return passtest;
    }

private:

};

TEST_F(ScalingFFAGMagnetTest, ConstructorTest) {
    ScalingFFAGMagnet* test = new ScalingFFAGMagnet("test");
    std::vector<int> data(15);
    size_t i = 0;
    test->setTanDelta(++i);
    test->setFieldIndex(++i);
    test->setDipoleConstant(++i);
    test->setR0(++i);
    double x = ++i;
    test->setCentre(Vector_t(x, x, x));
    x = ++i;
    endfieldmodel::Tanh* tanh = new endfieldmodel::Tanh(x, x, i);
    test->setEndField(tanh);
    test->setMaxOrder(++i);
    test->setPhiStart(++i);
    test->setPhiEnd(++i);
    test->setRMin(++i);
    test->setRMax(++i);
    test->setAzimuthalExtent(++i);

    std::vector<ScalingFFAGMagnet*> magnets(2);
    magnets[0] = test;
    magnets[1] = dynamic_cast<ScalingFFAGMagnet*>(test->clone());
    for (size_t j = 0; j < magnets.size(); ++j) {
        i = 0;
        test = magnets[j];
        EXPECT_NEAR(test->getTanDelta(), ++i, 1e-9);
        EXPECT_EQ(test->getFieldIndex(), ++i);
        EXPECT_NEAR(test->getDipoleConstant(), ++i, 1e-9);
        EXPECT_NEAR(test->getR0(), ++i, 1e-9);
        EXPECT_NEAR(test->getCentre()[0], ++i, 1e-9);
        EXPECT_NEAR(test->getCentre()[1], i, 1e-9);
        EXPECT_NEAR(test->getCentre()[2], i, 1e-9);
        ++i;
        EXPECT_NEAR(test->getEndField()->function(x, 0),
                    tanh->function(x, 0), 1e-9);
        EXPECT_EQ(test->getMaxOrder(), ++i);
        EXPECT_NEAR(test->getPhiStart(), ++i, 1e-9);
        EXPECT_NEAR(test->getPhiEnd(), ++i, 1e-9);
        EXPECT_NEAR(test->getRMin(), ++i, 1e-9);
        EXPECT_NEAR(test->getRMax(), ++i, 1e-9);
        EXPECT_NEAR(test->getAzimuthalExtent(), ++i, 1e-9);
    }

    endfieldmodel::Tanh* tanh2 = new endfieldmodel::Tanh(-10., -10., 10);
    magnets[0]->setEndField(tanh2);
    EXPECT_NEAR(magnets[0]->getEndField()->function(10., 0),
                tanh2->function(10., 0), 1e-9);
    delete magnets[0];
    delete magnets[1];
}

TEST_F(ScalingFFAGMagnetTest, PlacementTest) {
    double phiTest[] = {0., psi0_m, psi0_m*2.};
    double byTest[] = {0.5, 1., 0.5};
    for (size_t i = 0; i < 3; ++i) {
        Vector_t mom, E, B;
        double t = 0;
        Vector_t posCart(r0_m*sin(phiTest[i]), 0., r0_m*cos(phiTest[i]));
        sector_m->apply(posCart, mom, t, E, B);
        EXPECT_NEAR(B[1], byTest[i], 1e-3);
    }
}

TEST_F(ScalingFFAGMagnetTest, BTwoDTest) {
    std::ofstream fout("/tmp/b_twod.out");
    bool passtest = true;
    for (double y = 0.; y < 0.025; y += 0.015) {
        for (double r = r0_m; r < r0_m+1; r += 0.01) {
            for (double psi = -2.*psi0_m; psi < 4.00001*psi0_m; psi += psi0_m/10.) {
                passtest &= printLine(Vector_t(r, y, psi), 0., fout, 1e-1);
            }
        }
    }
    EXPECT_TRUE(passtest);
}

TEST_F(ScalingFFAGMagnetTest, ConvergenceYTest) {
    std::ofstream fout("/tmp/convergence_y.out");
    bool passtest = true;
    for (double y = 0.00001; y < 0.02; y *= 2.) {
        passtest &= printLine(Vector_t(r0_m, y, 2.*psi0_m), 0., fout, 1e-3);
    }
    for (double y = 0.051; y < 0.1; y += 0.002) {
        passtest &= printLine(Vector_t(r0_m, y, 2.*psi0_m), 0., fout, 1e-3);
    }
    EXPECT_TRUE(passtest);
}

TEST_F(ScalingFFAGMagnetTest, ConvergenceOrderTest) {
    for (double y = 0.05; y > 0.0002; y /= 10.) {
        std::vector<double> divBVec(3);
        double delta = y/100.;
        for (size_t i = 0; i < divBVec.size(); ++i) {
            sector_m->setMaxOrder((i+1)*2);
            sector_m->initialise();
            Vector_t pos(r0_m, y, psi0_m*2);
            double divB = getDivBCylindrical(pos, Vector_t(delta, delta, delta/r0_m));
            divB = fabs(divB);
            divBVec[i] = divB;
            if (i > 0) {
                EXPECT_LT(divBVec[i], divBVec[i-1]) << " with i "
                                                    << i << std::endl;
            }
        }
    }
    sector_m->setMaxOrder(10);
}

TEST_F(ScalingFFAGMagnetTest, ConvergenceEndLengthTest) {
    std::ofstream fout("/tmp/convergence_endlength.out");
    bool passtest = true;
    for (double endLength = 1.; endLength < 10.1; endLength += 1.) {
        endfieldmodel::Tanh* tanh = new endfieldmodel::Tanh(psi0_m, psi0_m/endLength, 20);
        sector_m->setEndField(tanh);
        sector_m->initialise();
        passtest &= printLine(Vector_t(r0_m, 0.05, 2.*psi0_m), psi0_m/endLength*r0_m, fout, 1e-3);
    }
    EXPECT_TRUE(passtest);
}

TEST_F(ScalingFFAGMagnetTest, VerticalBoundingBoxTest) {
    sector_m->setVerticalExtent(0.1);
    Vector_t mom, E, B;
    double t = 0;
    Vector_t pos(r0_m*sin(psi0_m), 0.09, r0_m*cos(psi0_m));

    EXPECT_FALSE(sector_m->apply(pos, mom, t, E, B));
    pos[1] = 0.11;
    EXPECT_TRUE(sector_m->apply(pos, mom, t, E, B));
    pos[1] = -0.11;
    EXPECT_TRUE(sector_m->apply(pos, mom, t, E, B));
    pos[1] = -0.09;
    EXPECT_FALSE(sector_m->apply(pos, mom, t, E, B));
}

TEST_F(ScalingFFAGMagnetTest, RadialBoundingBoxTest) {
    sector_m->setRMin(r0_m-0.1);
    Vector_t mom, E, B;
    double t = 0;
    double r1 = r0_m-0.09;
    Vector_t pos1(r1*sin(psi0_m), 0.0, r1*cos(psi0_m));
    double r2 = r0_m-0.11;
    Vector_t pos2(r2*sin(psi0_m), 0.0, r2*cos(psi0_m));
    EXPECT_FALSE(sector_m->apply(pos1, mom, t, E, B));
    EXPECT_TRUE(sector_m->apply(pos2, mom, t, E, B));

    sector_m->setRMax(r0_m+0.1);
    double r3 = r0_m+0.09;
    Vector_t pos3(r3*sin(psi0_m), 0.0, r3*cos(psi0_m));
    double r4 = r0_m+0.11;
    Vector_t pos4(r4*sin(psi0_m), 0.0, r4*cos(psi0_m));
    EXPECT_FALSE(sector_m->apply(pos3, mom, t, E, B));
    EXPECT_TRUE(sector_m->apply(pos4, mom, t, E, B));
}

TEST_F(ScalingFFAGMagnetTest, AzimuthalBoundingBoxTest) {
    sector_m->setAzimuthalExtent(psi0_m*5.);
    sector_m->setPhiStart(psi0_m*3.);
    Vector_t mom, E, B;
    double t = 0;
    double phi[] = {-2.1*psi0_m, -1.9*psi0_m, 7.9*psi0_m, 8.1*psi0_m};
    bool bb[] = {true, false, false, true};
    for(size_t i = 0; i < 4; ++i) {
        Vector_t pos(r0_m*sin(phi[i]), 0.0, r0_m*cos(phi[i]));
        EXPECT_EQ(sector_m->apply(pos, mom, t, E, B), bb[i]) << i << " " << pos;
    }
}


TEST_F(ScalingFFAGMagnetTest, GeometryTest) {
    Euclid3D delta = sector_m->getGeometry().getTotalTransform();
    Vector3D vec = delta.getVector();
    Vector3D rot = delta.getRotation().getAxis();
    EXPECT_EQ(vec(0), r0_m*(cos(psi0_m*4.)-1));
    EXPECT_EQ(vec(1), 0.);
    EXPECT_EQ(vec(2), r0_m*sin(psi0_m*4.));

    EXPECT_EQ(rot(0), 0.);
    EXPECT_EQ(rot(1), -psi0_m*4);
    EXPECT_EQ(rot(2), 0.);
}

