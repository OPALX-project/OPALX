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
        sector_m->setPhiEnd(M_PI/12.); // 15 degrees
        sector_m->setAzimuthalExtent(M_PI);
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
    double psi0_m = magnetLength_m/r0_m/2.; // radians

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

    void printLine(Vector_t posCyl, double aux, std::ofstream& fout) {
        double r = posCyl[0];
        double y = posCyl[1];
        double phi = posCyl[2];
        double x = r*cos(phi);
        double z = r*sin(phi);
        Vector_t posCart(x, y, z);
        Vector_t mom(0., 0., 0.);
        Vector_t E(0., 0., 0.);
        Vector_t B(0., 0., 0.);
        Vector_t BCart(0., 0., 0.);
        double t = 0;
        
        sector_m->getFieldValueCylindrical(posCyl, B);
        sector_m->apply(posCart, mom, t, E, BCart);
        double delta = y/1000.;
        if (delta > 1e-3 or delta < 1e-9) {
            delta = 1e-3;
        }
        double divBCyl =  getDivBCylindrical(posCyl, Vector_t(delta, delta, delta/r0_m));
        Vector_t curlBCyl = getCurlBCyl(posCyl, Vector_t(delta, delta, delta/r0_m)); //curlBCyl[0] should be cancelled by 2nd term
        double divBCart = getDivBCart(posCart, Vector_t(delta, delta, delta));
        Vector_t curlBCart = getCurlBCart(posCart, Vector_t(delta, delta, delta)); //curlBCyl[0] should be cancelled by 2nd term
        fout << "order " << sector_m->getMaxOrder()
             << " (x,y,z)  " << x << "    " << y << "    " << z
             << " (r,phi) " << r << "    " << phi
             << " B(r,z,phi) " << B[0] << "    " << B[1] << "    " << B[2]
             << " DivB_Cyl: " << divBCyl
             << " CurlB_Cyl: " << curlBCyl
             << " B(x,y,z) " << BCart[0] << "    " << BCart[1] << "    " << BCart[2]
             << " DivB_Cart: " << divBCart
             << " CurlB_Cart: " << curlBCart
             << " Aux: " << aux
             << std::endl;
    }

private:

};

TEST_F(ScalingFFAGMagnetTest, BTwoDTest) {
    std::ofstream fout("/tmp/b_twod.out");
    for (double y = 0.; y < 0.021; y += 0.01) {
        for (double r = r0_m; r < r0_m+1; r += 0.005) {
            for (double psi = -5.*psi0_m; psi < 5.00001*psi0_m; psi += psi0_m/20.) {
                printLine(Vector_t(r, y, psi), 0., fout);
            }
        }
    }
}

TEST_F(ScalingFFAGMagnetTest, ConvergenceYTest) {
    std::ofstream fout("/tmp/convergence_y.out");
    for (double y = 0.00001; y < 0.02; y *= 2.) {
        printLine(Vector_t(r0_m, y, psi0_m), 0., fout);
    }
    for (double y = 0.051; y < 0.1; y += 0.002) {
        printLine(Vector_t(r0_m, y, psi0_m), 0., fout);
    }
}

TEST_F(ScalingFFAGMagnetTest, ConvergenceOrderTest) {
    std::ofstream fout("/tmp/convergence_order.out");
    for (double y = 0.05; y > 0.0002; y /= 10.) {
        for (int order = 1; order < 21; ++order) {
            sector_m->setMaxOrder(order);
            sector_m->initialise();
            printLine(Vector_t(r0_m, y, psi0_m), 0., fout);
        }
    }
    sector_m->setMaxOrder(10);
}

TEST_F(ScalingFFAGMagnetTest, ConvergenceEndLengthTest) {
    std::ofstream fout("/tmp/convergence_endlength.out");

    for (double endLength = 1.; endLength < 10.1; endLength += 1.) {
        endfieldmodel::Tanh* tanh = new endfieldmodel::Tanh(psi0_m, psi0_m/endLength, 20);
        sector_m->setEndField(tanh);
        sector_m->initialise();
        printLine(Vector_t(r0_m, 0.05, psi0_m), psi0_m/endLength*r0_m, fout);
    }
}

TEST_F(ScalingFFAGMagnetTest, GeometryTest) {
    // Sucked geometry information from
    //     Classic/AbsBeamline/Ring.cpp::appendElement
    // Transform in OPAL-T coords
    Offset refOffset = Offset::localCylindricalOffset("test", M_PI/24., M_PI/24., 24.*M_PI/12.);
    Euclid3D refDelta = refOffset.getGeometry().getTotalTransform();   
    Vector3D refV = refDelta.getVector();
    Vector3D refR = refDelta.getRotation().getAxis();
    Euclid3D refEuclid(refV(2), refV(0), -refV(1), refR(2), refR(0), -refR(1));
    refDelta = refEuclid;
    Vector_t refDeltaPos(refDelta.getVector()(0), refDelta.getVector()(1), 0);
    double refEndRot = refDelta.getRotation().getAxis()(2);
    Vector_t refDeltaNorm(cos(refEndRot), sin(refEndRot), 0.);

    Euclid3D delta = sector_m->getGeometry().getTotalTransform();
    Vector3D v = delta.getVector();
    Vector3D r = delta.getRotation().getAxis();

    // Transform to cycl coordinates
    Euclid3D euclid(v(2), v(0), -v(1), r(2), r(0), -r(1));
    delta = euclid;
    // Calculate change in position
    Vector_t deltaPos(delta.getVector()(0), delta.getVector()(1), 0);
    double endRot = delta.getRotation().getAxis()(2);
    Vector_t deltaNorm(cos(endRot), sin(endRot), 0.);

    std::cerr << 24.*(1-cos(M_PI/12.)) << " " << 24.*sin(M_PI/12.) << " ** " << cos(M_PI/12.) << " " << sin(M_PI/12.) << " ** " << M_PI/12. << std::endl;
    std::cerr << deltaPos << " ** " << deltaNorm << " ** " << endRot << std::endl;
    std::cerr << refDeltaPos << " ** " << refDeltaNorm << " ** " << refEndRot << std::endl;
    EXPECT_TRUE(false) << "I need to check that the beam and fields are aligned";
    EXPECT_TRUE(false) << "Units should be degrees not rads";
}

