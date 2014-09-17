/* 
 *  Copyright (c) 2012, Chris Rogers
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

#include <exception>

#include "gtest/gtest.h"
#include "Utilities/OpalException.h"
#include "Structure/BoundaryGeometry.h"
#include "AbsBeamline/SBend3D.h"

SBend3D* loadFieldMap() {
    try {
        SBend3D* sbend3D = new SBend3D("name");
        sbend3D->setLengthUnits(10.); // cm
        sbend3D->setFieldUnits(1.e-4); // ?
        sbend3D->setFieldMapFileName("src/unit_tests/classic_src/AbsBeamline/SBend3D_map.field");
        return sbend3D;
    } catch (std::exception& exc) {
        return NULL;
    }
}

TEST(SBend3DTest, SBend3DGeometryTest) {
    SBend3D* field = loadFieldMap();
    if (field == NULL)
        return; // skip the test
    Vector_t B, E, centroid;
    double radius = 2110.;
    for (double phi = -Physics::pi/4.+1e-3; phi < Physics::pi/2.; phi += Physics::pi/20)
        for (double r = -9.9; r < 31.; r += 5.)
            for (double y = -6.1; y < 6.; y += 1.) {
                Vector_t B(0., 0., 0.);
                Vector_t pos(
                    (radius+r)*cos(phi)-radius,
                    y,
                    (radius+r)*sin(phi)
                );
                field->apply(pos, centroid, 0, E, B);
                if (r > -10. && r < 10. &&
                    phi > 0. && phi < Physics::pi/4. &&
                    y > -5. && y < 5.) {
                    EXPECT_GT(fabs(B(1)), 1.e-4) // 1 Gauss
                          << "Pol:  " << r+radius << " " << y << " " << phi
                          << " Pos: " << pos(0) << " " << pos(1) << " " << pos(2)
                          << "   B: " << B(0) << " " << B(1) << " " << B(2)
                          << std::endl;
                } else {
                    EXPECT_LT(fabs(B(1)), 1.)
                          << "Pol:  " << r+radius << " " << y << " " << phi
                          << " Pos: " << pos(0) << " " << pos(1) << " " << pos(2)
                          << "   B: " << B(0) << " " << B(1) << " " << B(2)
                          << std::endl;
                }
    }
}

void testField(double r, double y, double phi,
               double bx, double by, double bz,
               double tol) {
    SBend3D* field = loadFieldMap();
    if (field == NULL) {
        std::cerr << "SKIPPING SBEND3D TEST - FAILED TO LOAD FIELD" << std::endl;
        return; // skip the test
    }
    Vector_t B, E, centroid, pos;
    double radius = 2110.;
    pos = Vector_t(
        (radius+r)*cos(phi)-radius,
        y,
        (radius+r)*sin(phi)
    );
    field->apply(pos, centroid, 0, E, B);
    // the field map is rotated through pi/8. (into start at z=0. coordinates)
    double sR = sin(Physics::pi/8.);
    double cR = cos(Physics::pi/8.);
    double bxR = bx*cR-bz*sR;
    double bzR = bz*cR+bx*sR;
    EXPECT_NEAR(B(0), bxR, tol) << "R_p: " << r << " " << y << " " << phi
                                << " B: " << bx << " " << by << " " << bz;
    EXPECT_NEAR(B(1), by, tol) << "R_p: " << r << " " << y << " " << phi
                               << " B: " << bx << " " << by << " " << bz;
    EXPECT_NEAR(B(2), bzR, tol) << "R_p: " << r << " " << y << " " << phi
                                << " B: " << bx << " " << by << " " << bz;
}

TEST(SBend3DTest, SBend3DFieldTest) {
    // 211.000000000       0.00000000000       0.00000000000       0.00000000000              0.00000000000    
    testField(0., 0., Physics::pi/8., 0., 5346.80707376*1e-4, 0., 1e-6);

    // 211.000000000      0.250000000000       0.00000000000       71.6525328925       5351.82535885     -0.156196844700E-02
    testField(0., 2.5, Physics::pi/8.,
               71.65253*1e-4, 5351.8253*1e-4, -0.1561E-02*1e-4, 1e-7);
}

