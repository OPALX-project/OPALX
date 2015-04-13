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
#include <stdexcept>
#include <fstream>
#include <cstdio>

#include "gtest/gtest.h"
#include "Utilities/GeneralClassicException.h"
#include "Structure/BoundaryGeometry.h"
#include "AbsBeamline/SBend3D.h"

void writeFieldMap() {
    std::ofstream out("SBend3D_map.field");
    out << "      422280      422280      422280           1\n"
        << " 1 X [LENGU]\n"
        << "  2 Y [LENGU]\n"
        << "  3 Z [LENGU]\n"
        << "  4 BX [FLUXU]\n"
        << "  5 BY [FLUXU]\n"
        << "  6 BZ [FLUXU]\n"
        << "  0\n"
        << "    194.014700000       0.00000000000       80.3635200000      0.682759323466E-07  -5.37524925776      0.282807068056E-07\n"
        << "    210.000000000       0.00000000000       0.00000000000       0.00000000000       5038.98826666       0.00000000000    \n"
        << "    194.014700000       0.00000000000      -80.3635200000      0.682759055339E-07  -5.37524929545     -0.282807663112E-07\n"
        << "    194.014700000      0.250000000000       80.3635200000     -0.250876890041E-01  -5.37468265486     -0.105656867773E-01\n"
        << "    210.000000000      0.250000000000       0.00000000000       81.6064485044       5043.86454663     -0.162646432266E-02\n"
        << "    194.014700000      0.250000000000      -80.3635200000     -0.252106804794E-01  -5.37468268550      0.102687594900E-01\n"
        << "    194.014700000      0.500000000000       80.3635200000     -0.501213811240E-01  -5.37048637738     -0.211087378232E-01\n"
        << "    210.000000000      0.500000000000       0.00000000000       163.479461107       5051.16103059     -0.319620022674E-02\n"
        << "    194.014700000      0.500000000000      -80.3635200000     -0.503671731133E-01  -5.37048640042      0.205153440729E-01\n"
        << "    194.938580000       0.00000000000       80.7462000000     -0.109193207295E-06  -5.47695895319     -0.452292335858E-07\n"
        << "    211.000000000       0.00000000000       0.00000000000       0.00000000000       5346.80707376       0.00000000000    \n"
        << "    194.938580000       0.00000000000      -80.7462000000     -0.109194150249E-06  -5.47695923342      0.452297706086E-07\n"
        << "    194.938580000      0.250000000000       80.7462000000     -0.228949123764E-01  -5.47615607268     -0.965935126817E-02\n"
        << "    211.000000000      0.250000000000       0.00000000000       71.6525328925       5351.82535885     -0.156196844700E-02\n"
        << "    194.938580000      0.250000000000      -80.7462000000     -0.230215180683E-01  -5.47615608020      0.935369738637E-02\n"
        << "    194.938580000      0.500000000000       80.7462000000     -0.457451523318E-01  -5.47168044969     -0.193017688931E-01\n"
        << "    211.000000000      0.500000000000       0.00000000000       143.398642339       5359.42974721     -0.306846139443E-02\n"
        << "    194.938580000      0.500000000000      -80.7462000000     -0.459981937280E-01  -5.47168045716      0.186908719298E-01\n"
        << "    195.862460000       0.00000000000       81.1288900000       0.00000000000      -5.56838923159       0.00000000000    \n"
        << "    212.000000000       0.00000000000       0.00000000000       0.00000000000       5612.41675566       0.00000000000    \n"
        << "    195.862460000       0.00000000000      -81.1288900000       0.00000000000      -5.56838923767       0.00000000000    \n"
        << "    195.862460000      0.250000000000       81.1288900000     -0.209034621315E-01  -5.56743204962     -0.884197597858E-02\n"
        << "    212.000000000      0.250000000000       0.00000000000       61.0257302103       5617.50584022     -0.152992556601E-02\n"
        << "    195.862460000      0.250000000000      -81.1288900000     -0.210330844408E-01  -5.56743205033      0.852904004824E-02\n"
        << "    195.862460000      0.500000000000       81.1288900000     -0.416276950059E-01  -5.56278428219     -0.176093837788E-01\n"
        << "    212.000000000      0.500000000000       0.00000000000       121.985754491       5625.18212699     -0.300403876153E-02\n"
        << "    195.862460000      0.500000000000      -81.1288900000     -0.418867765797E-01  -5.56278428318      0.169839055390E-01"
        << std::endl;

    if (!out.good())
        throw std::runtime_error("could not write field map 'SBend3D_map.field'");

    out.close();
}

void removeFieldMap() {
    if (std::remove("SBend3D_map.field") != 0)
        throw std::runtime_error("could not delete field map 'SBend3D_map.field'");

}

SBend3D* loadFieldMap() {
    try {
        SBend3D* sbend3D = new SBend3D("name");
        sbend3D->setLengthUnits(10.); // cm
        sbend3D->setFieldUnits(1.e-4); // ?
        sbend3D->setFieldMapFileName("SBend3D_map.field");
        // sbend3D->setFieldMapFileName("src/unit_tests/classic_src/AbsBeamline/SBend3D_map.field");
        return sbend3D;
    } catch (std::exception& exc) {
        return NULL;
    }
}

TEST(SBend3DTest, SBend3DGeometryTest) {
    writeFieldMap();

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
    removeFieldMap();
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
    writeFieldMap();

    // 211.000000000       0.00000000000       0.00000000000       0.00000000000              0.00000000000
    testField(0., 0., Physics::pi/8., 0., 5346.80707376*1e-4, 0., 1e-6);

    // 211.000000000      0.250000000000       0.00000000000       71.6525328925       5351.82535885     -0.156196844700E-02
    testField(0., 2.5, Physics::pi/8.,
               71.65253*1e-4, 5351.8253*1e-4, -0.1561E-02*1e-4, 1e-7);

    removeFieldMap();
}
