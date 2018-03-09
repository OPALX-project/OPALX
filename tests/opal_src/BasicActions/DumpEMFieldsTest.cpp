/*
 *  Copyright (c) 2014, Chris Rogers
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

#include <fstream>
#include <iostream>

#include "gtest/gtest.h"

#include "opal_src/Utilities/MockComponent.h"

#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"
#include "BasicActions/DumpEMFields.h"

#include "opal_test_utilities/SilenceTest.h"

namespace DumpEMFieldsTest {

void setOneAttribute(DumpEMFields* dump, std::string name, double value) {
    Attributes::setReal(*dump->findAttribute(name), value);
}

void setAttributesCart(DumpEMFields* dump,
                   double x0, double dx, double nx,
                   double y0, double dy, double ny,
                   double z0, double dz, double nz,
                   double t0, double dt, double nt,
                   std::string filename, bool defaultCoords = true) {
    setOneAttribute(dump, "X_START", x0);
    setOneAttribute(dump, "DX", dx);
    setOneAttribute(dump, "X_STEPS", nx);
    setOneAttribute(dump, "Y_START", y0);
    setOneAttribute(dump, "DY", dy);
    setOneAttribute(dump, "Y_STEPS", ny);
    setOneAttribute(dump, "Z_START", z0);
    setOneAttribute(dump, "DZ", dz);
    setOneAttribute(dump, "Z_STEPS", nz);
    setOneAttribute(dump, "T_START", t0);
    setOneAttribute(dump, "DT", dt);
    setOneAttribute(dump, "T_STEPS", nt);
    Attributes::setString(*dump->findAttribute("FILE_NAME"), filename);
    if (!defaultCoords) {
        Attributes::setString(*dump->findAttribute("COORDINATE_SYSTEM"), "cARtesiAN");
    }
}

void setAttributesCyl(DumpEMFields* dump,
                   double x0, double dx, double nx,
                   double y0, double dy, double ny,
                   double z0, double dz, double nz,
                   double t0, double dt, double nt,
                   std::string filename) {
    setOneAttribute(dump, "R_START", x0);
    setOneAttribute(dump, "DR", dx);
    setOneAttribute(dump, "R_STEPS", nx);
    setOneAttribute(dump, "PHI_START", y0);
    setOneAttribute(dump, "DPHI", dy);
    setOneAttribute(dump, "PHI_STEPS", ny);
    setOneAttribute(dump, "Z_START", z0);
    setOneAttribute(dump, "DZ", dz);
    setOneAttribute(dump, "Z_STEPS", nz);
    setOneAttribute(dump, "T_START", t0);
    setOneAttribute(dump, "DT", dt);
    setOneAttribute(dump, "T_STEPS", nt);
    Attributes::setString(*dump->findAttribute("FILE_NAME"), filename);
    Attributes::setString(*dump->findAttribute("COORDINATE_SYSTEM"), "cYLindriCAL");
}
TEST(DumpEMFieldsTest, ConstructorDestructor) {
    OpalTestUtilities::SilenceTest silencer;

    // neither in the set and grid is null
    DumpEMFields* dump1 = new DumpEMFields();
    delete dump1;
    // grid is not null and it is in the set
    DumpEMFields* dump2 = new DumpEMFields();
    setAttributesCart(dump2, 1., 1., 1.,   1., 1., 1.,   1., 1., 1.,   1., 1., 1.,  "/dev/null");
    dump2->execute();
    delete dump2;
}

void execute_throws(DumpEMFields* dump, std::string reason) {
    try {
        dump->execute();
        EXPECT_TRUE(false) << reason;
    } catch (OpalException& exc) {
        // pass;
    }
}

TEST(DumpEMFieldsTest, executeTest) {
    OpalTestUtilities::SilenceTest silencer;

    // dump the fields
    DumpEMFields dump1;
    execute_throws(&dump1, "should throw due to nsteps < 1");
    setAttributesCart(&dump1, 1., 1., 1.,   1., 1., 1.,   1., 1., 1.,   1., 1., 1.,  "/dev/null", true);
    dump1.execute();  // should be okay (normal)
    setAttributesCart(&dump1, 1., 1., 1.,   1., 1., 1.,   1., 1., 1.,   1., 1., 1.,  "/dev/null", false);
    dump1.execute();  // should be okay (normal)
    setAttributesCart(&dump1, -1., -1., 1.,   -1., -1., 1.,   -1., -1., 1.,   1., 1., 1.,  "/dev/null");
    dump1.execute();  // should be okay (-ve step is okay)
    setAttributesCart(&dump1, -1., -1., 0.,   -1., -1., 1.,   -1., -1., 1.,   1., 1., 1.,  "/dev/null");
    execute_throws(&dump1, "should throw due to nsteps x < 1");
    setAttributesCart(&dump1, -1., -1., 1.,   -1., -1., 0.,   -1., -1., 1.,   1., 1., 1.,  "/dev/null");
    execute_throws(&dump1, "should throw due to nsteps y < 1");
    setAttributesCart(&dump1, -1., -1., 1.,   -1., -1., 1.,   -1., -1., 0.,   1., 1., 1.,  "/dev/null");
    execute_throws(&dump1, "should throw due to nsteps z < 1");
    setAttributesCart(&dump1, -1., -1., 1.,   -1., -1., 1.,   -1., -1., 1.,   1., 1., 0.,  "/dev/null");
    execute_throws(&dump1, "should throw due to nsteps t < 1");
    setAttributesCart(&dump1, -1., -1., 1.,   -1., -1., 1.,   -1., -1., 1.5,   1., 1., 1.,  "/dev/null");
    execute_throws(&dump1, "should throw due to nsteps not integer");
}

void clear_files() {
    size_t n_str_array = 0;
    std::string str_array[5] = {"test1", "test2", "test3", "test4", "testCyl"};
    for (size_t i = 0; i < n_str_array; ++i) {
        if (fopen(str_array[i].c_str(), "r") != NULL) {
            remove(str_array[i].c_str());
        }
    }
}

TEST(DumpEMFieldsTest, writeFieldsCartTest) {
    OpalTestUtilities::SilenceTest silencer;

    clear_files();
    DumpEMFields dump1;
    setAttributesCart(&dump1, 1., 1., 1.,   1., 1., 1.,   1., 1., 1.,   1., 1., 1.,  "test1");
    dump1.execute();
    DumpEMFields dump2;
    setAttributesCart(&dump2, 1., 1., 1.,   1., 1., 1.,   1., 1., 1.,   1., 1., 1.,  "test2");
    dump2.execute();
    DumpEMFields dump3;
    setAttributesCart(&dump3, 1., 1., 1.,   1., 1., 1.,   1., 1., 1.,   1., 1., 1.,  "test3");
    // note we don't execute dump3; so it should not be written
    DumpEMFields dump4;
    setAttributesCart(&dump4, 0.1, 0.1, 3.,   -0.1, 0.2, 2.,   0.2, 0.3, 2.,   1., 1., 2.,  "test4");
    dump4.execute();
    MockComponent comp;
    try {
        DumpEMFields::writeFields(&comp);
    } catch (OpalException& exc) {
        EXPECT_TRUE(false) << "Threw OpalException on writefields: " << exc.what() << std::endl;;
    }
    std::ifstream fin1("test1");
    EXPECT_TRUE(fin1.good());
    std::ifstream fin2("test2");
    EXPECT_TRUE(fin2.good());
    std::ifstream fin3("test3");
    EXPECT_FALSE(fin3.good());  // does not exist
    std::ifstream fin4("test4");
    EXPECT_TRUE(fin4.good());
    int n_lines;
    fin4 >> n_lines;
    EXPECT_EQ(n_lines, 24);
    std::string test_line;
    for (size_t i = 0; i < 12; ++i) {
        std::getline(fin4, test_line);
    }
    std::vector<double> line(10, 0.);
    double tol = 1e-9;
    for (size_t line_index = 0; line_index < 24; ++line_index) {
        for (size_t i = 0; i < 10; ++i) {
            fin4 >> line[i];
        }
        if (line_index == 0) {
            EXPECT_NEAR(line[0], 0.1, tol);
            EXPECT_NEAR(line[1], -0.1, tol);
            EXPECT_NEAR(line[2], 0.2, tol);
            EXPECT_NEAR(line[3], 1., tol);
        }
        if (line[1] < 0.) {
            EXPECT_NEAR(line[4], line[0], tol);
            EXPECT_NEAR(line[5], line[1], tol);
            EXPECT_NEAR(line[6], line[2], tol);
        } else {
            EXPECT_NEAR(line[4], 0., tol);
            EXPECT_NEAR(line[5], 0., tol);
            EXPECT_NEAR(line[6], 0., tol);
        }
        EXPECT_NEAR(line[7], -line[4], tol);
        EXPECT_NEAR(line[8], -line[5], tol);
        EXPECT_NEAR(line[9], -line[6], tol);
    }
    // clear_files();
}

TEST(DumpEMFieldsTest, writeFieldsCylTest) {
    OpalTestUtilities::SilenceTest silencer;

    clear_files();
    DumpEMFields dump;
    setAttributesCyl(&dump, 0.1, 0.1, 3.,   90., 45., 16,   0.2, 0.3, 2.,   1., 1., 2.,  "testCyl");
    dump.execute();
    // depending on execution order, this might write cartesian tests as well... never mind
    MockComponent comp;
    try {
        DumpEMFields::writeFields(&comp);
    } catch (OpalException& exc) {
        EXPECT_TRUE(false) << "Threw OpalException on writefields: " << exc.what() << std::endl;;
    }
    std::ifstream fin("testCyl");
    EXPECT_TRUE(fin.good());
    int n_lines;
    fin >> n_lines;
    EXPECT_EQ(n_lines, 192);
    std::string test_line;
    for (size_t i = 0; i < 12; ++i) {
        std::getline(fin, test_line);
    }
    std::vector<double> line(10, 0.);
    double tol = 1e-9;
    for (size_t line_index = 0; line_index < 24; ++line_index) {
        for (size_t i = 0; i < 10; ++i) {
            fin >> line[i];
        }
        if (line_index == 0) {
            EXPECT_NEAR(line[0], 0.1, tol);
            EXPECT_NEAR(line[1], 90, tol);
            EXPECT_NEAR(line[2], 0.2, tol);
            EXPECT_NEAR(line[3], 1., tol);
        }
        while (line[1] > 360.) {
            line[1] -= 360.;
        }
        if (line[1] < 90. || line[1] > 270.) {
            EXPECT_NEAR(line[4]*line[4]+line[5]*line[5], line[0]*line[0], tol);
            EXPECT_NEAR(line[6], line[2], tol);
        } else {
            EXPECT_NEAR(line[4], 0., tol);
            EXPECT_NEAR(line[5], 0., tol);
            EXPECT_NEAR(line[6], 0., tol);
        }
        EXPECT_NEAR(line[7], -line[4], tol);
        EXPECT_NEAR(line[8], -line[5], tol);
        EXPECT_NEAR(line[9], -line[6], tol);
    }
    clear_files();

    // EXPECT_TRUE(false) << "Do DumpEMFields cylindrical documentation!";
}
}
