//
// Class DumpEMFields
//   DumpEMFields dumps the dynamically changing fields of a Ring in a user-
//   defined grid.
//
// Copyright (c) 2017, Chris Rogers
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include "AbsBeamline/Component.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BasicActions/DumpEMFields.h"
#include "BeamlineGeometry/NullGeometry.h"
#include "Elements/OpalBeamline.h"
#include "Fields/NullField.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"
#include "gtest/gtest.h"

namespace {

    /** Mockup for an Opal Component (e.g. field object). The idea is to test
     *  field lookup routines and placement routines and the like by generating a
     *  "fake" component.
     */
    class MockComponent : public Component {
    public:
        MockComponent() : Component("MockComponent"), field_m(std::make_unique<NullField>()) {}
        MockComponent(const MockComponent& /*rhs*/)
            : Component("MockComponent"), field_m(std::make_unique<NullField>()) {}
        ~MockComponent() override = default;
        void accept(BeamlineVisitor&) const override {}
        ElementBase* clone() const override { return new MockComponent(*this); }
        EMField& getField() override { return *field_m; }
        EMField& getField() const override { return *field_m; }
        bool apply() override { return false; }
        bool apply(
                const size_t& /*i*/, const double& /*t*/, Vector_t<double, 3>& /*E*/,
                Vector_t<double, 3>& /*B*/) override {
            return false;
        }
        bool apply(
                const Vector_t<double, 3>& r, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
                Vector_t<double, 3>& E, Vector_t<double, 3>& B) override {
            if (r(0) < 0. || r(0) > 1. || r(1) < -1. || r(1) > 0. || r(2) < 0. || r(2) > 1.) {
                return true;  // isOutOfBounds
            }
            B(0) = r(0);
            B(1) = r(1);
            B(2) = r(2);
            E(0) = -r(0);
            E(1) = -r(1);
            E(2) = -r(2);
            return false;  // NOT isOutOfBounds
        }
        bool applyToReferenceParticle(
                const Vector_t<double, 3>& /*r*/, const Vector_t<double, 3>& /*P*/,
                const double& /*t*/, Vector_t<double, 3>& /*E*/,
                Vector_t<double, 3>& /*B*/) override {
            return false;
        }
        void initialise(PartBunch_t*, double&, double&) override {}
        void finalise() override {}
        bool bends() const override { return true; }
        void getDimensions(double&, double&) const override {}

        BGeometryBase& getGeometry() override { return geometry_m; }
        const BGeometryBase& getGeometry() const override { return geometry_m; }

        NullGeometry geometry_m;
        std::unique_ptr<NullField> field_m;
    };

    void setOneAttribute(DumpEMFields* dump, const std::string& name, const double value) {
        Attributes::setReal(*dump->findAttribute(name), value);
    }

    void setAttributesCart(
            DumpEMFields* dump, const double x0, const double dx, const double nx, const double y0,
            const double dy, const double ny, const double z0, const double dz, const double nz,
            const double t0, const double dt, const double nt, const std::string& filename) {
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
    }

    void setOriginCyl(DumpEMFields* dump, const double xc, const double yc, const double zc) {
        setOneAttribute(dump, "CYL_ORIGIN_X", xc);
        setOneAttribute(dump, "CYL_ORIGIN_Y", yc);
        setOneAttribute(dump, "CYL_ORIGIN_Z", zc);
    }

    void setAttributesCyl(
            DumpEMFields* dump, const double r0, const double dr, const double nr,
            const double phi0, const double dphi, const double nphi, const double y0,
            const double dy, const double ny, const double t0, const double dt, const double nt,
            const std::string& filename) {
        setOneAttribute(dump, "R_START", r0);
        setOneAttribute(dump, "DR", dr);
        setOneAttribute(dump, "R_STEPS", nr);
        setOneAttribute(dump, "PHI_START", phi0);
        setOneAttribute(dump, "DPHI", dphi);
        setOneAttribute(dump, "PHI_STEPS", nphi);
        setOneAttribute(dump, "Y_START", y0);
        setOneAttribute(dump, "DY", dy);
        setOneAttribute(dump, "Y_STEPS", ny);
        setOneAttribute(dump, "T_START", t0);
        setOneAttribute(dump, "DT", dt);
        setOneAttribute(dump, "T_STEPS", nt);
        Attributes::setString(*dump->findAttribute("FILE_NAME"), filename);
        Attributes::setPredefinedString(*dump->findAttribute("COORDINATE_SYSTEM"), "cYLindriCAL");
    }

    TEST(TestDumpEMFields, ConstructorDestructor) {
        // neither in the set and grid is null
        const auto* dump1 = new DumpEMFields();
        delete dump1;
        // grid is not null and it is in the set
        auto* dump2 = new DumpEMFields();
        setAttributesCart(dump2, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., "/dev/null");
        dump2->execute();
        delete dump2;
        DumpEMFields::clearDumps();
    }

    void execute_throws(DumpEMFields* dump, const std::string& reason) {
        try {
            dump->execute();
            EXPECT_TRUE(false) << reason;
        } catch (OpalException&) {
            // pass;
        }
    }

    TEST(TestDumpEMFields, executeTest) {
        gmsg = new Inform(nullptr, -1);
        // dump the fields
        DumpEMFields dump1;
        execute_throws(&dump1, "should throw due to nsteps < 1");
        setAttributesCart(&dump1, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., "/dev/null");
        dump1.execute();  // should be okay (normal)
        std::set<std::shared_ptr<Component>> elements;
        elements.insert(std::make_shared<MockComponent>());
        EXPECT_NO_THROW(DumpEMFields::writeFields(elements));
        setAttributesCart(
                &dump1, -1., -1., 1., -1., -1., 1., -1., -1., 1., 1., 1., 1., "/dev/null");
        dump1.execute();  // should be okay (-ve step is okay)
        setAttributesCart(
                &dump1, -1., -1., 0., -1., -1., 1., -1., -1., 1., 1., 1., 1., "/dev/null");
        execute_throws(&dump1, "should throw due to nsteps x < 1");
        setAttributesCart(
                &dump1, -1., -1., 1., -1., -1., 0., -1., -1., 1., 1., 1., 1., "/dev/null");
        execute_throws(&dump1, "should throw due to nsteps y < 1");
        setAttributesCart(
                &dump1, -1., -1., 1., -1., -1., 1., -1., -1., 0., 1., 1., 1., "/dev/null");
        execute_throws(&dump1, "should throw due to nsteps z < 1");
        setAttributesCart(
                &dump1, -1., -1., 1., -1., -1., 1., -1., -1., 1., 1., 1., 0., "/dev/null");
        execute_throws(&dump1, "should throw due to nsteps t < 1");
        setAttributesCart(
                &dump1, -1., -1., 1., -1., -1., 1., -1., -1., 1.5, 2., 2., 1., "/dev/null");
        execute_throws(&dump1, "should throw due to nsteps not integer");
        DumpEMFields::clearDumps();
    }

    void clear_files(std::set<std::string> const& files) {
        std::string auxDirectory = OpalData::getInstance()->getAuxiliaryOutputDirectory();

        for (const std::string& fname : files) {
            std::filesystem::remove(Util::combineFilePath({auxDirectory, fname}));
        }
    }

    TEST(TestDumpEMFields, writeFieldsCartTest) {
        std::string auxDirectory = OpalData::getInstance()->getAuxiliaryOutputDirectory();
        std::filesystem::create_directory(auxDirectory);
        gmsg               = new Inform(nullptr, -1);
        std::string fname1 = "test5";
        std::string fname2 = "test6";
        std::string fname3 = "test7";
        std::string fname4 = "test8";

        clear_files({fname1, fname2, fname3, fname4});
        DumpEMFields dump1;
        setAttributesCart(&dump1, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., fname1);
        dump1.execute();
        DumpEMFields dump2;
        setAttributesCart(&dump2, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., fname2);
        dump2.execute();
        DumpEMFields dump3;
        setAttributesCart(&dump3, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., fname3);
        // note we don't execute dump3; so it should not be written
        DumpEMFields dump4;
        setAttributesCart(&dump4, 0.1, 0.1, 3., -0.1, 0.2, 2., 0.2, 0.3, 2., 1., 1., 2., fname4);
        dump4.execute();
        std::set<std::shared_ptr<Component>> elements;
        elements.insert(std::make_shared<MockComponent>());
        try {
            DumpEMFields::writeFields(elements);
        } catch (OpalException& exc) {
            EXPECT_TRUE(false) << "Threw OpalException on writefields: " << exc.what() << std::endl;
        }
        std::ifstream fin1(Util::combineFilePath({auxDirectory, fname1}));
        EXPECT_TRUE(fin1.good());
        std::ifstream fin2(Util::combineFilePath({auxDirectory, fname2}));
        EXPECT_TRUE(fin2.good());
        std::ifstream fin3(Util::combineFilePath({auxDirectory, fname3}));
        EXPECT_FALSE(fin3.good());  // does not exist
        std::ifstream fin4(Util::combineFilePath({auxDirectory, fname4}));
        EXPECT_TRUE(fin4.good());
        int n_lines;
        fin4 >> n_lines;
        EXPECT_EQ(n_lines, 24);
        std::string test_line;
        for (size_t i = 0; i < 12; ++i) {
            std::getline(fin4, test_line);
        }
        std::vector line(10, 0.);
        for (size_t line_index = 0; line_index < 24; ++line_index) {
            double tol = 1e-9;
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
        clear_files({fname1, fname2, fname3, fname4});
        DumpEMFields::clearDumps();
    }

    TEST(TestDumpEMFields, writeFieldsCylTest) {
        std::string auxDirectory = OpalData::getInstance()->getAuxiliaryOutputDirectory();
        std::filesystem::create_directory(auxDirectory);
        gmsg = new Inform(nullptr, -1);

        std::string fnameCyl = "testCyl";

        clear_files({fnameCyl});
        DumpEMFields dump;
        setAttributesCyl(
                &dump, 0.1, 0.1, 3., 90. * Units::deg2rad, 45. * Units::deg2rad, 16, 0.2, 0.3, 2.,
                1., 1., 2., fnameCyl);
        dump.execute();
        // depending on execution order, this might write cartesian tests as well... never mind
        std::set<std::shared_ptr<Component>> elements;
        elements.insert(std::make_shared<MockComponent>());
        try {
            DumpEMFields::writeFields(elements);
        } catch (OpalException& exc) {
            EXPECT_TRUE(false) << "Threw OpalException on writefields: " << exc.what() << std::endl;
        }
        std::ifstream fin(Util::combineFilePath({auxDirectory, fnameCyl}));
        EXPECT_TRUE(fin.good());
        int n_lines;
        fin >> n_lines;
        EXPECT_EQ(n_lines, 192);
        std::string test_line;
        for (size_t i = 0; i < 12; ++i) {
            std::getline(fin, test_line);
        }
        std::vector line(10, 0.);
        for (size_t line_index = 0; line_index < 24; ++line_index) {
            constexpr double tol = 1e-9;
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
            if (line[1] < 90.) {
                EXPECT_NEAR(line[4] * line[4] + line[5] * line[5], line[0] * line[0], tol);
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
        clear_files({fnameCyl});

        // EXPECT_TRUE(false) << "Do DumpEMFields cylindrical documentation!";
        DumpEMFields::clearDumps();
    }

    TEST(TestDumpEMFields, writeFieldsCylOriginTest) {
        std::string auxDirectory = OpalData::getInstance()->getAuxiliaryOutputDirectory();
        std::filesystem::create_directory(auxDirectory);
        gmsg = new Inform(nullptr, -1);

        std::string fnameCyl1 = "testCyl1";
        std::string fnameCyl2 = "testCyl2";

        clear_files({fnameCyl1, fnameCyl2});
        DumpEMFields dump1;
        setAttributesCyl(
                &dump1, 0.1, 0., 1., 45. * Units::deg2rad, 0., 1, -0.1, 0., 1., 1., 0., 1.,
                fnameCyl1);
        dump1.execute();
        DumpEMFields dump2;
        setAttributesCyl(
                &dump2, 0.1, 0., 1., 45. * Units::deg2rad, 0., 1, -0.1, 0., 1., 1., 0., 1.,
                fnameCyl2);
        setOriginCyl(&dump2, 0.01, 0.02, 0.03);
        dump2.execute();
        // depending on execution order, this might write cartesian tests as well... never mind
        std::set<std::shared_ptr<Component>> elements;
        elements.insert(std::make_shared<MockComponent>());
        try {
            DumpEMFields::writeFields(elements);
        } catch (OpalException& exc) {
            EXPECT_TRUE(false) << "Threw OpalException on writefields: " << exc.what() << std::endl;
        }
        // Read the origin (0,0,0) file
        std::ifstream fin1(Util::combineFilePath({auxDirectory, fnameCyl1}));
        EXPECT_TRUE(fin1.good());
        int n_lines;
        fin1 >> n_lines;
        EXPECT_EQ(n_lines, 1);
        std::string test_line;
        for (size_t i = 0; i < 12; ++i) {
            std::getline(fin1, test_line);
        }
        std::vector line1(10, 0.);
        for (size_t i = 0; i < 10; ++i) {
            fin1 >> line1[i];
        }
        // Read the origin (1,2,3) file
        std::ifstream fin2(Util::combineFilePath({auxDirectory, fnameCyl2}));
        EXPECT_TRUE(fin2.good());
        fin2 >> n_lines;
        EXPECT_EQ(n_lines, 1);
        for (size_t i = 0; i < 12; ++i) {
            std::getline(fin2, test_line);
        }
        std::vector line2(10, 0.);
        for (size_t i = 0; i < 10; ++i) {
            fin2 >> line2[i];
        }
        EXPECT_EQ(line1[0], line2[0]);
        EXPECT_EQ(line1[1], line2[1]);
        EXPECT_EQ(line1[2], line2[2]);
        EXPECT_EQ(line1[3], line2[3]);
        EXPECT_NE(line1[4], line2[4]);
        EXPECT_NE(line1[5], line2[5]);
        EXPECT_NE(line1[6], line2[6]);
        EXPECT_NE(line1[7], line2[7]);
        EXPECT_NE(line1[8], line2[8]);
        EXPECT_NE(line1[9], line2[9]);
        clear_files({fnameCyl1, fnameCyl2});
        DumpEMFields::clearDumps();
    }

    TEST(TestDumpEMFields, BadGrid) {
        std::string auxDirectory = OpalData::getInstance()->getAuxiliaryOutputDirectory();
        std::filesystem::create_directory(auxDirectory);
        gmsg = new Inform(nullptr, -1);

        std::string fnameCyl = "testCyl";

        clear_files({fnameCyl});
        DumpEMFields dump;
        setAttributesCyl(
                &dump, 0.1, 0.1, 3., 90. * Units::deg2rad, 45. * Units::deg2rad, 16, 0.2, 0.3, 2.,
                1., 1., 2., fnameCyl);
        dump.execute();
        DumpEMFields::failGrid();
        std::set<std::shared_ptr<Component>> elements;
        elements.insert(std::make_shared<MockComponent>());
        EXPECT_ANY_THROW(DumpEMFields::writeFields(elements));
        DumpEMFields::clearDumps();
    }

    TEST(TestDumpEMFields, CouldNotOpenFile) {
        std::string auxDirectory = OpalData::getInstance()->getAuxiliaryOutputDirectory();
        std::filesystem::create_directory(auxDirectory);
        gmsg                 = new Inform(nullptr, -1);
        std::string fnameCyl = ".";  // Opening for write to a directory always fails
        DumpEMFields dump;
        setAttributesCyl(
                &dump, 0.1, 0.1, 3., 90. * Units::deg2rad, 45. * Units::deg2rad, 16, 0.2, 0.3, 2.,
                1., 1., 2., fnameCyl);
        dump.execute();
        std::set<std::shared_ptr<Component>> elements;
        elements.insert(std::make_shared<MockComponent>());
        EXPECT_ANY_THROW(DumpEMFields::writeFields(elements));
        DumpEMFields::clearDumps();
    }

    TEST(TestDumpEMFields, WriteFailure) {
        std::string auxDirectory = OpalData::getInstance()->getAuxiliaryOutputDirectory();
        std::filesystem::create_directory(auxDirectory);
        gmsg                 = new Inform(nullptr, -1);
        std::string fnameCyl = "testCyl";
        DumpEMFields dump;
        setAttributesCyl(
                &dump, 0.1, 0.1, 3., 90. * Units::deg2rad, 45. * Units::deg2rad, 16, 0.2, 0.3, 2.,
                1., 1., 2., fnameCyl);
        dump.execute();
        DumpEMFields::failWrite();
        std::set<std::shared_ptr<Component>> elements;
        elements.insert(std::make_shared<MockComponent>());
        EXPECT_ANY_THROW(DumpEMFields::writeFields(elements));
        DumpEMFields::clearDumps();
    }
}  // namespace
