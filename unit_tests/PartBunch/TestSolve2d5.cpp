//
// Copyright (c) 2008 - 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
//
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

#include <gtest/gtest.h>
#include <cmath>
#include <csignal>
#include <filesystem>
#include <memory>
#include <random>
#include <string>
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "Lines/Line.h"
#include "PartBunch/PartBunch.h"
#include "PartBunch/Solve2d5.h"
#include "PartBunch/SubField.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

namespace {

    constexpr bool VerboseTest = true;

    class TestableFieldSolverCmd : public FieldSolverCmd {
    public:
        void setType(const std::string& t) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::TYPE], t);
        }
    };

    using PartBunch_t         = PartBunch<double, 3>;
    using Solve2d5_t          = Solve2d5<double>;
    using ParticleContainer_t = PartBunch_t::ParticleContainer_t;
    using Point_t             = ippl::Vector<double, 3>;

    class ParticleBaseDiagnostic : public Solve2d5_t::NullDiagnostic {
    public:
        explicit ParticleBaseDiagnostic(const ParticleContainer_t* pc) {
            r_m       = Solve2d5_t::VectorView_t("fsr", pc->R.getParticleCount());
            p_m       = Solve2d5_t::VectorView_t("fsp", pc->R.getParticleCount());
            invalid_m = Solve2d5_t::BooleanView_t("fsinv", pc->R.getParticleCount());
        }

        std::tuple<std::vector<Vector_t<double, 3>>, std::vector<Vector_t<double, 3>>>
        getParticles() const {
            const auto rHost       = Kokkos::create_mirror(r_m);
            const auto pHost       = Kokkos::create_mirror(p_m);
            const auto invalidHost = Kokkos::create_mirror(invalid_m);
            Kokkos::deep_copy(rHost, r_m);
            Kokkos::deep_copy(pHost, p_m);
            Kokkos::deep_copy(invalidHost, invalid_m);
            std::vector<Vector_t<double, 3>> r;
            std::vector<Vector_t<double, 3>> p;
            for (size_t i = 0; i < r_m.extent(0); ++i) {
                if (!invalidHost(i)) {
                    r.push_back(rHost(i));
                    p.push_back(pHost(i));
                }
            }
            return std::make_tuple(r, p);
        }

        Solve2d5_t::VectorView_t r_m;
        Solve2d5_t::VectorView_t p_m;
        Solve2d5_t::BooleanView_t invalid_m;
    };

    class FrenetSerretScatterDiagnostic : public ParticleBaseDiagnostic {
    public:
        explicit FrenetSerretScatterDiagnostic(const ParticleContainer_t* pc)
            : ParticleBaseDiagnostic(pc) {}

        KOKKOS_FUNCTION void frenetSerretScatter(
                const size_t n, const Solve2d5<double>::Vector3D_t& r,
                const Solve2d5<double>::Vector3D_t& p, const bool invalid) const {
            r_m[n]       = r;
            p_m[n]       = p;
            invalid_m[n] = invalid;
        }
    };

    class FrenetSerretGatherDiagnostic : public ParticleBaseDiagnostic {
    public:
        explicit FrenetSerretGatherDiagnostic(const ParticleContainer_t* pc)
            : ParticleBaseDiagnostic(pc) {}

        KOKKOS_FUNCTION void frenetSerretGather(
                const size_t n, const Solve2d5<double>::Vector3D_t& r,
                const Solve2d5<double>::Vector3D_t& p, const bool invalid) const {
            r_m[n]       = r;
            p_m[n]       = p;
            invalid_m[n] = invalid;
        }
    };

    class BoostToBeamDiagnostic : public ParticleBaseDiagnostic {
    public:
        explicit BoostToBeamDiagnostic(const ParticleContainer_t* pc)
            : ParticleBaseDiagnostic(pc) {}

        KOKKOS_FUNCTION void boostToBeam(
                const size_t n, const Solve2d5<double>::Vector3D_t& r,
                const Solve2d5<double>::Vector3D_t& p, const bool invalid) const {
            r_m[n]       = r;
            p_m[n]       = p;
            invalid_m[n] = invalid;
        }
    };

    class DeboostFromBeamDiagnostic : public ParticleBaseDiagnostic {
    public:
        explicit DeboostFromBeamDiagnostic(const ParticleContainer_t* pc)
            : ParticleBaseDiagnostic(pc) {}

        KOKKOS_FUNCTION void deboostFromBeam(
                const size_t n, const Solve2d5<double>::Vector3D_t& r,
                const Solve2d5<double>::Vector3D_t& p, const bool invalid) const {
            r_m[n]       = r;
            p_m[n]       = p;
            invalid_m[n] = invalid;
        }
    };

    class GatheredFieldBaseDiagnostic : public Solve2d5_t::NullDiagnostic {
    public:
        explicit GatheredFieldBaseDiagnostic(const ParticleContainer_t* pc) {
            e_m       = Solve2d5_t::VectorView_t("e", pc->R.getParticleCount());
            b_m       = Solve2d5_t::VectorView_t("b", pc->R.getParticleCount());
            invalid_m = Solve2d5_t::BooleanView_t("fsinv", pc->R.getParticleCount());
        }

        std::tuple<std::vector<Vector_t<double, 3>>, std::vector<Vector_t<double, 3>>>
        getParticleFields() const {
            const auto eHost       = Kokkos::create_mirror(e_m);
            const auto bHost       = Kokkos::create_mirror(b_m);
            const auto invalidHost = Kokkos::create_mirror(invalid_m);
            Kokkos::deep_copy(eHost, e_m);
            Kokkos::deep_copy(bHost, b_m);
            Kokkos::deep_copy(invalidHost, invalid_m);
            std::vector<Vector_t<double, 3>> e;
            std::vector<Vector_t<double, 3>> b;
            for (size_t i = 0; i < e_m.extent(0); ++i) {
                if (!invalidHost(i)) {
                    e.push_back(eHost(i));
                    b.push_back(bHost(i));
                }
            }
            return std::make_tuple(e, b);
        }

        Solve2d5_t::VectorView_t e_m;
        Solve2d5_t::VectorView_t b_m;
        Solve2d5_t::BooleanView_t invalid_m;
    };

    class GatherEFieldDiagnostic : public GatheredFieldBaseDiagnostic {
    public:
        explicit GatherEFieldDiagnostic(const ParticleContainer_t* pc)
            : GatheredFieldBaseDiagnostic(pc) {}

        KOKKOS_FUNCTION void gatherEField(
                const size_t n, const Solve2d5<double>::Vector3D_t& e,
                const Solve2d5<double>::Vector3D_t& b, const bool invalid) const {
            e_m[n]       = e;
            b_m[n]       = b;
            invalid_m[n] = invalid;
        }
    };

    class DeboostedDiagnostic : public GatheredFieldBaseDiagnostic {
    public:
        explicit DeboostedDiagnostic(const ParticleContainer_t* pc)
            : GatheredFieldBaseDiagnostic(pc) {}

        KOKKOS_FUNCTION void deboostFromBeam(
                const size_t n, const Solve2d5<double>::Vector3D_t& e,
                const Solve2d5<double>::Vector3D_t& b, const bool invalid) const {
            e_m[n]       = e;
            b_m[n]       = b;
            invalid_m[n] = invalid;
        }
    };

    class LongitudinalFieldDiagnostic : public GatheredFieldBaseDiagnostic {
    public:
        explicit LongitudinalFieldDiagnostic(const ParticleContainer_t* pc)
            : GatheredFieldBaseDiagnostic(pc) {}

        KOKKOS_FUNCTION void longitudinalField(
                const size_t n, const Solve2d5<double>::Vector3D_t& e,
                const Solve2d5<double>::Vector3D_t& b, const bool invalid) const {
            e_m[n]       = e;
            b_m[n]       = b;
            invalid_m[n] = invalid;
        }
    };

    class LabFrameFieldsDiagnostic : public GatheredFieldBaseDiagnostic {
    public:
        explicit LabFrameFieldsDiagnostic(const ParticleContainer_t* pc)
            : GatheredFieldBaseDiagnostic(pc) {}

        KOKKOS_FUNCTION void labFrameFields(
                const size_t n, const Solve2d5<double>::Vector3D_t& e,
                const Solve2d5<double>::Vector3D_t& b, const bool invalid) const {
            e_m[n]       = e;
            b_m[n]       = b;
            invalid_m[n] = invalid;
        }
    };

    class ScatterChargeDiagnostic : public Solve2d5_t::NullDiagnostic {
    public:
        explicit ScatterChargeDiagnostic(const Solve2d5_t::ScalarGridView3D_t& rho)
            : NullDiagnostic() {
            Kokkos::resize(rhoView, rho.extent(0), rho.extent(1), rho.extent(2));
        }

        KOKKOS_FUNCTION void scatterCharge(const Solve2d5_t::ScalarGridView3D_t& rho) const {
            Kokkos::deep_copy(rhoView, rho);
        }

        Field_t<3>::view_type rhoView;
    };

    class ScatterChargeDensityDiagnostic : public Solve2d5_t::NullDiagnostic {
    public:
        explicit ScatterChargeDensityDiagnostic(const Solve2d5_t::ScalarGridView3D_t& rho)
            : NullDiagnostic() {
            Kokkos::resize(rhoView, rho.extent(0), rho.extent(1), rho.extent(2));
        }

        KOKKOS_FUNCTION void scatterChargeDensity(const Solve2d5_t::ScalarGridView3D_t& rho) const {
            Kokkos::deep_copy(rhoView, rho);
        }

        Field_t<3>::view_type rhoView;
    };

    class TotalDensityDiagnostic : public Solve2d5_t::NullDiagnostic {
    public:
        explicit TotalDensityDiagnostic(const Solve2d5_t::LineDensityView_t& lineDensity)
            : NullDiagnostic() {
            lineDensityView_m =
                    Solve2d5_t::LineDensityView_t(lineDensity.label(), lineDensity.extent(0));
        }

        KOKKOS_FUNCTION void totalDensity(const Solve2d5_t::LineDensityView_t& lineDensity) const {
            Kokkos::deep_copy(lineDensityView_m, lineDensity);
        }

        Solve2d5_t::LineDensityView_t lineDensityView_m;
    };

    class LineDensityDiagnostic : public Solve2d5_t::NullDiagnostic {
    public:
        explicit LineDensityDiagnostic(const Solve2d5_t::LineDensityView_t& lineDensity)
            : NullDiagnostic() {
            Kokkos::resize(lineDensityView_m, lineDensity.extent(0));
        }

        KOKKOS_FUNCTION void lineDensity(const Solve2d5_t::LineDensityView_t& lineDensity) const {
            Kokkos::deep_copy(lineDensityView_m, lineDensity);
        }

        Solve2d5_t::LineDensityView_t lineDensityView_m;
    };

    class LineDensityGradientDiagnostic : public Solve2d5_t::NullDiagnostic {
    public:
        explicit LineDensityGradientDiagnostic(const Solve2d5_t::LineDensityView_t& lineDensity)
            : NullDiagnostic() {
            Kokkos::resize(lineDensityView_m, lineDensity.extent(0));
        }

        KOKKOS_FUNCTION void lineDensityGradient(
                const Solve2d5_t::LineDensityView_t& lineDensity) const {
            Kokkos::deep_copy(lineDensityView_m, lineDensity);
        }

        Solve2d5_t::LineDensityView_t lineDensityView_m;
    };

    class EFieldDiagnostic : public Solve2d5_t::NullDiagnostic {
    public:
        explicit EFieldDiagnostic(const VField_t<T, 3>::view_type& eField) : NullDiagnostic() {
            Kokkos::resize(eFieldView_m, eField.extent(0), eField.extent(1), eField.extent(2));
        }

        KOKKOS_FUNCTION void eField(const VField_t<T, 3>::view_type& eField) const {
            Kokkos::deep_copy(eFieldView_m, eField);
        }

        VField_t<T, 3>::view_type eFieldView_m;
    };

    class TestSolve2d5 : public testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);

            // DataSink requires a basename to create *.stat / *.lbal writers.
            OpalData::getInstance()->storeInputFn("unit_test.opal");

            // Many OPAL writers assume `gmsg` is initialized (see SDDSWriter/StatWriter).
            // Unit tests normally don't set this up via Main().
            gmsg = new Inform(nullptr, -1);

            // DataSink::DataSink() constructs HDF5 writers when enabled, but the unit
            // test doesn't have an H5PartWrapper. Disable HDF5 for this test.
            Options::enableHDF5 = false;
        }

        static void TearDownTestSuite() {
            delete gmsg;
            gmsg = nullptr;
            ippl::finalize();
        }

        static constexpr double nx                 = 8;
        static constexpr double ny                 = 9;
        static constexpr double nz                 = 10;
        static constexpr size_t kDefaultNParticles = 64;

        void SetUp() override {
            // Remove any reference path file
            if (std::filesystem::exists("data/unit_test_DesignPath.dat")) {
                std::filesystem::remove("data/unit_test_DesignPath.dat");
            }
            // Keep mesh small so scatter/solve/gather are quick.
            fsCmd = std::make_shared<TestableFieldSolverCmd>();
            fsCmd->setType("NONE");
            fsCmd->setNX(nx);
            fsCmd->setNY(ny);
            fsCmd->setNZ(nz);

            dataSink       = std::make_shared<DataSink>();
            beam           = std::make_shared<Beam>();
            Beam* testBeam = Beam::find("UNNAMED_BEAM");

            // qi/mi/lbt are used by rho scaling; but with NullSolver we mostly validate
            // "no-throw" and deterministic zero E behavior.
            bunch = std::make_shared<PartBunch_t>(
                    /*qi=*/std::vector{1.0},
                    /*mi=*/std::vector{1.0},
                    /*beams=*/std::vector{testBeam},
                    /*totalParticlesPerBeam=*/std::vector{kDefaultNParticles},
                    /*lbt=*/1.0,
                    /*integration_method=*/"LF2", fsCmd.get(), dataSink.get());
            pc = bunch->getParticleContainer();
        }

        void TearDown() override {
            // Ensure device allocations can be freed between tests.
            bunch.reset();
            dataSink.reset();
            fsCmd.reset();
            pc.reset();
            // Remove reference path file
            if (std::filesystem::exists("data/unit_test_DesignPath.dat")) {
                std::filesystem::remove("data/unit_test_DesignPath.dat");
            }
        }

        void createParticles(
                const std::vector<Vector_t<double, 3>>& r,
                const std::vector<Vector_t<double, 3>>& p) const {
            pc->createParticles(r.size());
            const auto R_host       = pc->R.getHostMirror();
            const auto P_host       = pc->P.getHostMirror();
            const auto dt_host      = pc->dt.getHostMirror();
            const auto E_host       = pc->E.getHostMirror();
            const auto invalid_host = pc->InvalidMask.getHostMirror();
            const double dt         = bunch->getdT();
            const double qi         = pc->getChargePerParticle();
            for (size_t i = 0; i < r.size(); ++i) {
                R_host(i)       = r[i];
                P_host(i)       = p[i];
                dt_host(i)      = dt;
                E_host(i)       = {0.0, 0.0, 0.0};
                invalid_host(i) = false;
            }
            for (size_t i = r.size(); i < invalid_host.extent(0); ++i) {
                invalid_host(i) = true;
            }
            Kokkos::deep_copy(pc->R.getView(), R_host);
            Kokkos::deep_copy(pc->P.getView(), P_host);
            Kokkos::deep_copy(pc->dt.getView(), dt_host);
            Kokkos::deep_copy(pc->E.getView(), E_host);
            Kokkos::deep_copy(pc->InvalidMask.getView(), invalid_host);
            pc->setQ(qi);
            ippl::Comm->barrier();
            Kokkos::fence();
        }

        std::tuple<std::vector<Vector_t<double, 3>>, std::vector<Vector_t<double, 3>>>
        getParticles() const {
            const auto R_host       = pc->R.getHostMirror();
            const auto P_host       = pc->P.getHostMirror();
            const auto invalid_host = pc->InvalidMask.getHostMirror();
            Kokkos::deep_copy(R_host, pc->R.getView());
            Kokkos::deep_copy(P_host, pc->P.getView());
            Kokkos::deep_copy(invalid_host, pc->InvalidMask.getView());
            std::vector<Vector_t<double, 3>> r;
            std::vector<Vector_t<double, 3>> p;
            for (size_t i = 0; i < R_host.extent(0); ++i) {
                if (!invalid_host(i)) {
                    r.push_back(R_host(i));
                    p.push_back(P_host(i));
                }
            }
            return std::make_tuple(r, p);
        }

        void rebuildBunch() {
            Beam* testBeam = Beam::find("UNNAMED_BEAM");
            bunch          = std::make_shared<PartBunch_t>(
                    /*qi=*/std::vector{1.0},
                    /*mi=*/std::vector{1.0},
                    /*beams=*/std::vector{testBeam},
                    /*totalParticlesPerBeam=*/std::vector{kDefaultNParticles},
                    /*lbt=*/1.0,
                    /*integration_method=*/"LF2", fsCmd.get(), dataSink.get());
            pc = bunch->getParticleContainer();
        }

        void makeReferencePathFile(
                const std::filesystem::path& fileName, const std::vector<Point_t>& points,
                const bool noTimeColumn = false) {
            // Emulates the file created by OrbitThreader.  Use noTimeColumn = true to create
            // an invalid file.
            // Make the directory path exist
            std::filesystem::create_directories(fileName.parent_path());
            // Create the file
            std::ofstream file(fileName);
            if (!file.is_open()) {
                throw std::runtime_error("Failed to open file: " + fileName.string());
            }
            // The header
            file << "#    1 – s    2 – Rx    3 - Ry    4 - Rz    5 - Px    6 - Py    7 - Pz    "
                    "8 - Efx    9 - Efy    10 - Efz    11 - Bfx    12 - Bfy    13 - Bfz    "
                    "14 – Ekin   15 - t"
                 << std::endl;
            // The lines
            for (auto& point : points) {
                file << 0 << " " << point.data_m[0] << " " << point.data_m[1] << " "
                     << point.data_m[2] << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0
                     << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0;
                if (!noTimeColumn) {
                    file << " " << 0;
                }
                file << std::endl;
            }
        }

        static void expectParticle(
                size_t index, const std::vector<Vector_t<double, 3>>& rs,
                const std::vector<Vector_t<double, 3>>& ps, const Vector_t<double, 3>& r,
                const Vector_t<double, 3>& p, const double tolerance = 1e-6) {
            SCOPED_TRACE(std::format("Index = {}", index));
            EXPECT_NEAR(r[0], rs[index].data_m[0], tolerance);
            EXPECT_NEAR(r[1], rs[index].data_m[1], tolerance);
            EXPECT_NEAR(r[2], rs[index].data_m[2], tolerance);
            EXPECT_NEAR(p[0], ps[index].data_m[0], tolerance);
            EXPECT_NEAR(p[1], ps[index].data_m[1], tolerance);
            EXPECT_NEAR(p[2], ps[index].data_m[2], tolerance);
        }

        static void expectParticleFields(
                size_t index, const std::vector<Vector_t<double, 3>>& es,
                const std::vector<Vector_t<double, 3>>& bs, const Vector_t<double, 3>& e,
                const Vector_t<double, 3>& b, const double eTolerance = 1e-6,
                const double bTolerance = 1e-16) {
            SCOPED_TRACE(std::format("Index = {}", index));
            EXPECT_NEAR(e[0], es[index].data_m[0], eTolerance);
            EXPECT_NEAR(e[1], es[index].data_m[1], eTolerance);
            EXPECT_NEAR(e[2], es[index].data_m[2], eTolerance);
            EXPECT_NEAR(b[0], bs[index].data_m[0], bTolerance);
            EXPECT_NEAR(b[1], bs[index].data_m[1], bTolerance);
            EXPECT_NEAR(b[2], bs[index].data_m[2], bTolerance);
        }

        struct RhoValue {
            size_t i;
            size_t j;
            size_t k;
            double value;
        };
        static void expectChargeDensity(
                const Solve2d5_t::ScalarGridView3D_t& rho, const std::vector<RhoValue>& expected) {
            // Transfer to host
            const auto hostView = Kokkos::create_mirror_view(rho);
            Kokkos::deep_copy(hostView, rho);
            // Print
            if (VerboseTest) {
                for (size_t k = 0; k < rho.extent(2); ++k) {
                    for (size_t j = 0; j < rho.extent(1); ++j) {
                        std::cout << k << "," << j << ": ";
                        for (size_t i = 0; i < rho.extent(0); ++i) {
                            std::cout << hostView(i, j, k) << " ";
                        }
                        std::cout << std::endl;
                    }
                }
            }
            // Check values
            for (const auto& e : expected) {
                SCOPED_TRACE(std::format("Index: {},{},{}", e.i, e.j, e.k));
                EXPECT_NEAR(hostView(e.i, e.j, e.k), e.value, 1e-6);
            }
        }

        static void expectLineDensity(
                const Solve2d5_t::LineDensityView_t& lineDensity,
                const std::vector<double> expected) {
            // Transfer to host
            const auto hostView = Kokkos::create_mirror_view(lineDensity);
            Kokkos::deep_copy(hostView, lineDensity);
            // Print
            if (VerboseTest) {
                for (size_t k = 0; k < lineDensity.extent(0); ++k) {
                    std::cout << hostView(k) << " ";
                }
            }
            std::cout << std::endl;
            // Check values
            if (!expected.empty()) {
                ASSERT_EQ(hostView.extent(0), expected.size());
                for (size_t i = 0; i < expected.size(); ++i) {
                    SCOPED_TRACE(std::format("Index: {}", i));
                    EXPECT_NEAR(hostView(i), expected[i], 1e-6);
                }
            }
        }

        struct EValue {
            size_t i;
            size_t j;
            size_t k;
            double valueX;
            double valueY;
            double valueZ;
        };
        static void expectEField(
                const Solve2d5_t::VectorGridView3D_t& eField, const std::vector<EValue>& expected) {
            // Transfer to host
            const auto hostView = Kokkos::create_mirror_view(eField);
            Kokkos::deep_copy(hostView, eField);
            // Print the magnitudes
            if (VerboseTest) {
                for (size_t k = 0; k < hostView.extent(2); ++k) {
                    for (size_t j = 0; j < hostView.extent(1); ++j) {
                        std::cout << k << "," << j << ": ";
                        for (size_t i = 0; i < hostView.extent(0); ++i) {
                            std::cout << hostView(i, j, k).Pnorm() << " ";
                        }
                        std::cout << std::endl;
                    }
                }
            }
            // Check values
            for (const auto& e : expected) {
                SCOPED_TRACE(std::format("Index: {},{},{}", e.i, e.j, e.k));
                EXPECT_NEAR(hostView(e.i, e.j, e.k)[0], e.valueX, 1e3);
                EXPECT_NEAR(hostView(e.i, e.j, e.k)[1], e.valueY, 1e3);
                EXPECT_NEAR(hostView(e.i, e.j, e.k)[2], e.valueZ, 1e3);
            }
        }

        std::shared_ptr<TestableFieldSolverCmd> fsCmd;
        std::shared_ptr<DataSink> dataSink;
        std::shared_ptr<Beam> beam;
        std::shared_ptr<PartBunch_t> bunch;
        std::shared_ptr<ParticleContainer_t> pc;
    };

    TEST_F(TestSolve2d5, LoadReferencePath_Missing) {
        fsCmd->setType("FFT2D5");
        EXPECT_ANY_THROW(rebuildBunch());
    }

    TEST_F(TestSolve2d5, LoadReferencePath_ReadFail) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}},
                true);
        fsCmd->setType("FFT2D5");
        EXPECT_ANY_THROW(rebuildBunch());
    }

    TEST_F(TestSolve2d5, LoadReferencePath_Success) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        ASSERT_NO_THROW(rebuildBunch());
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        EXPECT_EQ(solver->getReferencePath().size(), 4);
        EXPECT_EQ(solver->getReferencePath()[0].data_m[0], 0);
        EXPECT_EQ(solver->getReferencePath()[0].data_m[1], 0);
        EXPECT_EQ(solver->getReferencePath()[0].data_m[2], 0);
        EXPECT_EQ(solver->getReferencePath()[1].data_m[0], 0);
        EXPECT_EQ(solver->getReferencePath()[1].data_m[1], 0);
        EXPECT_EQ(solver->getReferencePath()[1].data_m[2], 1);
        EXPECT_EQ(solver->getReferencePath()[2].data_m[0], 1);
        EXPECT_EQ(solver->getReferencePath()[2].data_m[1], 0);
        EXPECT_EQ(solver->getReferencePath()[2].data_m[2], 2);
        EXPECT_EQ(solver->getReferencePath()[3].data_m[0], 0);
        EXPECT_EQ(solver->getReferencePath()[3].data_m[1], 0);
        EXPECT_EQ(solver->getReferencePath()[3].data_m[2], 3);
    }

    TEST_F(TestSolve2d5, SliceSolverSetup) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        ASSERT_NO_THROW(rebuildBunch());
        const auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        EXPECT_EQ(solver->getNumSlices(), nz);
        ASSERT_FALSE(solver->getSliceMesh() == nullptr);
        EXPECT_EQ(solver->getSliceMesh()->getGridsize(0), 8);
        EXPECT_EQ(solver->getSliceMesh()->getGridsize(1), 9);
        EXPECT_EQ(solver->getSliceMesh()->getOrigin()[0], 0);
        EXPECT_EQ(solver->getSliceMesh()->getOrigin()[1], 0);
        EXPECT_DOUBLE_EQ(solver->getSliceMesh()->getMeshSpacing()[0], 0.125);
        EXPECT_NEAR(solver->getSliceMesh()->getMeshSpacing()[1], 0.111111, 1e-6);
        ASSERT_FALSE(solver->getSliceLayout() == nullptr);
        EXPECT_EQ(solver->getSliceLayout()->getDomain()[0], 8);
        EXPECT_EQ(solver->getSliceLayout()->getDomain()[1], 9);
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_NoParticles) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        ASSERT_NO_THROW(rebuildBunch());
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        EXPECT_NO_THROW(solver->scatterToGrid(*bunch));
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_NoReferencePath) {
        fsCmd->setType("FFT2D5");
        ASSERT_NO_THROW(rebuildBunch());
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{1, 2, 3}}, {{4, 5, 6}});
        const FrenetSerretScatterDiagnostic info(pc.get());
        solver->scatterToGrid<FrenetSerretScatterDiagnostic>(*bunch, info);
        auto [r, p] = getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {1, 2, 3}, {4, 5, 6});
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_ShortReferencePath) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}});
        fsCmd->setType("FFT2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{1, 2, 3}}, {{4, 5, 6}});
        solver->loadReferencePath();
        const FrenetSerretScatterDiagnostic info(pc.get());
        solver->scatterToGrid<FrenetSerretScatterDiagnostic>(*bunch, info);
        auto [r, p] = getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {1, 2, 3}, {4, 5, 6});
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_TrivialReferencePath) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}});
        fsCmd->setType("FFT2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 0}, {0, 0, 0.5}, {0, 0, 1}}, {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}});
        solver->loadReferencePath();
        const FrenetSerretScatterDiagnostic info(pc.get());
        solver->scatterToGrid<FrenetSerretScatterDiagnostic>(*bunch, info);
        auto [r, p] = info.getParticles();
        ASSERT_EQ(r.size(), 3);
        ASSERT_EQ(p.size(), 3);
        expectParticle(0, r, p, {0, 0, 0}, {0, 0, 1});
        expectParticle(1, r, p, {0, 0, 0.5}, {0, 0, 1});
        expectParticle(2, r, p, {0, 0, 1}, {0, 0, 1});
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}});
        fsCmd->setType("FFT2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 2}}, {{0, 0, 1}});
        solver->loadReferencePath();
        const FrenetSerretScatterDiagnostic info(pc.get());
        solver->scatterToGrid<FrenetSerretScatterDiagnostic>(*bunch, info);
        auto [r, p] = info.getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {-0.7071068, 0, 1.7071068}, {-0.7071068, 0, 0.7071068});
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_AtRefCorner) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{2, 0, 2}}, {{0, 0, 1}});
        solver->loadReferencePath();
        const FrenetSerretScatterDiagnostic info(pc.get());
        solver->scatterToGrid<FrenetSerretScatterDiagnostic>(*bunch, info);
        auto [r, p] = info.getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {0.7071068, 0, 2.4142136}, {-0.7071068, 0, 0.7071068});
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_DegenerateRefPath) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 0}, {0, 0, 1}, {1, 0, 2}});
        fsCmd->setType("FFT2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 2}}, {{0, 0, 1}});
        solver->loadReferencePath();
        const FrenetSerretScatterDiagnostic info(pc.get());
        solver->scatterToGrid<FrenetSerretScatterDiagnostic>(*bunch, info);
        auto [r, p] = info.getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {-0.7071068, 0, 1.7071068}, {-0.7071068, 0, 0.7071068});
    }

    TEST_F(TestSolve2d5, BoostToBeamFrame_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0.5, 2}}, {{0.001, 0.002, 0.577}});
        solver->loadReferencePath();
        const BoostToBeamDiagnostic info(pc.get());
        solver->scatterToGrid<BoostToBeamDiagnostic>(*bunch, info);
        auto [r, p] = info.getParticles();
        ASSERT_EQ(r.size(), 1);
        expectParticle(0, r, p, {0, 0.5, 2.0}, {0.001, 0.002, 0.0});
    }

    TEST_F(TestSolve2d5, ScatterToGrid_SimpleCharge) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 0}}, {{0, 0, 0}});
        solver->loadReferencePath();
        const ScatterChargeDiagnostic info(solver->getRho()->getView());
        solver->scatterToGrid<ScatterChargeDiagnostic>(*bunch, info);
        expectChargeDensity(
                info.rhoView, {{6, 6, 6, 0.00520833},
                               {6, 6, 7, 0.00520833},
                               {6, 7, 6, 0.00520833},
                               {6, 7, 7, 0.00520833},
                               {7, 6, 6, 0.00520833},
                               {7, 6, 7, 0.00520833},
                               {7, 7, 6, 0.00520833},
                               {7, 7, 7, 0.00520833},
                               {8, 7, 7, 0.00000}});
    }

    TEST_F(TestSolve2d5, ScatterToGrid_ClosedRingCharge) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, -3}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        fsCmd->setClosedRing(true);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, -2}}, {{0, 0, 0}});
        solver->loadReferencePath();
        const ScatterChargeDiagnostic info(solver->getRho()->getView());
        solver->scatterToGrid<ScatterChargeDiagnostic>(*bunch, info);
        expectChargeDensity(
                info.rhoView, {{6, 6, 6, 0.00520833},
                               {6, 6, 7, 0.00520833},
                               {6, 7, 6, 0.00520833},
                               {6, 7, 7, 0.00520833},
                               {7, 6, 6, 0.00520833},
                               {7, 6, 7, 0.00520833},
                               {7, 7, 6, 0.00520833},
                               {7, 7, 7, 0.00520833},
                               {8, 7, 7, 0.00000}});
    }

    TEST_F(TestSolve2d5, ScatterToGrid_SimpleChargeDensity) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 0}}, {{0, 0, 0}});
        solver->loadReferencePath();
        const ScatterChargeDensityDiagnostic info(solver->getRho()->getView());
        solver->scatterToGrid<ScatterChargeDensityDiagnostic>(*bunch, info);
        expectChargeDensity(
                info.rhoView, {{6, 6, 6, 1.00000},
                               {6, 6, 7, 1.00000},
                               {6, 7, 6, 1.00000},
                               {6, 7, 7, 1.00000},
                               {7, 6, 6, 1.00000},
                               {7, 6, 7, 1.00000},
                               {7, 7, 6, 1.00000},
                               {7, 7, 7, 1.00000},
                               {8, 7, 7, 0.00000}});
    }

    TEST_F(TestSolve2d5, TotalDensity_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 0}}, {{0, 0, 0}});
        solver->loadReferencePath();
        solver->scatterToGrid(*bunch);
        const TotalDensityDiagnostic totalInfo(solver->getLineDensity());
        solver->calculateLineDensity<TotalDensityDiagnostic>(totalInfo);
        expectLineDensity(
                totalInfo.lineDensityView_m, {0, 0, 0, 0, 0, 0, 0, 4, 4, 0, 0, 0, 0, 0});
    }

    TEST_F(TestSolve2d5, LineDensity_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 0}}, {{0, 0, 0}});
        solver->loadReferencePath();
        solver->scatterToGrid(*bunch);
        const LineDensityDiagnostic lineInfo(solver->getLineDensity());
        solver->calculateLineDensity<LineDensityDiagnostic>(lineInfo);
        expectLineDensity(
                lineInfo.lineDensityView_m, {0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0});
    }

    TEST_F(TestSolve2d5, LineDensityGradient_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 0}}, {{0, 0, 0}});
        solver->loadReferencePath();
        solver->scatterToGrid(*bunch);
        const LineDensityGradientDiagnostic lineInfo(solver->getLineDensityGradient());
        solver->calculateLineDensity<LineDensityGradientDiagnostic>(lineInfo);
        expectLineDensity(
                lineInfo.lineDensityView_m,
                {0, 0, 0, 0, 0, 1, 1, -1, -1, 0, 0, 0});
    }

    TEST_F(TestSolve2d5, LineDensity_RingNotClosed) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        fsCmd->setClosedRing(false);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 2}}, {{0, 0, 0}});
        solver->loadReferencePath();
        solver->scatterToGrid(*bunch);
        const LineDensityDiagnostic lineInfo(solver->getLineDensity());
        solver->calculateLineDensity<LineDensityDiagnostic>(lineInfo);
        expectLineDensity(
                lineInfo.lineDensityView_m, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0});
    }

    TEST_F(TestSolve2d5, LineDensity_RingClosed) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        fsCmd->setClosedRing(true);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 2}}, {{0, 0, 0}});
        solver->loadReferencePath();
        solver->scatterToGrid(*bunch);
        const LineDensityDiagnostic lineInfo(solver->getLineDensity());
        solver->calculateLineDensity<LineDensityDiagnostic>(lineInfo);
        expectLineDensity(
                lineInfo.lineDensityView_m, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0});
    }

    TEST_F(TestSolve2d5, SolvePoissons_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 0}}, {{0, 0, 0}});
        solver->loadReferencePath();
        solver->scatterToGrid(*bunch);
        const EFieldDiagnostic info(solver->getEField()->getView());
        solver->solvePoissons<EFieldDiagnostic>(info);
        expectEField(
                info.eFieldView_m, {{1, 1, 6, -3.250090326e9, -3.250090326e9, 0},
                                    {2, 1, 6, -3.218709596e9, -3.896406064e9, 0},
                                    {4, 4, 6, -7.214724490e9, -7.214724490e9, 0},
                                    {4, 8, 6, -10.565567441e9, 6.311498630e9, 0}});
    }

    TEST_F(TestSolve2d5, ToFrenetSerretGather_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 2}}, {{0, 0, 1}});
        solver->loadReferencePath();
        solver->scatterToGrid(*bunch);
        solver->solvePoissons();
        const FrenetSerretGatherDiagnostic info(pc.get());
        solver->gatherFromGrid<FrenetSerretGatherDiagnostic>(*bunch, info);
        auto [r, p] = info.getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {-0.7071068, 0, 1.7071068}, {-0.7071068, 0, 0.7071068});
    }

    TEST_F(TestSolve2d5, GatherEField_TwoParticles) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{2, 0, 0}, {-2, 0, 0}}, {{0, 0, 0}, {0, 0, 0}});
        solver->loadReferencePath();
        solver->scatterToGrid(*bunch);
        solver->solvePoissons();
        const GatherEFieldDiagnostic gatherInfo(pc.get());
        solver->gatherFromGrid<GatherEFieldDiagnostic>(*bunch, gatherInfo);
        auto [e, b] = gatherInfo.getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(0, e, b, {4.493908875e9, 0, 0}, {0, 0, 0}, 1e3);
        expectParticleFields(1, e, b, {-4.493908875e9, 0, 0}, {0, 0, 0}, 1e3);
    }

    TEST_F(TestSolve2d5, Deboost_TwoStationaryParticles) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{2, 0, 0}, {-2, 0, 0}}, {{0, 0, 0}, {0, 0, 0}});
        solver->loadReferencePath();
        solver->scatterToGrid(*bunch);
        solver->solvePoissons();
        const DeboostedDiagnostic info(pc.get());
        solver->gatherFromGrid<DeboostedDiagnostic>(*bunch, info);
        auto [e, b] = info.getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(0, e, b, {4.493908875e9, 0, 0}, {0, 0, 0}, 1e3);
        expectParticleFields(1, e, b, {-4.493908875e9, 0, 0}, {0, 0, 0}, 1e3);
    }

    TEST_F(TestSolve2d5, Deboost_TwoRelativisticParticles) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{2, 0, 0}, {-2, 0, 0}}, {{0, 0, 1}, {0, 0, 1}});
        solver->loadReferencePath();
        solver->scatterToGrid(*bunch);
        solver->solvePoissons();
        const DeboostedDiagnostic info(pc.get());
        solver->gatherFromGrid<DeboostedDiagnostic>(*bunch, info);
        auto [e, b] = info.getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(0, e, b, {6.355346880e9, 0, 0}, {0, -14.9900, 0}, 1e3, 1e-4);
        expectParticleFields(1, e, b, {-6.355346880e9, 0, 0}, {0, 14.9900, 0}, 1e3, 1e-4);
    }

    TEST_F(TestSolve2d5, LongitudinalField_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{2, 0, 0}, {-2, 0, 0}}, {{0, 0, 1}, {0, 0, 1}});
        solver->loadReferencePath();
        solver->scatterToGrid(*bunch);
        solver->solvePoissons();
        const LineDensityGradientDiagnostic lineInfo(solver->getLineDensityGradient());
        solver->calculateLineDensity<LineDensityGradientDiagnostic>(lineInfo);
        expectLineDensity(lineInfo.lineDensityView_m,{});
        const LongitudinalFieldDiagnostic info(pc.get());
        solver->gatherFromGrid<LongitudinalFieldDiagnostic>(*bunch, info);
        auto [e, b] = info.getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(0, e, b, {6.355346880e9, 0, 57.160829398e9}, {0, -14.9900, 0}, 1e3, 1e-4);
        expectParticleFields(1, e, b, {-6.355346880e9, 0, 57.160829398e9}, {0, 14.9900, 0}, 1e3, 1e-4);
    }

    TEST_F(TestSolve2d5, LabFrameFields_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("FFT2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{2, 0, 0}, {-2, 0, 0}}, {{0, 0, 1}, {0, 0, 1}});
        solver->loadReferencePath();
        solver->scatterToGrid(*bunch);
        solver->solvePoissons();
        solver->calculateLineDensity();
        const LabFrameFieldsDiagnostic info(pc.get());
        solver->gatherFromGrid<LabFrameFieldsDiagnostic>(*bunch, info);
        auto [e, b] = info.getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(0, e, b, {6.355346880e9, 0, 57.160829398e9}, {0, -14.9900, 0}, 1e3, 1e-4);
        expectParticleFields(1, e, b, {-6.355346880e9, 0, 57.160829398e9}, {0, 14.9900, 0}, 1e3, 1e-4);
    }
}  // namespace
