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
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

namespace {

    constexpr bool VerboseTest = false;

    class TestableFieldSolverCmd : public FieldSolverCmd {
    public:
        void setType(const std::string& t) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::TYPE], t);
        }

        void setBinsName(const std::string& name) {
            Attributes::setString(this->itsAttr[FIELDSOLVER::BINS], name);
        }

        void setParallelDimensions(bool parallelX, bool parallelY, bool parallelZ) {
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTX], parallelX);
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTY], parallelY);
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTZ], parallelZ);
        }
    };

    using PartBunch_t         = PartBunch<double, 3>;
    using Solve2d5_t          = Solve2d5<double>;
    using ParticleContainer_t = PartBunch_t::ParticleContainer_t;
    using Point_t             = ippl::Vector<double, 3>;

    class Info : public Solve2d5_t::NullDiagnostic {
    public:
        explicit Info(const Kind kind) : kind_m(kind) {}

        void initialise(
                const PartBunch_t& partBunch, const Field_t<3>& rho,
                const Solve2d5<double>::LineDensityView_t& lineDensity,
                const Solve2d5<double>::LineDensityView_t& lineDensityGradient,
                const VField_t<T, 3>& eField) {
            auto& pcs = partBunch.getParticleContainers();
            if (!pcs.empty()) {
                auto& pc  = pcs.front();
                r_m       = Solve2d5_t::VectorView_t("fsr", pc->R.getParticleCount());
                p_m       = Solve2d5_t::VectorView_t("fsp", pc->R.getParticleCount());
                e_m       = Solve2d5_t::VectorView_t("e", pc->R.getParticleCount());
                b_m       = Solve2d5_t::VectorView_t("b", pc->R.getParticleCount());
                invalid_m = Solve2d5_t::BooleanView_t("fsinv", pc->R.getParticleCount());
            }
            Kokkos::resize(
                    rhoView_m, rho.getView().extent(0), rho.getView().extent(1),
                    rho.getView().extent(2));
            Kokkos::resize(lineDensityView_m, lineDensity.extent(0));
            Kokkos::resize(lineDensityGradientView_m, lineDensityGradient.extent(0));
            Kokkos::resize(
                    eFieldView_m, eField.getView().extent(0), eField.getView().extent(1),
                    eField.getView().extent(2));
        }

        KOKKOS_FUNCTION void frenetSerretScatter(
                const size_t n, const Solve2d5<double>::Vector3D_t& r,
                const Solve2d5<double>::Vector3D_t& p, const bool invalid) const {
            if (kind_m == Kind::FrenetSerretScatter) {
                r_m[n]       = r;
                p_m[n]       = p;
                invalid_m[n] = invalid;
            }
        }

        KOKKOS_FUNCTION void boostToBeam(
                const size_t n, const Solve2d5<double>::Vector3D_t& r,
                const Solve2d5<double>::Vector3D_t& p, const bool invalid) const {
            if (kind_m == Kind::BoostToBeam) {
                r_m[n]       = r;
                p_m[n]       = p;
                invalid_m[n] = invalid;
            }
        }

        KOKKOS_FUNCTION void scatterCharge(const Solve2d5_t::ScalarGridView3D_t& rho) const {
            if (kind_m == Kind::ScatterCharge) {
                Kokkos::deep_copy(rhoView_m, rho);
            }
        }

        KOKKOS_FUNCTION void scatterChargeDensity(const Solve2d5_t::ScalarGridView3D_t& rho) const {
            if (kind_m == Kind::ScatterChargeDensity) {
                Kokkos::deep_copy(rhoView_m, rho);
            }
        }

        KOKKOS_FUNCTION void totalDensity(const Solve2d5_t::LineDensityView_t& lineDensity) const {
            if (kind_m == Kind::TotalDensity) {
                Kokkos::deep_copy(lineDensityView_m, lineDensity);
            }
        }

        KOKKOS_FUNCTION void lineDensity(const Solve2d5_t::LineDensityView_t& lineDensity) const {
            if (kind_m == Kind::LineDensity) {
                Kokkos::deep_copy(lineDensityView_m, lineDensity);
            }
        }

        KOKKOS_FUNCTION void lineDensityGradient(
                const Solve2d5_t::LineDensityView_t& lineDensity) const {
            if (kind_m == Kind::LineDensityGradient) {
                Kokkos::deep_copy(lineDensityGradientView_m, lineDensity);
            }
        }

        KOKKOS_FUNCTION void potential(const Field_t<2U>::view_type& phi, const size_t z) const {
            if (kind_m == Kind::Potential) {
                const auto destSlice = Kokkos::subview(rhoView_m, Kokkos::ALL(), Kokkos::ALL(), z);
                Kokkos::deep_copy(destSlice, phi);
            }
        }

        KOKKOS_FUNCTION void eField(const VField_t<T, 3>::view_type& eField) const {
            if (kind_m == Kind::EField) {
                Kokkos::deep_copy(eFieldView_m, eField);
            }
        }

        KOKKOS_FUNCTION void frenetSerretGather(
                const size_t n, const Solve2d5<double>::Vector3D_t& r,
                const Solve2d5<double>::Vector3D_t& p, const bool invalid) const {
            if (kind_m == Kind::FrenetSerretGather) {
                r_m[n]       = r;
                p_m[n]       = p;
                invalid_m[n] = invalid;
            }
        }

        KOKKOS_FUNCTION void gatherEField(
                const size_t n, const Solve2d5<double>::Vector3D_t& e,
                const Solve2d5<double>::Vector3D_t& b, const bool invalid) const {
            if (kind_m == Kind::GatherEField) {
                e_m[n]       = e;
                b_m[n]       = b;
                invalid_m[n] = invalid;
            }
        }

        KOKKOS_FUNCTION void deboostFromBeam(
                const size_t n, const Solve2d5<double>::Vector3D_t& e,
                const Solve2d5<double>::Vector3D_t& b, const bool invalid) const {
            if (kind_m == Kind::Deboosted) {
                e_m[n]       = e;
                b_m[n]       = b;
                invalid_m[n] = invalid;
            }
        }

        KOKKOS_FUNCTION void longitudinalField(
                const size_t n, const Solve2d5<double>::Vector3D_t& e,
                const Solve2d5<double>::Vector3D_t& b, const bool invalid) const {
            if (kind_m == Kind::LongitudinalField) {
                e_m[n]       = e;
                b_m[n]       = b;
                invalid_m[n] = invalid;
            }
        }

        KOKKOS_FUNCTION void labFrameFields(
                const size_t n, const Solve2d5<double>::Vector3D_t& e,
                const Solve2d5<double>::Vector3D_t& b, const bool invalid) const {
            if (kind_m == Kind::LabFrameFields) {
                e_m[n]       = e;
                b_m[n]       = b;
                invalid_m[n] = invalid;
            }
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

        Field_t<3>::view_type rhoView_m;
        Solve2d5_t::VectorView_t r_m;
        Solve2d5_t::VectorView_t p_m;
        Solve2d5_t::VectorView_t e_m;
        Solve2d5_t::VectorView_t b_m;
        Solve2d5_t::BooleanView_t invalid_m;
        Solve2d5_t::LineDensityView_t lineDensityView_m;
        Solve2d5_t::LineDensityView_t lineDensityGradientView_m;
        VField_t<T, 3>::view_type eFieldView_m;
        Kind kind_m;
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
            fsCmd_m = std::make_shared<TestableFieldSolverCmd>();
            fsCmd_m->setType("NONE");
            fsCmd_m->setNX(nx);
            fsCmd_m->setNY(ny);
            fsCmd_m->setNZ(nz);
            fsCmd_m->setParallelDimensions(true, true, true);
            fsCmd_m->setFieldSolverCmdType();
            fsCmd_m->setDomainDecomposition();

            dataSink_m     = std::make_shared<DataSink>();
            beam_m         = std::make_shared<Beam>();
            Beam* testBeam = Beam::find("UNNAMED_BEAM");

            // qi/mi/lbt are used by rho scaling; but with NullSolver we mostly validate
            // "no-throw" and deterministic zero E behavior.
            bunch_m = std::make_shared<PartBunch_t>(
                    /*qi=*/std::vector{1.0},
                    /*mi=*/std::vector{1.0},
                    /*beams=*/std::vector{testBeam},
                    /*totalParticlesPerBeam=*/std::vector{kDefaultNParticles},
                    /*lbt=*/1.0,
                    /*integration_method=*/"LF2", fsCmd_m.get(), dataSink_m.get());
            pc_m = bunch_m->getParticleContainer();
        }

        void TearDown() override {
            // Ensure device allocations can be freed between tests.
            bunch_m.reset();
            dataSink_m.reset();
            fsCmd_m.reset();
            pc_m.reset();
            // Remove reference path file
            if (std::filesystem::exists("data/unit_test_DesignPath.dat")) {
                std::filesystem::remove("data/unit_test_DesignPath.dat");
            }
        }

        void createParticles(
                const std::vector<Vector_t<double, 3>>& r,
                const std::vector<Vector_t<double, 3>>& p,
                const std::vector<bool>& invalid = {}) const {
            pc_m->createParticles(r.size());
            const auto R_host       = pc_m->R.getHostMirror();
            const auto P_host       = pc_m->P.getHostMirror();
            const auto dt_host      = pc_m->dt.getHostMirror();
            const auto E_host       = pc_m->E.getHostMirror();
            const auto invalid_host = pc_m->InvalidMask.getHostMirror();
            const double dt         = bunch_m->getdT();
            const double qi         = pc_m->getChargePerParticle();
            for (size_t i = 0; i < r.size(); ++i) {
                R_host(i)  = r[i];
                P_host(i)  = p[i];
                dt_host(i) = dt;
                E_host(i)  = {0.0, 0.0, 0.0};
                if (i < invalid.size()) {
                    invalid_host(i) = invalid[i];
                } else {
                    invalid_host(i) = false;
                }
            }
            for (size_t i = r.size(); i < invalid_host.extent(0); ++i) {
                invalid_host(i) = true;
            }
            Kokkos::deep_copy(pc_m->R.getView(), R_host);
            Kokkos::deep_copy(pc_m->P.getView(), P_host);
            Kokkos::deep_copy(pc_m->dt.getView(), dt_host);
            Kokkos::deep_copy(pc_m->E.getView(), E_host);
            Kokkos::deep_copy(pc_m->InvalidMask.getView(), invalid_host);
            pc_m->setQ(qi);
            ippl::Comm->barrier();
            Kokkos::fence();
        }

        [[nodiscard]] std::tuple<std::vector<Vector_t<double, 3>>, std::vector<Vector_t<double, 3>>>
        getParticles() const {
            const auto R_host       = pc_m->R.getHostMirror();
            const auto P_host       = pc_m->P.getHostMirror();
            const auto invalid_host = pc_m->InvalidMask.getHostMirror();
            Kokkos::deep_copy(R_host, pc_m->R.getView());
            Kokkos::deep_copy(P_host, pc_m->P.getView());
            Kokkos::deep_copy(invalid_host, pc_m->InvalidMask.getView());
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
            bunch_m        = std::make_shared<PartBunch_t>(
                    /*qi=*/std::vector{1.0},
                    /*mi=*/std::vector{1.0},
                    /*beams=*/std::vector{testBeam},
                    /*totalParticlesPerBeam=*/std::vector{kDefaultNParticles},
                    /*lbt=*/1.0,
                    /*integration_method=*/"LF2", fsCmd_m.get(), dataSink_m.get());
            pc_m = bunch_m->getParticleContainer();
            bunch_m->getFieldSolver()->orbitThreadersReady();
        }

        void makeReferencePathFile(
                const std::filesystem::path& fileName, const std::vector<Point_t>& points,
                const bool noTimeColumn = false) {
            // Emulates the file created by OrbitThreader.  Use noTimeColumn = true to create
            // an invalid file.
            // Make the directory path exist
            auto parentPath = fileName.parent_path();
            if (!parentPath.empty()) {
                std::filesystem::create_directories(parentPath);
            }
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
            SCOPED_TRACE("Index = " + std::to_string(index));
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
            SCOPED_TRACE("Index = " + std::to_string(index));
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
                SCOPED_TRACE(
                        "Index: " + std::to_string(e.i) + "," + std::to_string(e.j) + ","
                        + std::to_string(e.k));
                EXPECT_NEAR(hostView(e.i, e.j, e.k), e.value, 1e-6);
            }
        }

        static void expectTotalCharge(
                const Solve2d5_t::ScalarGridView3D_t& rho, const double dt,
                const std::vector<double>& expectedCharge, const double tolerance) {
            // Transfer to host
            const auto hostView = Kokkos::create_mirror_view(rho);
            Kokkos::deep_copy(hostView, rho);
            // Print
            if (VerboseTest) {
                for (size_t k = 0; k < rho.extent(2); ++k) {
                    for (size_t j = 0; j < rho.extent(1); ++j) {
                        std::cout << k << "," << j << ": ";
                        for (size_t i = 0; i < rho.extent(0); ++i) {
                            std::cout << hostView(i, j, k) / dt << " ";
                        }
                        std::cout << std::endl;
                    }
                }
            }
            // Check values
            for (size_t slice = 0; slice < expectedCharge.size(); ++slice) {
                // Total charge in this slice
                double total{};
                for (size_t j = 0; j < rho.extent(1); ++j) {
                    for (size_t i = 0; i < rho.extent(0); ++i) {
                        total += hostView(i, j, slice);
                    }
                }
                // Was this expected?
                // Note that the charge reported is multiplied by the timestep dt.
                SCOPED_TRACE("Slice: " + std::to_string(slice));
                EXPECT_NEAR(total / dt, expectedCharge[slice], tolerance);
            }
        }

        static void expectTotalChargeDensity(
                const Solve2d5_t::ScalarGridView3D_t& rho,
                const std::vector<double>& expectedChargeDensity, const double tolerance) {
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
            for (size_t slice = 0; slice < expectedChargeDensity.size(); ++slice) {
                // Total charge density in this slice
                double total{};
                for (size_t j = 0; j < rho.extent(1); ++j) {
                    for (size_t i = 0; i < rho.extent(0); ++i) {
                        total += hostView(i, j, slice);
                    }
                }
                // Was this expected?
                // Note that the charge density reported is multiplied by the timestep dt.
                SCOPED_TRACE("Slice: " + std::to_string(slice));
                EXPECT_NEAR(total, expectedChargeDensity[slice], tolerance);
            }
        }

        struct PhiValue {
            size_t i;
            size_t j;
            size_t k;
            double value;
        };
        static void expectPotential(
                const Solve2d5_t::ScalarGridView3D_t& phi, const std::vector<PhiValue>& expected,
                const double error = 1e-6) {
            // Transfer to host
            const auto hostView = Kokkos::create_mirror_view(phi);
            Kokkos::deep_copy(hostView, phi);
            // Print
            if (VerboseTest) {
                for (size_t k = 0; k < phi.extent(2); ++k) {
                    for (size_t j = 0; j < phi.extent(1); ++j) {
                        std::cout << k << "," << j << ": ";
                        for (size_t i = 0; i < phi.extent(0); ++i) {
                            std::cout << hostView(i, j, k) << " ";
                        }
                        std::cout << std::endl;
                    }
                }
            }
            // Check values
            for (const auto& e : expected) {
                SCOPED_TRACE(
                        "Index: " + std::to_string(e.i) + "," + std::to_string(e.j) + ","
                        + std::to_string(e.k));
                EXPECT_NEAR(hostView(e.i, e.j, e.k), e.value, error);
            }
        }

        static void expectLineDensity(
                const Solve2d5_t::LineDensityView_t& lineDensity,
                const std::vector<double>& expected) {
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
                    SCOPED_TRACE("Index: " + std::to_string(i));
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
                std::cout << "Ex" << std::endl;
                for (size_t k = 0; k < hostView.extent(2); ++k) {
                    for (size_t j = 0; j < hostView.extent(1); ++j) {
                        std::cout << k << "," << j << ": ";
                        for (size_t i = 0; i < hostView.extent(0); ++i) {
                            std::cout << hostView(i, j, k).data_m[0] << " ";
                        }
                        std::cout << std::endl;
                    }
                }
                std::cout << "Ey" << std::endl;
                for (size_t k = 0; k < hostView.extent(2); ++k) {
                    for (size_t j = 0; j < hostView.extent(1); ++j) {
                        std::cout << k << "," << j << ": ";
                        for (size_t i = 0; i < hostView.extent(0); ++i) {
                            std::cout << hostView(i, j, k).data_m[1] << " ";
                        }
                        std::cout << std::endl;
                    }
                }
                std::cout << "Ez" << std::endl;
                for (size_t k = 0; k < hostView.extent(2); ++k) {
                    for (size_t j = 0; j < hostView.extent(1); ++j) {
                        std::cout << k << "," << j << ": ";
                        for (size_t i = 0; i < hostView.extent(0); ++i) {
                            std::cout << hostView(i, j, k).data_m[2] << " ";
                        }
                        std::cout << std::endl;
                    }
                }
            }
            // Check values
            for (const auto& e : expected) {
                SCOPED_TRACE(
                        "Index: " + std::to_string(e.i) + "," + std::to_string(e.j) + ","
                        + std::to_string(e.k));
                EXPECT_NEAR(hostView(e.i, e.j, e.k)[0], e.valueX, 1e3);
                EXPECT_NEAR(hostView(e.i, e.j, e.k)[1], e.valueY, 1e3);
                EXPECT_NEAR(hostView(e.i, e.j, e.k)[2], e.valueZ, 1e3);
            }
        }

        std::shared_ptr<TestableFieldSolverCmd> fsCmd_m;
        std::shared_ptr<DataSink> dataSink_m;
        std::shared_ptr<Beam> beam_m;
        std::shared_ptr<PartBunch_t> bunch_m;
        std::shared_ptr<ParticleContainer_t> pc_m;
    };

    TEST_F(TestSolve2d5, LoadReferencePath_Missing) {
        fsCmd_m->setType("FFT2D5");
        EXPECT_ANY_THROW(rebuildBunch());
    }

    TEST_F(TestSolve2d5, Configuration_SerialLayoutAccepted) {
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setParallelDimensions(false, false, false);

        ASSERT_NO_THROW(fsCmd_m->execute());
        const auto decomposition = fsCmd_m->getDomainDecomposition();
        EXPECT_FALSE(decomposition[0]);
        EXPECT_FALSE(decomposition[1]);
        EXPECT_FALSE(decomposition[2]);
    }

    TEST_F(TestSolve2d5, Configuration_ParallelLayoutRejected) {
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setParallelDimensions(false, false, true);

        EXPECT_THROW(fsCmd_m->execute(), OpalException);
    }

    TEST_F(TestSolve2d5, Configuration_BinningRejected) {
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setBinsName("UNUSED_BINNING");
        fsCmd_m->setParallelDimensions(false, false, false);

        EXPECT_THROW(fsCmd_m->execute(), OpalException);
    }

    TEST_F(TestSolve2d5, Configuration_MultipleRanksRejected) {
        if (ippl::Comm->size() == 1) {
            GTEST_SKIP() << "This validation requires multiple MPI ranks.";
        }

        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setParallelDimensions(false, false, false);
        try {
            fsCmd_m->execute();
            FAIL() << "Expected FFT2D5 input validation to reject multiple MPI ranks.";
        } catch (const OpalException& exception) {
            const std::string message = exception.what();
            EXPECT_NE(message.find("exactly one MPI rank"), std::string::npos);
        }

        try {
            rebuildBunch();
            FAIL() << "Expected FFT2D5 construction to reject multiple MPI ranks.";
        } catch (const OpalException& exception) {
            const std::string message = exception.what();
            EXPECT_NE(message.find("exactly one MPI rank"), std::string::npos);
        }
    }

    TEST_F(TestSolve2d5, LoadReferencePath_ReadFail) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}},
                true);
        fsCmd_m->setType("FFT2D5");
        EXPECT_ANY_THROW(rebuildBunch());
    }

    TEST_F(TestSolve2d5, LoadReferencePath_Success) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}});
        fsCmd_m->setType("FFT2D5");
        ASSERT_NO_THROW(rebuildBunch());
        auto* solver     = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        auto hostRefPath = Kokkos::create_mirror_view(solver->getReferencePath());
        Kokkos::deep_copy(hostRefPath, solver->getReferencePath());
        EXPECT_EQ(hostRefPath.size(), 4);
        EXPECT_EQ(hostRefPath[0].data_m[0], 0);
        EXPECT_EQ(hostRefPath[0].data_m[1], 0);
        EXPECT_EQ(hostRefPath[0].data_m[2], 0);
        EXPECT_EQ(hostRefPath[1].data_m[0], 0);
        EXPECT_EQ(hostRefPath[1].data_m[1], 0);
        EXPECT_EQ(hostRefPath[1].data_m[2], 1);
        EXPECT_EQ(hostRefPath[2].data_m[0], 1);
        EXPECT_EQ(hostRefPath[2].data_m[1], 0);
        EXPECT_EQ(hostRefPath[2].data_m[2], 2);
        EXPECT_EQ(hostRefPath[3].data_m[0], 0);
        EXPECT_EQ(hostRefPath[3].data_m[1], 0);
        EXPECT_EQ(hostRefPath[3].data_m[2], 3);
    }

    TEST_F(TestSolve2d5, LoadReferencePath_SpecifiedFile) {
        makeReferencePathFile(
                "Specified_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setRefPathFileName("Specified_DesignPath.dat");
        ASSERT_NO_THROW(rebuildBunch());
        auto* solver     = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        auto hostRefPath = Kokkos::create_mirror_view(solver->getReferencePath());
        Kokkos::deep_copy(hostRefPath, solver->getReferencePath());
        EXPECT_EQ(hostRefPath.size(), 4);
        EXPECT_EQ(hostRefPath[0].data_m[0], 0);
        EXPECT_EQ(hostRefPath[0].data_m[1], 0);
        EXPECT_EQ(hostRefPath[0].data_m[2], 0);
        EXPECT_EQ(hostRefPath[1].data_m[0], 0);
        EXPECT_EQ(hostRefPath[1].data_m[1], 0);
        EXPECT_EQ(hostRefPath[1].data_m[2], 1);
        EXPECT_EQ(hostRefPath[2].data_m[0], 1);
        EXPECT_EQ(hostRefPath[2].data_m[1], 0);
        EXPECT_EQ(hostRefPath[2].data_m[2], 2);
        EXPECT_EQ(hostRefPath[3].data_m[0], 0);
        EXPECT_EQ(hostRefPath[3].data_m[1], 0);
        EXPECT_EQ(hostRefPath[3].data_m[2], 3);
    }

    TEST_F(TestSolve2d5, SliceSolverSetup) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}});
        fsCmd_m->setType("FFT2D5");
        ASSERT_NO_THROW(rebuildBunch());
        const auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        EXPECT_EQ(solver->getNumSlices(), nz);
        ASSERT_FALSE(solver->getSliceMesh() == nullptr);
        EXPECT_EQ(solver->getSliceMesh()->getGridsize(0), 8);
        EXPECT_EQ(solver->getSliceMesh()->getGridsize(1), 9);
        EXPECT_EQ(solver->getSliceMesh()->getOrigin()[0], -0.5);
        EXPECT_EQ(solver->getSliceMesh()->getOrigin()[1], -0.5);
        EXPECT_DOUBLE_EQ(solver->getSliceMesh()->getMeshSpacing()[0], 0.125);
        EXPECT_NEAR(solver->getSliceMesh()->getMeshSpacing()[1], 0.111111, 1e-6);
        ASSERT_FALSE(solver->getSliceLayout() == nullptr);
        EXPECT_EQ(solver->getSliceLayout()->getDomain()[0], 8);
        EXPECT_EQ(solver->getSliceLayout()->getDomain()[1], 9);
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_NoParticles) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}});
        fsCmd_m->setType("FFT2D5");
        ASSERT_NO_THROW(rebuildBunch());
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        EXPECT_NO_THROW(solver->scatterToGrid(*bunch_m));
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_ShortReferencePath) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}});
        fsCmd_m->setType("FFT2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{1, 2, 3}}, {{4, 5, 6}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::FrenetSerretScatter);
        solver->scatterToGrid<Info>(*bunch_m, *info);
        auto [r, p] = getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {1, 2, 3}, {4, 5, 6});
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_TrivialReferencePath) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}});
        fsCmd_m->setType("FFT2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 0}, {0, 0, 0.5}, {0, 0, 1}}, {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::FrenetSerretScatter);
        solver->scatterToGrid<Info>(*bunch_m, *info);
        auto [r, p] = info->getParticles();
        ASSERT_EQ(r.size(), 3);
        ASSERT_EQ(p.size(), 3);
        expectParticle(0, r, p, {0, 0, 0}, {0, 0, 1});
        expectParticle(1, r, p, {0, 0, 0.5}, {0, 0, 1});
        expectParticle(2, r, p, {0, 0, 1}, {0, 0, 1});
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}});
        fsCmd_m->setType("FFT2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 2}}, {{0, 0, 1}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::FrenetSerretScatter);
        solver->scatterToGrid<Info>(*bunch_m, *info);
        auto [r, p] = info->getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {-0.7071068, 0, 1.7071068}, {-0.7071068, 0, 0.7071068});
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_AtRefCorner) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}});
        fsCmd_m->setType("FFT2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{2, 0, 2}}, {{0, 0, 1}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::FrenetSerretScatter);
        solver->scatterToGrid<Info>(*bunch_m, *info);
        auto [r, p] = info->getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {0.7071068, 0, 2.4142136}, {-0.7071068, 0, 0.7071068});
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_DegenerateRefPath) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 0}, {0, 0, 1}, {1, 0, 2}});
        fsCmd_m->setType("FFT2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 2}}, {{0, 0, 1}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::FrenetSerretScatter);
        solver->scatterToGrid<Info>(*bunch_m, *info);
        auto [r, p] = info->getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {-0.7071068, 0, 1.7071068}, {-0.7071068, 0, 0.7071068});
    }

    TEST_F(TestSolve2d5, BoostToBeamFrame_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd_m->setType("FFT2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0.5, 2}}, {{0.001, 0.002, 0.577}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::BoostToBeam);
        solver->scatterToGrid<Info>(*bunch_m, *info);
        auto [r, p] = info->getParticles();
        ASSERT_EQ(r.size(), 1);
        expectParticle(0, r, p, {0, 0.5, 2.0}, {0.001, 0.002, 0.0});
    }

    TEST_F(TestSolve2d5, BoostToBeamFrame_SimpleInvalid) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd_m->setType("FFT2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles(
                {{0, 0.5, 2}, {0, 0, 2}}, {{0.001, 0.002, 0.577}, {0.001, 0.002, 0.577}},
                {false, true});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::BoostToBeam);
        solver->scatterToGrid<Info>(*bunch_m, *info);
        auto [r, p] = info->getParticles();
        ASSERT_EQ(r.size(), 2);
        expectParticle(0, r, p, {0, 0.5, 2.0}, {0.001, 0.002, 0.0});
        expectParticle(1, r, p, {0, 0, 0}, {0, 0, 0});  // The invalid particle
    }

    TEST_F(TestSolve2d5, ScatterToGrid_SimpleCharge) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 3}}, {{0, 0, 0}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterCharge);
        solver->scatterToGrid<Info>(*bunch_m, *info);
        expectChargeDensity(
                info->rhoView_m, {{6, 6, 6, 0.00520833},
                                  {6, 6, 7, 0.00520833},
                                  {6, 7, 6, 0.00520833},
                                  {6, 7, 7, 0.00520833},
                                  {7, 6, 6, 0.00520833},
                                  {7, 6, 7, 0.00520833},
                                  {7, 7, 6, 0.00520833},
                                  {7, 7, 7, 0.00520833},
                                  {8, 7, 7, 0.00000}});
    }

    TEST_F(TestSolve2d5, ScatterToGrid_ClosedRingChargeLo) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        fsCmd_m->setClosedRing(true);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 0}}, {{0, 0, 0}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterCharge);
        solver->scatterToGrid<Info>(*bunch_m, *info);
        expectChargeDensity(
                info->rhoView_m, {{6, 6, 0, 0.00520833},
                                  {6, 6, 1, 0.00520833},
                                  {6, 7, 0, 0.00520833},
                                  {6, 7, 1, 0.00520833},
                                  {7, 6, 0, 0.00520833},
                                  {7, 6, 1, 0.00520833},
                                  {7, 7, 0, 0.00520833},
                                  {7, 7, 1, 0.00520833},
                                  {8, 7, 1, 0.00000}});
    }

    TEST_F(TestSolve2d5, ScatterToGrid_ClosedRingChargeHi) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        fsCmd_m->setClosedRing(true);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 6}}, {{0, 0, 0}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterCharge);
        solver->scatterToGrid<Info>(*bunch_m, *info);
        expectChargeDensity(
                info->rhoView_m, {{6, 6, 12, 0.00520833},
                                  {6, 6, 13, 0.00520833},
                                  {6, 7, 12, 0.00520833},
                                  {6, 7, 13, 0.00520833},
                                  {7, 6, 12, 0.00520833},
                                  {7, 6, 13, 0.00520833},
                                  {7, 7, 12, 0.00520833},
                                  {7, 7, 13, 0.00520833},
                                  {8, 7, 13, 0.00000}});
    }

    TEST_F(TestSolve2d5, ScatterToGrid_SimpleChargeDensity) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 3}}, {{0, 0, 0}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterChargeDensity);
        solver->scatterToGrid<Info>(*bunch_m, *info);
        expectChargeDensity(
                info->rhoView_m, {{6, 6, 6, 1.00000},
                                  {6, 6, 7, 1.00000},
                                  {6, 7, 6, 1.00000},
                                  {6, 7, 7, 1.00000},
                                  {7, 6, 6, 1.00000},
                                  {7, 6, 7, 1.00000},
                                  {7, 7, 6, 1.00000},
                                  {7, 7, 7, 1.00000},
                                  {8, 7, 7, 0.00000}});
    }

    TEST_F(TestSolve2d5, ScatterToGrid_ClosedRingChargeDensityLo) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        fsCmd_m->setClosedRing(true);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 0}}, {{0, 0, 0}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterChargeDensity);
        solver->scatterToGrid<Info>(*bunch_m, *info);
        expectChargeDensity(
                info->rhoView_m, {{6, 6, 12, 1.00000},
                                  {6, 6, 1, 1.00000},
                                  {6, 7, 12, 1.00000},
                                  {6, 7, 1, 1.00000},
                                  {7, 6, 12, 1.00000},
                                  {7, 6, 1, 1.00000},
                                  {7, 7, 12, 1.00000},
                                  {7, 7, 1, 1.00000},
                                  {8, 7, 1, 0.00000}});
    }

    TEST_F(TestSolve2d5, ScatterToGrid_ClosedRingChargeDensityHi) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        fsCmd_m->setClosedRing(true);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 6}}, {{0, 0, 0}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterChargeDensity);
        solver->scatterToGrid<Info>(*bunch_m, *info);
        expectChargeDensity(
                info->rhoView_m, {{6, 6, 12, 1.00000},
                                  {6, 6, 1, 1.00000},
                                  {6, 7, 12, 1.00000},
                                  {6, 7, 1, 1.00000},
                                  {7, 6, 12, 1.00000},
                                  {7, 6, 1, 1.00000},
                                  {7, 7, 12, 1.00000},
                                  {7, 7, 1, 1.00000},
                                  {8, 7, 1, 0.00000}});
    }

    TEST_F(TestSolve2d5, TotalDensity_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 3}}, {{0, 0, 0}});
        solver->scatterToGrid(*bunch_m);
        const auto totalInfo = solver->createDiagnostic<Info>(Info::Kind::TotalDensity);
        solver->calculateLineDensity<Info>(*totalInfo);
        expectLineDensity(totalInfo->lineDensityView_m, {0, 0, 0, 0, 0, 0, 4, 4, 0, 0, 0, 0, 0, 0});
    }

    TEST_F(TestSolve2d5, LineDensity_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 3}}, {{0, 0, 0}});
        solver->scatterToGrid(*bunch_m);
        const auto lineInfo = solver->createDiagnostic<Info>(Info::Kind::LineDensity);
        solver->calculateLineDensity<Info>(*lineInfo);
        expectLineDensity(lineInfo->lineDensityView_m, {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0});
    }

    TEST_F(TestSolve2d5, LineDensityGradient_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 3}}, {{0, 0, 0}});
        solver->scatterToGrid(*bunch_m);
        const auto lineInfo = solver->createDiagnostic<Info>(Info::Kind::LineDensityGradient);
        solver->calculateLineDensity<Info>(*lineInfo);
        expectLineDensity(
                lineInfo->lineDensityGradientView_m, {0, 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0});
    }

    TEST_F(TestSolve2d5, LineDensity_RingNotClosed) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        fsCmd_m->setClosedRing(false);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 5}}, {{0, 0, 0}});
        solver->scatterToGrid(*bunch_m);
        const auto lineInfo = solver->createDiagnostic<Info>(Info::Kind::LineDensity);
        solver->calculateLineDensity<Info>(*lineInfo);
        expectLineDensity(lineInfo->lineDensityView_m, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0});
    }

    TEST_F(TestSolve2d5, LineDensity_RingClosed) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        fsCmd_m->setClosedRing(true);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 5.5}}, {{0, 0, 0}});
        solver->scatterToGrid(*bunch_m);
        const auto lineInfo = solver->createDiagnostic<Info>(Info::Kind::LineDensity);
        solver->calculateLineDensity<Info>(*lineInfo);
        expectLineDensity(lineInfo->lineDensityView_m, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0});
    }

    TEST_F(TestSolve2d5, SolvePoissons_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 3}}, {{0, 0, 0}});
        solver->scatterToGrid(*bunch_m);
        const auto info = solver->createDiagnostic<Info>(Info::Kind::EField);
        solver->solvePoissons<Info>(*info);
        expectEField(
                info->eFieldView_m, {{1, 1, 6, -3.250090326e9, -3.250090326e9, 0},
                                     {2, 1, 6, -3.218709596e9, -3.896406064e9, 0},
                                     {4, 4, 6, -7.214724490e9, -7.214724490e9, 0},
                                     {4, 8, 6, -10.565567441e9, 6.311498630e9, 0}});
    }

    TEST_F(TestSolve2d5, ToFrenetSerretGather_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, 2}}, {{0, 0, 1}});
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::FrenetSerretGather);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [r, p] = info->getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {-0.7071068, 0, 1.7071068}, {-0.7071068, 0, 0.7071068});
    }

    TEST_F(TestSolve2d5, GatherEField_TwoParticles) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 0}, {0, 0, 0}});
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        const auto gatherInfo = solver->createDiagnostic<Info>(Info::Kind::GatherEField);
        solver->gatherFromGrid<Info>(*bunch_m, *gatherInfo);
        auto [e, b] = gatherInfo->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(0, e, b, {8.987817751e9, 0, 0}, {0, 0, 0}, 1e3);
        expectParticleFields(1, e, b, {-8.987817751e9, 0, 0}, {0, 0, 0}, 1e3);
    }

    TEST_F(TestSolve2d5, Deboost_TwoStationaryParticles) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 0}, {0, 0, 0}});
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::Deboosted);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(0, e, b, {8.987817751e9, 0, 0}, {0, 0, 0}, 1e3);
        expectParticleFields(1, e, b, {-8.987817751e9, 0, 0}, {0, 0, 0}, 1e3);
    }

    TEST_F(TestSolve2d5, Deboost_TwoRelativisticParticles) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::Deboosted);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(0, e, b, {12.710693760e9, 0, 0}, {0, -29.9801, 0}, 1e3, 1e-4);
        expectParticleFields(1, e, b, {-12.710693760e9, 0, 0}, {0, 29.9801, 0}, 1e3, 1e-4);
    }

    TEST_F(TestSolve2d5, LongitudinalField_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        const auto lineInfo = solver->createDiagnostic<Info>(Info::Kind::LineDensityGradient);
        solver->calculateLineDensity<Info>(*lineInfo);
        expectLineDensity(lineInfo->lineDensityGradientView_m, {});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LongitudinalField);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(
                0, e, b, {12.710693760e9, 0, 57.160829398e9}, {0, -29.9801, 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-12.710693760e9, 0, 57.160829398e9}, {0, 29.9801, 0}, 1e3, 1e-4);
    }

    TEST_F(TestSolve2d5, LabFrameFields_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(
                0, e, b, {12.710693760e9, 0, 57.160829398e9}, {0, -29.9801, 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-12.710693760e9, 0, 57.160829398e9}, {0, 29.9801, 0}, 1e3, 1e-4);
    }

    TEST_F(TestSolve2d5, LabFrameFields_SimpleInvalid) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles(
                {{2, 0, 3}, {-2, 0, 3}, {-2, 0, 3}}, {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}},
                {false, false, true});
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 3);
        ASSERT_EQ(b.size(), 3);
        expectParticleFields(
                0, e, b, {12.710693760e9, 0, 57.160829398e9}, {0, -29.9801, 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-12.710693760e9, 0, 57.160829398e9}, {0, 29.9801, 0}, 1e3, 1e-4);
        expectParticleFields(2, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
    }

    TEST_F(TestSolve2d5, LabFrameFields_SimpleOutOfBoundsX) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{4, 0, 3}, {-4, 0, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(0, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
        expectParticleFields(1, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
    }

    TEST_F(TestSolve2d5, LabFrameFields_SimpleOutOfBoundsY) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 4, 3}, {0, -4, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(0, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
        expectParticleFields(1, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
    }

    TEST_F(TestSolve2d5, LabFrameFields_SimpleOutOfBoundsZ) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{0, 0, -1}, {0, 0, 7}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(0, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
        expectParticleFields(1, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
    }

    TEST_F(TestSolve2d5, LabFrameFields_ClosedRing) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        fsCmd_m->setClosedRing(true);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{2, 0, 5}, {-2, 0, 5}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(*bunch_m);
        const auto efieldInfo = solver->createDiagnostic<Info>(Info::Kind::EField);
        solver->solvePoissons<Info>(*efieldInfo);
        expectEField(efieldInfo->eFieldView_m, {});
        const auto lineInfo = solver->createDiagnostic<Info>(Info::Kind::LineDensityGradient);
        solver->calculateLineDensity<Info>(*lineInfo);
        expectLineDensity(lineInfo->lineDensityGradientView_m, {});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 3);
        ASSERT_EQ(b.size(), 3);
        expectParticleFields(
                0, e, b, {12.710693760e9, 0, 28.580515699e9}, {0, -29.9801, 0}, 1e6, 1e-4);
        expectParticleFields(
                1, e, b, {-12.710693760e9, 0, 28.580515699e9}, {0, 29.9801, 0}, 1e6, 1e-4);
        expectParticleFields(2, e, b, {0, 0, 28.580515699e9}, {0, 0, 0}, 1e6, 1e-4);
    }

    TEST_F(TestSolve2d5, LabFrameFields_OpenRing) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        fsCmd_m->setClosedRing(false);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{2, 0, 5}, {-2, 0, 5}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(*bunch_m);
        const auto efieldInfo = solver->createDiagnostic<Info>(Info::Kind::EField);
        solver->solvePoissons<Info>(*efieldInfo);
        expectEField(efieldInfo->eFieldView_m, {});
        const auto lineInfo = solver->createDiagnostic<Info>(Info::Kind::LineDensityGradient);
        solver->calculateLineDensity<Info>(*lineInfo);
        expectLineDensity(lineInfo->lineDensityGradientView_m, {});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 3);
        ASSERT_EQ(b.size(), 3);
        expectParticleFields(
                0, e, b, {12.710693760e9, 0, 57.160829398e9}, {0, -29.9801, 0}, 1e6, 1e-4);
        expectParticleFields(
                1, e, b, {-12.710693760e9, 0, 57.160829398e9}, {0, 29.9801, 0}, 1e6, 1e-4);
        expectParticleFields(2, e, b, {0, 0, 0}, {0, 0, 0}, 1e6, 1e-4);
    }

    TEST_F(TestSolve2d5, LabFrameFields_Pipe) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        fsCmd_m->setBeamRadius(0.5);
        fsCmd_m->setPipeMode("CIRCULAR");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(
                0, e, b, {12.710693760e9, 0, 50.688114128e9}, {0, -29.9801, 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-12.710693760e9, 0, 50.688114128e9}, {0, 29.9801, 0}, 1e3, 1e-4);
    }

    TEST_F(TestSolve2d5, LabFrameFields_Plates) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        fsCmd_m->setBeamRadius(0.5);
        fsCmd_m->setPipeMode("PLATES");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(
                0, e, b, {12.710693760e9, 0, 55.030260593e9}, {0, -29.9801, 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-12.710693760e9, 0, 55.030260593e9}, {0, 29.9801, 0}, 1e3, 1e-4);
    }

    TEST_F(TestSolve2d5, LoadReferencePath_Empty) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        fsCmd_m->setClosedRing(false);
        ASSERT_NO_THROW(rebuildBunch());
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        EXPECT_EQ(solver->getReferencePath().size(), 0);
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(0, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
        expectParticleFields(1, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
    }

    TEST_F(TestSolve2d5, LoadReferencePath_EmptyClosedRing) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        fsCmd_m->setClosedRing(true);
        ASSERT_NO_THROW(rebuildBunch());
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        EXPECT_EQ(solver->getReferencePath().size(), 0);
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        solver->calculateLineDensity();
        fsCmd_m->setClosedRing(false);
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(0, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
        expectParticleFields(1, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
    }

    TEST_F(TestSolve2d5, LabFrameFields_NoParticles) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        solver->scatterToGrid(*bunch_m);
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(*bunch_m, *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 0);
        ASSERT_EQ(b.size(), 0);
    }

    TEST_F(TestSolve2d5, LabFrameFields_MainApi) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        fsCmd_m->setClosedRing(true);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles({{2, 0, 5}, {-2, 0, 5}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}});
        solver->runSolver();
        // Check the fields through the original particle bunch object
        auto& pc         = *bunch_m->getParticleContainers().front();
        const auto eHost = pc.E.getHostMirror();
        const auto bHost = pc.B.getHostMirror();
        Kokkos::deep_copy(eHost, pc.E.getView());
        Kokkos::deep_copy(bHost, pc.B.getView());
        EXPECT_NEAR(eHost(0).data_m[0], 12.710693760e9, 1e6);
        EXPECT_NEAR(eHost(0).data_m[1], 0, 1e6);
        EXPECT_NEAR(eHost(0).data_m[2], 28.580515699e9, 1e6);
        EXPECT_NEAR(eHost(1).data_m[0], -12.710693760e9, 1e6);
        EXPECT_NEAR(eHost(1).data_m[1], 0, 1e6);
        EXPECT_NEAR(eHost(1).data_m[2], 28.580515699e9, 1e6);
        EXPECT_NEAR(eHost(2).data_m[0], 0, 1e6);
        EXPECT_NEAR(eHost(2).data_m[1], 0, 1e6);
        EXPECT_NEAR(eHost(2).data_m[2], 28.580515699e9, 1e6);
        EXPECT_NEAR(bHost(0).data_m[0], 0, 1e-4);
        EXPECT_NEAR(bHost(0).data_m[1], -29.9801, 1e-4);
        EXPECT_NEAR(bHost(0).data_m[2], 0, 1e-4);
        EXPECT_NEAR(bHost(1).data_m[0], 0, 1e-4);
        EXPECT_NEAR(bHost(1).data_m[1], 29.9801, 1e-4);
        EXPECT_NEAR(bHost(1).data_m[2], 0, 1e-4);
        EXPECT_NEAR(bHost(2).data_m[0], 0, 1e-4);
        EXPECT_NEAR(bHost(2).data_m[1], 0, 1e-4);
        EXPECT_NEAR(bHost(2).data_m[2], 0, 1e-4);
    }

    const std::vector<Vector_t<double, 3>> kvR{
            {0, 0, 0.015},
            {0.001297, 0.000313, 0.015},
            {-0.002474, -0.000306, 0.015},
            {0.003779, -0.000992, 0.015},
            {0.003249, 0.002463, 0.015},
            {0.000431, 0.003054, 0.015},
            {0.000922, -0.000395, 0.015},
            {-0.00152, -0.001949, 0.015},
            {-0.000457, -0.001906, 0.015},
            {0.001279, -0.001958, 0.015},
            {-0.000838, 0.001237, 0.015},
            {0.00247, -0.000702, 0.015},
            {0.002301, 0.002154, 0.015},
            {0.004221, -0.002662, 0.015},
            {-0.000098, -0.004888, 0.015},
            {-0.003878, 0.003027, 0.015},
            {0.00359, 0.002568, 0.015},
            {0.000592, -0.001122, 0.015},
            {0.002341, -0.00093, 0.015},
            {0.00347, 0.003078, 0.015},
            {0.003485, 0.001949, 0.015},
            {0.001053, 0.004782, 0.015},
            {-0.001846, -0.002673, 0.015},
            {-0.004112, 0.00095, 0.015},
            {-0.003921, 0.000305, 0.015},
            {-0.001103, 0.002645, 0.015},
            {-0.004743, -0.001235, 0.015},
            {0.001427, -0.001055, 0.015},
            {0.000732, -0.000865, 0.015},
            {0.002362, -0.000268, 0.015},
            {-0.002485, 0.00323, 0.015},
            {0, -0.003618, 0.015},
            {0.000772, -0.001941, 0.015},
            {-0.000304, 0.001598, 0.015},
            {-0.000154, 0.001234, 0.015},
            {0.001468, -0.000561, 0.015},
            {0.004444, -0.000237, 0.015},
            {-0.001599, 0.003921, 0.015},
            {-0.004281, 0.001765, 0.015},
            {-0.00027, -0.000969, 0.015},
            {0.002707, 0.00219, 0.015},
            {0.002008, -0.001695, 0.015},
            {0.002548, -0.002384, 0.015},
            {0.003123, 0.001459, 0.015},
            {0.002301, -0.00299, 0.015},
            {0.000115, 0.002636, 0.015},
            {0.000956, 0.001245, 0.015},
            {0.002491, 0.004326, 0.015},
            {-0.000022, -0.002197, 0.015},
            {0.000936, -0.00379, 0.015},
            {-0.001211, 0.001373, 0.015},
            {-0.001251, -0.002267, 0.015},
            {-0.002655, 0.001932, 0.015},
            {-0.001909, -0.0021, 0.015},
            {0.001448, -0.001055, 0.015},
            {-0.000656, -0.004163, 0.015},
            {-0.001028, 0.001322, 0.015},
            {-0.000823, 0.000064, 0.015},
            {-0.001577, 0.002143, 0.015},
            {-0.003156, 0.002715, 0.015},
            {-0.000144, -0.000627, 0.015},
            {-0.000337, 0.000516, 0.015},
            {0.004408, -0.000911, 0.015},
            {-0.001765, -0.004383, 0.015},
            {0.001956, -0.00401, 0.015},
            {0.00047, 0.001193, 0.015},
            {-0.001969, -0.001025, 0.015},
            {-0.000721, -0.002182, 0.015},
            {0.002709, 0.003587, 0.015},
            {-0.000945, 0.004862, 0.015},
            {0.004687, -0.001026, 0.015},
            {-0.004885, -0.000268, 0.015},
            {0.001136, 0.003461, 0.015},
            {0.003936, 0.000516, 0.015},
            {-0.000152, 0.001287, 0.015},
            {0.000891, -0.001938, 0.015},
            {-0.002197, -0.001791, 0.015},
            {0.003507, -0.000884, 0.015},
            {0.000804, -0.002557, 0.015},
            {0.004456, -0.001058, 0.015},
            {0.003731, -0.003172, 0.015},
            {-0.0001, 0.003518, 0.015},
            {0.003609, 0.002787, 0.015},
            {-0.001987, -0.004217, 0.015},
            {0.002069, 0.002363, 0.015},
            {-0.000427, 0.003715, 0.015},
            {0.000223, 0.003423, 0.015},
            {0.000745, 0.00195, 0.015},
            {-0.000893, -0.004917, 0.015},
            {-0.000708, 0.001333, 0.015},
            {0.000664, -0.001195, 0.015},
            {0.001506, -0.001444, 0.015},
            {0.000078, -0.001752, 0.015},
            {0.002761, -0.003631, 0.015},
            {0.003172, 0.000248, 0.015},
            {-0.002463, 0.000414, 0.015},
            {0.001406, 0.000663, 0.015},
            {0.001928, -0.003953, 0.015},
            {-0.003837, -0.002558, 0.015},
            {0.002148, 0.00192, 0.015}};
    const std::vector<Vector_t<double, 3>> kvP{
            {0, 0, 0.4},
            {0.000154, 0.00001, 0.4},
            {-0.000139, -0.000005, 0.4},
            {-0.000008, 0.0001, 0.4},
            {-0.000092, -0.000006, 0.4},
            {0.000121, 0.000036, 0.4},
            {-0.000045, -0.00015, 0.4},
            {-0.000063, -0.000124, 0.4},
            {0.000052, 0.000138, 0.4},
            {-0.000088, 0.000111, 0.4},
            {-0.000071, 0.000135, 0.4},
            {-0.000083, 0.00011, 0.4},
            {-0.000077, 0.000097, 0.4},
            {0.000003, -0.000009, 0.4},
            {0.000026, -0.000021, 0.4},
            {0.000029, 0.000002, 0.4},
            {0.000049, -0.000057, 0.4},
            {-0.00012, -0.000098, 0.4},
            {-0.000078, 0.000114, 0.4},
            {0.000039, 0.000045, 0.4},
            {-0.000037, -0.000089, 0.4},
            {0.000006, 0.000032, 0.4},
            {0.000067, -0.000102, 0.4},
            {-0.000072, 0.000046, 0.4},
            {-0.000045, -0.000088, 0.4},
            {-0.000055, -0.000119, 0.4},
            {-0.000032, 0.000001, 0.4},
            {0.00006, -0.000137, 0.4},
            {-0.000006, -0.000156, 0.4},
            {-0.000049, -0.000132, 0.4},
            {-0.000009, 0.000092, 0.4},
            {-0.000065, 0.000089, 0.4},
            {-0.000131, -0.000063, 0.4},
            {0.000095, -0.000117, 0.4},
            {0.000155, 0.00001, 0.4},
            {0.000032, 0.000148, 0.4},
            {-0.00002, 0.00007, 0.4},
            {0.000072, 0.000045, 0.4},
            {-0.000054, 0.000026, 0.4},
            {0.000151, 0.000044, 0.4},
            {-0.000005, -0.000115, 0.4},
            {-0.000125, 0.000053, 0.4},
            {-0.000104, 0.000048, 0.4},
            {-0.000087, -0.000077, 0.4},
            {0.000068, 0.00008, 0.4},
            {-0.000033, -0.000132, 0.4},
            {0.000109, -0.000105, 0.4},
            {-0.000002, -0.000009, 0.4},
            {0.000129, 0.000063, 0.4},
            {0.000023, -0.000097, 0.4},
            {-0.000149, -0.000009, 0.4},
            {0.000074, -0.000115, 0.4},
            {0.000118, -0.000025, 0.4},
            {0.000004, 0.000132, 0.4},
            {0.000081, 0.000125, 0.4},
            {-0.00007, 0.00005, 0.4},
            {-0.000113, -0.0001, 0.4},
            {0.000129, -0.00009, 0.4},
            {0.00006, -0.000122, 0.4},
            {0.000077, 0.000044, 0.4},
            {0.000058, -0.000148, 0.4},
            {0.000124, -0.0001, 0.4},
            {0.000037, -0.000059, 0.4},
            {0.000012, -0.000051, 0.4},
            {0.000072, 0.000002, 0.4},
            {-0.000105, 0.000114, 0.4},
            {0.000021, 0.000142, 0.4},
            {-0.000045, -0.000135, 0.4},
            {0.000009, 0.00007, 0.4},
            {-0.000019, 0.00001, 0.4},
            {0.000045, 0.000006, 0.4},
            {-0.000011, -0.000031, 0.4},
            {0.000034, -0.000104, 0.4},
            {-0.000082, 0.000053, 0.4},
            {-0.000152, 0.00003, 0.4},
            {-0.000008, 0.000144, 0.4},
            {0.000117, -0.00006, 0.4},
            {-0.000041, 0.000103, 0.4},
            {-0.000072, -0.000114, 0.4},
            {0.000002, 0.000064, 0.4},
            {0.000028, 0.000016, 0.4},
            {-0.00001, 0.000113, 0.4},
            {-0.000035, -0.000055, 0.4},
            {-0.000044, 0.000038, 0.4},
            {-0.000071, -0.000103, 0.4},
            {0.000031, 0.000101, 0.4},
            {0.000105, 0.000049, 0.4},
            {-0.000112, -0.000093, 0.4},
            {0.000003, 0.000003, 0.4},
            {0.000123, -0.00009, 0.4},
            {0.000122, -0.000094, 0.4},
            {0.000019, 0.000144, 0.4},
            {0.000029, 0.000147, 0.4},
            {0.000011, -0.000065, 0.4},
            {-0.000101, -0.000071, 0.4},
            {-0.000109, 0.000086, 0.4},
            {-0.00012, -0.000093, 0.4},
            {-0.000068, 0.000034, 0.4},
            {-0.000025, 0.000056, 0.4},
            {0.000105, -0.000078, 0.4}};

    TEST_F(TestSolve2d5, KvBeam_Charge) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 0.03}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(10);
        fsCmd_m->setNY(10);
        fsCmd_m->setNZ(3);
        fsCmd_m->setPipeSizeX(0.02);
        fsCmd_m->setPipeSizeY(0.02);
        fsCmd_m->setClosedRing(false);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles(kvR, kvP);
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterCharge);
        solver->doRunSolver<Info>(*info);
        // Check the charge distribution
        expectTotalCharge(info->rhoView_m, bunch_m->dt_m, {0.0, 0.0, 100.0, 0.0, 0.0}, 1e-6);
    }

    TEST_F(TestSolve2d5, KvBeam_ChargeDensity) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 0.03}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(10);
        fsCmd_m->setNY(10);
        fsCmd_m->setNZ(3);
        fsCmd_m->setPipeSizeX(0.02);
        fsCmd_m->setPipeSizeY(0.02);
        fsCmd_m->setClosedRing(false);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles(kvR, kvP);
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterChargeDensity);
        solver->doRunSolver<Info>(*info);
        // Check the charge distribution
        constexpr double vol = 0.002 * 0.002 * 0.01;
        expectTotalChargeDensity(info->rhoView_m, {0.0, 0.0, 100.0 / vol, 0.0, 0.0}, 1e2);
    }

    TEST_F(TestSolve2d5, KvBeam_Potential) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 0.03}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(10);
        fsCmd_m->setNY(10);
        fsCmd_m->setNZ(3);
        fsCmd_m->setPipeSizeX(0.02);
        fsCmd_m->setPipeSizeY(0.02);
        fsCmd_m->setClosedRing(false);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles(kvR, kvP);
        const auto info = solver->createDiagnostic<Info>(Info::Kind::Potential);
        solver->doRunSolver<Info>(*info);
        // Check the charge density
        expectPotential(               //i, j, k, phi
                info->rhoView_m, {{5, 5, 2, 963538234764917},
                                     {5, 6, 2, 947474778602442.75},
                                     {4, 6, 2, 954504662037936.12},
                                     {6, 7, 2, 948706293379785}}, 10.0);
    }

    TEST_F(TestSolve2d5, KvBeam_ElectricField) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 0.03}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(10);
        fsCmd_m->setNY(10);
        fsCmd_m->setNZ(3);
        fsCmd_m->setPipeSizeX(0.02);
        fsCmd_m->setPipeSizeY(0.02);
        fsCmd_m->setClosedRing(false);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch_m->getFieldSolver());
        createParticles(kvR, kvP);
        const auto info = solver->createDiagnostic<Info>(Info::Kind::EField);
        solver->doRunSolver<Info>(*info);
        // Check the charge density
        expectEField(               //i, j, k, Ex,                  Ey,                   Ez
                info->eFieldView_m, {{5, 5, 2, 22636807214983288.0, 3856068046813647.5,   0.0},
                                     {5, 6, 2, -927805928809394.75, 6103969523576944.0,   0.0},
                                     {4, 6, 2, -5967249946832035.0, 12198160561896932.0,  0.0},
                                     {6, 7, 2, 12256750350007952.0, 20338717682993576.0,  0.0}});
    }
}  // namespace
