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
#include <algorithm>
#include <cmath>
#include <csignal>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <random>
#include <string>
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "Lines/Line.h"
#include "PartBunch/PartBunch.h"
#include "SpaceCharge/FFT2D5/FFT2D5Algorithm.h"
#include "SpaceCharge/FFT2D5/FFT2D5Poisson.h"
#include "SpaceCharge/SpaceChargeConfigBuilder.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

namespace {
    using namespace opalx::spacecharge;

    constexpr bool VerboseTest = false;

    double referenceCellPotential(double x, double y, double hx, double hy);

    // Direct real-space convolution is independent of the production FFT and
    // corner primitive. The deposited source is tested separately by CIC tests.
    template <typename HostView>
    double referenceGridPotential(
            const HostView& source, size_t i, size_t j, size_t k, double hx, double hy) {
        double result = 0.0;
        for (size_t sx = 1; sx + 1 < source.extent(0); ++sx) {
            for (size_t sy = 1; sy + 1 < source.extent(1); ++sy) {
                if (source(sx, sy, k) != 0.0) {
                    result += source(sx, sy, k) * hx * hy / Physics::epsilon_0
                              * referenceCellPotential(
                                      (static_cast<double>(i) - sx) * hx,
                                      (static_cast<double>(j) - sy) * hy, hx, hy);
                }
            }
        }
        return result;
    }

    template <typename HostView>
    std::array<double, 2> referenceGridField(
            const HostView& source, size_t i, size_t j, size_t k, double hx, double hy) {
        return {(referenceGridPotential(source, i - 1, j, k, hx, hy)
                 - referenceGridPotential(source, i + 1, j, k, hx, hy))
                        / (2 * hx),
                (referenceGridPotential(source, i, j - 1, k, hx, hy)
                 - referenceGridPotential(source, i, j + 1, k, hx, hy))
                        / (2 * hy)};
    }

    double referenceTwoParticleForce() {
        // In the 0.5-m slices the two unit-charge particles at x=+/-2 m
        // each deposit line charge 1 C/m. Both transverse CIC shapes have
        // four equal weights. Directly integrate each uniform source cell,
        // central-difference its potential, then gather those four fields.
        static const double force = [] {
            double result = 0.0;
            for (const double tx : {1.75, 2.25}) {
                for (const double ty : {-0.25, 0.25}) {
                    for (const double center : {-2.0, 2.0}) {
                        for (const double sx : {center - 0.25, center + 0.25}) {
                            for (const double sy : {-0.25, 0.25}) {
                                result +=
                                        (referenceCellPotential(tx - sx - 0.5, ty - sy, 0.5, 0.5)
                                         - referenceCellPotential(tx - sx + 0.5, ty - sy, 0.5, 0.5))
                                        / (16.0 * Physics::epsilon_0);
                            }
                        }
                    }
                }
            }
            return result;
        }();
        return force;
    }

    class TestableFieldSolverCmd final : public FieldSolverCmd {
    public:
        void setType(const std::string& t) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::TYPE], t);
        }

        void setBinsName(const std::string& name) {
            Attributes::setString(this->itsAttr[FIELDSOLVER::BINS], name);
        }

        void setParallelDimensions(
                const bool parallelX, const bool parallelY, const bool parallelZ) {
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTX], parallelX);
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTY], parallelY);
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTZ], parallelZ);
        }
    };

    using PartBunch_t         = PartBunch<double, 3>;
    using Solve2d5_t          = FFT2D5Algorithm;
    using ParticleContainer_t = PartBunch_t::ParticleContainer_t;
    using Point_t             = ippl::Vector<double, 3>;

    class Info : public Solve2d5_t::NullDiagnostic {
    public:
        explicit Info(const Kind kind) : kind_m(kind) {}

        void initialise(
                std::span<ParticleContainer_t* const> particles, const Field_t<3>& rho,
                const FFT2D5Algorithm::LineDensityView_t& lineDensity,
                const FFT2D5Algorithm::LineDensityView_t& lineDensityGradient,
                const VField_t<T, 3>& eField) {
            if (!particles.empty()) {
                auto* pc  = particles.front();
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
                const size_t n, const FFT2D5Algorithm::Vector3D_t& r,
                const FFT2D5Algorithm::Vector3D_t& p, const bool invalid) const {
            if (kind_m == Kind::FrenetSerretScatter) {
                r_m[n]       = r;
                p_m[n]       = p;
                invalid_m[n] = invalid;
            }
        }

        KOKKOS_FUNCTION void boostToBeam(
                const size_t n, const FFT2D5Algorithm::Vector3D_t& r,
                const FFT2D5Algorithm::Vector3D_t& p, const bool invalid) const {
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
                const size_t n, const FFT2D5Algorithm::Vector3D_t& r,
                const FFT2D5Algorithm::Vector3D_t& p, const bool invalid) const {
            if (kind_m == Kind::FrenetSerretGather) {
                r_m[n]       = r;
                p_m[n]       = p;
                invalid_m[n] = invalid;
            }
        }

        KOKKOS_FUNCTION void gatherEField(
                const size_t n, const FFT2D5Algorithm::Vector3D_t& e,
                const FFT2D5Algorithm::Vector3D_t& b, const bool invalid) const {
            if (kind_m == Kind::GatherEField) {
                e_m[n]       = e;
                b_m[n]       = b;
                invalid_m[n] = invalid;
            }
        }

        KOKKOS_FUNCTION void deboostFromBeam(
                const size_t n, const FFT2D5Algorithm::Vector3D_t& e,
                const FFT2D5Algorithm::Vector3D_t& b, const bool invalid) const {
            if (kind_m == Kind::Deboosted) {
                e_m[n]       = e;
                b_m[n]       = b;
                invalid_m[n] = invalid;
            }
        }

        KOKKOS_FUNCTION void longitudinalField(
                const size_t n, const FFT2D5Algorithm::Vector3D_t& e,
                const FFT2D5Algorithm::Vector3D_t& b, const bool invalid) const {
            if (kind_m == Kind::LongitudinalField) {
                e_m[n]       = e;
                b_m[n]       = b;
                invalid_m[n] = invalid;
            }
        }

        KOKKOS_FUNCTION void labFrameFields(
                const size_t n, const FFT2D5Algorithm::Vector3D_t& e,
                const FFT2D5Algorithm::Vector3D_t& b, const bool invalid) const {
            if (kind_m == Kind::LabFrameFields) {
                e_m[n]       = e;
                b_m[n]       = b;
                invalid_m[n] = invalid;
            }
        }

        [[nodiscard]] std::tuple<std::vector<Vector_t<double, 3>>, std::vector<Vector_t<double, 3>>>
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

        [[nodiscard]] std::tuple<std::vector<Vector_t<double, 3>>, std::vector<Vector_t<double, 3>>>
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
            fsCmd_m->execute();

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
                    /*integration_method=*/"LF2",
                    makeCartesianDomainConfig(buildSpaceChargeConfig(*fsCmd_m, {})));
            pc_m = bunch_m->getParticleContainer();
        }

        void TearDown() override {
            // Ensure device allocations can be freed between tests.
            solver_m.reset();
            pc_m.reset();
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

        SpaceChargeSolveContext context() const {
            SpaceChargeStepState step;
            step.timeStep = bunch_m->getdT();
            step.mpiSize  = ippl::Comm->size();
            return {activity_m, step};
        }

        void rebuildBunch() {
            solver_m.reset();
            pc_m.reset();
            if (fsCmd_m->getType() == "FFT2D5") {
                fsCmd_m->setParallelDimensions(false, false, false);
            }
            fsCmd_m->execute();
            Beam* testBeam = Beam::find("UNNAMED_BEAM");
            bunch_m        = std::make_shared<PartBunch_t>(
                    /*qi=*/std::vector{1.0},
                    /*mi=*/std::vector{1.0},
                    /*beams=*/std::vector{testBeam},
                    /*totalParticlesPerBeam=*/std::vector{kDefaultNParticles},
                    /*lbt=*/1.0,
                    /*integration_method=*/"LF2",
                    makeCartesianDomainConfig(buildSpaceChargeConfig(*fsCmd_m, {})));
            pc_m        = bunch_m->getParticleContainer();
            auto config = std::get<FFT2D5Config>(buildSpaceChargeConfig(*fsCmd_m, {}));
            const std::array particles{pc_m.get()};
            solver_m = std::make_unique<FFT2D5Algorithm>(
                    config, particles, bunch_m->getBunchStateHandler());
            solver_m->ensureInitialized();
        }

        static void makeReferencePathFile(
                const std::filesystem::path& fileName, const std::vector<Point_t>& points,
                const bool noTimeColumn = false) {
            // Emulates the file created by ReferenceTracer.  Use noTimeColumn = true to create
            // an invalid file.
            // Make the directory path exist
            if (const auto parentPath = fileName.parent_path(); !parentPath.empty()) {
                std::filesystem::create_directories(parentPath);
            }
            // Create the file
            std::ofstream file(fileName);
            if (!file.is_open()) {
                throw std::runtime_error("Failed to open file: " + fileName.string());
            }
            // The header
            file << "#    1 - s    2 - Rx    3 - Ry    4 - Rz    5 - Px    6 - Py    7 - Pz    "
                    "8 - Efx    9 - Efy    10 - Efz    11 - Bfx    12 - Bfy    13 - Bfz    "
                    "14 - Ekin   15 - t"
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
                const size_t index, const std::vector<Vector_t<double, 3>>& rs,
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
                const size_t index, const std::vector<Vector_t<double, 3>>& es,
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

        std::array<std::uint8_t, 1> activity_m{1};
        std::unique_ptr<FFT2D5Algorithm> solver_m;
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

    TEST_F(TestSolve2d5, Configuration_TransverseParallelLayoutRejected) {
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setParallelDimensions(true, false, false);

        EXPECT_NO_THROW(fsCmd_m->execute());
        EXPECT_THROW((void)buildSpaceChargeConfig(*fsCmd_m, {}), OpalException);
    }

    TEST_F(TestSolve2d5, Configuration_BinningRejected) {
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setBinsName("UNUSED_BINNING");
        fsCmd_m->setParallelDimensions(false, false, false);

        EXPECT_NO_THROW(fsCmd_m->execute());
        EXPECT_THROW((void)buildSpaceChargeConfig(*fsCmd_m, {}), OpalException);
    }

    TEST_F(TestSolve2d5, Configuration_MultipleRanksRequireZDecomposition) {
        if (ippl::Comm->size() == 1) {
            GTEST_SKIP() << "This validation requires multiple MPI ranks.";
        }
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setParallelDimensions(false, false, false);
        EXPECT_NO_THROW(fsCmd_m->execute());
        EXPECT_THROW((void)buildSpaceChargeConfig(*fsCmd_m, {}), OpalException);
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
        auto* solver     = solver_m.get();
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
        auto* solver     = solver_m.get();
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
        const auto* solver = solver_m.get();
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
        fsCmd_m->setType("FFT2D5");
        ASSERT_NO_THROW(rebuildBunch());
        auto* solver = solver_m.get();
        EXPECT_NO_THROW(solver->scatterToGrid(context()));
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_ShortReferencePath) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setClosedRing(false);
        EXPECT_THROW(rebuildBunch(), OpalException);
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_TrivialReferencePath) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}});
        fsCmd_m->setType("FFT2D5");
        rebuildBunch();
        auto* solver = solver_m.get();
        createParticles({{0, 0, 0}, {0, 0, 0.5}, {0, 0, 1}}, {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::FrenetSerretScatter);
        solver->scatterToGrid<Info>(context(), *info);
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 2}}, {{0, 0, 1}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::FrenetSerretScatter);
        solver->scatterToGrid<Info>(context(), *info);
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
        auto* solver = solver_m.get();
        createParticles({{2, 0, 2}}, {{0, 0, 1}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::FrenetSerretScatter);
        solver->scatterToGrid<Info>(context(), *info);
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 2}}, {{0, 0, 1}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::FrenetSerretScatter);
        solver->scatterToGrid<Info>(context(), *info);
        auto [r, p] = info->getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {-0.7071068, 0, 1.7071068}, {-0.7071068, 0, 0.7071068});
    }

    TEST_F(TestSolve2d5, BoostToBeamFrame_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd_m->setType("FFT2D5");
        rebuildBunch();
        auto* solver = solver_m.get();
        createParticles({{0, 0.5, 2}}, {{0.001, 0.002, 0.577}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::BoostToBeam);
        solver->scatterToGrid<Info>(context(), *info);
        auto [r, p] = info->getParticles();
        ASSERT_EQ(r.size(), 1);
        expectParticle(0, r, p, {0, 0.5, 2.0}, {0.001, 0.002, 0.0});
    }

    TEST_F(TestSolve2d5, BoostToBeamFrame_SimpleInvalid) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd_m->setType("FFT2D5");
        rebuildBunch();
        auto* solver = solver_m.get();
        createParticles(
                {{0, 0.5, 2}, {0, 0, 2}}, {{0.001, 0.002, 0.577}, {0.001, 0.002, 0.577}},
                {false, true});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::BoostToBeam);
        solver->scatterToGrid<Info>(context(), *info);
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 3}}, {{0, 0, 0}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterCharge);
        solver->scatterToGrid<Info>(context(), *info);
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 0}}, {{0, 0, 0}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterCharge);
        solver->scatterToGrid<Info>(context(), *info);
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 6}}, {{0, 0, 0}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterCharge);
        solver->scatterToGrid<Info>(context(), *info);
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 3}}, {{0, 0, 0}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterChargeDensity);
        solver->scatterToGrid<Info>(context(), *info);
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 0}}, {{0, 0, 0}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterChargeDensity);
        solver->scatterToGrid<Info>(context(), *info);
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 6}}, {{0, 0, 0}});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterChargeDensity);
        solver->scatterToGrid<Info>(context(), *info);
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 3}}, {{0, 0, 0}});
        solver->scatterToGrid(context());
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 3}}, {{0, 0, 0}});
        solver->scatterToGrid(context());
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 3}}, {{0, 0, 0}});
        solver->scatterToGrid(context());
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 5}}, {{0, 0, 0}});
        solver->scatterToGrid(context());
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 5.5}}, {{0, 0, 0}});
        solver->scatterToGrid(context());
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 3}}, {{0, 0, 0}});
        solver->scatterToGrid(context());
#if !OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        const auto referenceSource = Kokkos::create_mirror_view_and_copy(
                Kokkos::HostSpace(), solver->getRho()->getView());
#endif
        const auto info = solver->createDiagnostic<Info>(Info::Kind::EField);
        solver->solvePoissons<Info>(*info);
#if OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        expectEField(
                info->eFieldView_m, {{1, 1, 6, -3250090326, -3250090326, 0},
                                     {2, 1, 6, -3218709596, -3896406064, 0},
                                     {4, 4, 6, -7214724490, -7214724490, 0},
                                     {4, 8, 6, -10565567441, 6311498630, 0}});
#else
        const auto field =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), info->eFieldView_m);
        for (const auto& point :
             std::array<std::array<size_t, 2>, 4>{{{1, 1}, {2, 1}, {4, 4}, {4, 8}}}) {
            const auto expected =
                    referenceGridField(referenceSource, point[0], point[1], 6, 0.5, 0.5);
            EXPECT_NEAR(field(point[0], point[1], 6)[0], expected[0], 1e3);
            EXPECT_NEAR(field(point[0], point[1], 6)[1], expected[1], 1e3);
        }
#endif
    }

    TEST_F(TestSolve2d5, SliceWorkspaceKeepsNeutralSourcesAndClearsVacatedSlices) {
        for (int batch : {0, 8}) {
            makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
            fsCmd_m->setType("FFT2D5");
            fsCmd_m->setNX(12);
            fsCmd_m->setNY(12);
            fsCmd_m->setNZ(12);
            fsCmd_m->setPipeSizeX(6);
            fsCmd_m->setPipeSizeY(6);
            Options::fft2d5BatchSize = batch;
            rebuildBunch();
            Options::fft2d5BatchSize = 0;
            const auto rho           = solver_m->getRho()->getView();
            const auto field         = solver_m->getEField()->getView();
            auto source              = Kokkos::create_mirror_view(rho);
            Kokkos::deep_copy(source, 0.0);
            // This slice has zero signed total charge but a nonzero dipole field.
            source(3, 3, 4) = 1.0;
            source(8, 7, 4) = -1.0;
            Kokkos::deep_copy(rho, source);
            Kokkos::deep_copy(field, Vector_t<double, 3>(7.0));
            solver_m->solvePoissons();
            const auto first = Kokkos::create_mirror(field);
            Kokkos::deep_copy(first, field);
            double maximum = 0.0;
            for (size_t x = 0; x < first.extent(0); ++x) {
                for (size_t y = 0; y < first.extent(1); ++y) {
                    for (unsigned d = 0; d < 3; ++d) {
                        maximum = std::max(maximum, std::abs(first(x, y, 4)[d]));
                        EXPECT_DOUBLE_EQ(first(x, y, 9)[d], 0.0);
                    }
                }
            }
            EXPECT_GT(maximum, 0.0);

            // Move the same source to another slice; reusing the solver must leave no old field.
            Kokkos::deep_copy(source, 0.0);
            source(3, 3, 9) = 1.0;
            source(8, 7, 9) = -1.0;
            Kokkos::deep_copy(rho, source);
            solver_m->solvePoissons();
            const auto moved = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), field);
            for (size_t x = 0; x < moved.extent(0); ++x) {
                for (size_t y = 0; y < moved.extent(1); ++y) {
                    for (size_t z = 0; z < moved.extent(2); ++z) {
                        for (unsigned d = 0; d < 3; ++d) {
                            if (z == 9) {
                                EXPECT_NEAR(
                                        moved(x, y, z)[d], first(x, y, 4)[d], maximum * 1.0e-13);
                            } else {
                                EXPECT_DOUBLE_EQ(moved(x, y, z)[d], 0.0);
                            }
                        }
                    }
                }
            }
            Kokkos::deep_copy(rho, 0.0);
            solver_m->solvePoissons();
            const auto empty = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), field);
            for (size_t x = 0; x < empty.extent(0); ++x) {
                for (size_t y = 0; y < empty.extent(1); ++y) {
                    for (size_t z = 0; z < empty.extent(2); ++z) {
                        for (unsigned d = 0; d < 3; ++d) {
                            EXPECT_DOUBLE_EQ(empty(x, y, z)[d], 0.0);
                        }
                    }
                }
            }
        }
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
        auto* solver = solver_m.get();
        createParticles({{0, 0, 2}}, {{0, 0, 1}});
        solver->scatterToGrid(context());
        solver->solvePoissons();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::FrenetSerretGather);
        solver->gatherFromGrid<Info>(context(), *info);
        auto [r, p] = info->getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {-0.7071068, 0, 1.7071068}, {-0.7071068, 0, 0.7071068});
    }

    TEST_F(TestSolve2d5, GatherEField_OneParticle) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        rebuildBunch();
        auto* solver = solver_m.get();
        createParticles({{0.1, 0, 3}}, {{0, 0, 0}});
        solver->scatterToGrid(context());
        solver->solvePoissons();
        const auto gatherInfo = solver->createDiagnostic<Info>(Info::Kind::GatherEField);
        solver->gatherFromGrid<Info>(context(), *gatherInfo);
        auto [e, b] = gatherInfo->getParticleFields();
        ASSERT_EQ(e.size(), 1);
        ASSERT_EQ(b.size(), 1);
        expectParticleFields(0, e, b, {0, 0, 0}, {0, 0, 0}, 1e3);
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
        auto* solver = solver_m.get();
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 0}, {0, 0, 0}});
        solver->scatterToGrid(context());
        solver->solvePoissons();
        const auto gatherInfo = solver->createDiagnostic<Info>(Info::Kind::GatherEField);
        solver->gatherFromGrid<Info>(context(), *gatherInfo);
        auto [e, b] = gatherInfo->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
#if OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        expectParticleFields(0, e, b, {4493908875, 0, 0}, {0, 0, 0}, 1e3);
        expectParticleFields(1, e, b, {-4493908875, 0, 0}, {0, 0, 0}, 1e3);
#else
        expectParticleFields(0, e, b, {referenceTwoParticleForce(), 0, 0}, {0, 0, 0}, 1e3);
        expectParticleFields(1, e, b, {-referenceTwoParticleForce(), 0, 0}, {0, 0, 0}, 1e3);
#endif
    }

    // Fixed laboratory charge and geometry must give the same transverse E at
    // every beam speed. The self magnetic force must cancel beta^2 of E.
    TEST_F(TestSolve2d5, TransverseLabDensityAndMagneticCancellation) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);
        fsCmd_m->setPipeMode("NONE");
        std::vector<Point_t> staticE;
        for (const double ps : {0.0, 0.4, 1.0, -0.4}) {
            SCOPED_TRACE(ps);
            rebuildBunch();
            createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, ps}, {0, 0, ps}});
            solver_m->doRunSolver(context());
            const auto e =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc_m->E.getView());
            const auto b =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc_m->B.getView());
            const double gamma = std::sqrt(1.0 + ps * ps);
            const double beta  = ps / gamma;
            if (ps == 0.0) staticE = {e(0), e(1)};
            for (std::size_t i = 0; i < 2; ++i) {
                for (unsigned d = 0; d < 2; ++d) {
                    EXPECT_NEAR(e(i)[d], staticE[i][d], 1e-10 * std::abs(staticE[i][0]));
                }
                EXPECT_NEAR(b(i)[0], -beta * e(i)[1] / Physics::c, 1e-10);
                EXPECT_NEAR(b(i)[1], beta * e(i)[0] / Physics::c, 1e-10);
                const double forceX = e(i)[0] - beta * Physics::c * b(i)[1];
                EXPECT_NEAR(
                        forceX, staticE[i][0] / (gamma * gamma), 1e-10 * std::abs(staticE[i][0]));
            }
        }
    }

    TEST_F(TestSolve2d5, UniformCylinderMatchesGaussLaw) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(128);
        fsCmd_m->setNY(128);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(2);
        fsCmd_m->setPipeSizeY(2);
        fsCmd_m->setPipeMode("NONE");
        rebuildBunch();
        constexpr double radius = 0.5;
        constexpr int rings = 32, angles = 128, slices = 12;
        std::vector<Point_t> positions, momenta;
        for (int k = 0; k < slices; ++k) {
            for (int j = 0; j < rings; ++j) {
                const double r = radius * std::sqrt((j + 0.5) / rings);
                for (int i = 0; i < angles; ++i) {
                    const double phi = 2 * Physics::pi * (i + 0.5) / angles;
                    positions.emplace_back(r * std::cos(phi), r * std::sin(phi), (k + 0.5) * 0.5);
                    momenta.emplace_back(0.0, 0.0, 0.4);
                }
            }
        }
        createParticles(positions, momenta);
        solver_m->doRunSolver(context());
        const auto e = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc_m->E.getView());
        const auto b = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc_m->B.getView());
        const double lambda  = positions.size() * pc_m->getChargePerParticle() / 6.0;
        const double slope   = lambda / (2 * Physics::pi * Physics::epsilon_0 * radius * radius);
        const double beta    = 0.4 / std::sqrt(1.16);
        double radiusSquared = 0, electricMoment = 0, forceMoment = 0;
        for (std::size_t i = 0; i < positions.size(); ++i) {
            const auto& r   = positions[i];
            const double r2 = r[0] * r[0] + r[1] * r[1];
            if (r[2] < 2 || r[2] > 4 || r2 > 0.16) continue;
            radiusSquared += r2;
            electricMoment += r[0] * e(i)[0] + r[1] * e(i)[1];
            forceMoment += r[0] * (e(i)[0] - beta * Physics::c * b(i)[1])
                           + r[1] * (e(i)[1] + beta * Physics::c * b(i)[0]);
        }
        ASSERT_GT(radiusSquared, 0);
        // 2% allows CIC/finite-ring discretization at 32 cells per edge radius;
        // it is well below the old 7.7% E and 42% Lorentz-force excesses.
        EXPECT_NEAR(electricMoment / (slope * radiusSquared), 1.0, 0.02);
        EXPECT_NEAR(forceMoment / (slope * radiusSquared / 1.16), 1.0, 0.02);
    }

    // Analytic grid density removes source-particle sampling noise. This checks
    // the existing zeroth-order slice model on straight and circular paths;
    // Qiang's O(sigma/R) curvature increment is recorded separately, not assumed
    // to be present in the independent flat transverse solves.
    TEST_F(TestSolve2d5, SmoothGaussianSlicesConvergeWithoutParticleNoise) {
        constexpr double sigma = 1e-3, lambda = 1e-9, ps = 0.4;
        constexpr double box = 12 * sigma, centerZ = 0.2, halfPath = 0.1;
        constexpr int slices   = 4;
        const double beta      = ps / std::sqrt(1 + ps * ps);
        const double fieldUnit = lambda / (4 * Physics::pi * Physics::epsilon_0 * sigma);
        std::ofstream csv;
        if (const char* filename = std::getenv("OPALX_FFT2D5_FORCE_CSV")) {
            csv.open(filename);
            ASSERT_TRUE(csv.is_open());
            csv << std::setprecision(17)
                << "mesh,radius_m,x_over_sigma,Ex_V_per_m,By_T,force_per_charge_V_per_m,"
                   "straight_Ex_V_per_m,qiang_delta_Ex_V_per_m\n";
        }
        // Differentiate the regularized increment in Qiang Eq. (98). A fixed
        // Simpson quadrature suffices for this smooth one-dimensional integrand.
        const auto curvatureCoefficient = [](const double x) {
            constexpr int intervals = 1024;
            double sum              = 0;
            for (int j = 1; j <= intervals; ++j) {
                const double w = static_cast<double>(j) / intervals;
                const double a = -0.5 * x * x * w * w;
                const double value =
                        (1 - w * w) / w * std::expm1(a) - (1 - w * w) * x * x * w * std::exp(a);
                sum += (j == intervals ? 1 : (j % 2 ? 4 : 2)) * value;
            }
            return sum / (3 * intervals);
        };
        for (const double radius : {0.0, 1.0, 10.0, 100.0}) {
            double previousError = 1;
            for (const int cells : {64, 128, 256, 512}) {
                SCOPED_TRACE(testing::Message() << "mesh=" << cells << ", R=" << radius);
                // The path is in laboratory coordinates; keep the midpoint
                // tangent along +z so this benchmark isolates the local slice.
                {
                    std::filesystem::create_directories("data");
                    std::ofstream path("data/unit_test_DesignPath.dat");
                    path << std::setprecision(17) << "s rx ry rz px py pz ex ey ez bx by bz ke t\n";
                    const int segments = radius == 0 ? 2 : 2000;
                    for (int j = 0; j <= segments; ++j) {
                        const double s     = -halfPath + 2 * halfPath * j / segments;
                        const double theta = radius == 0 ? 0 : s / radius;
                        const double x     = radius == 0 ? 0 : radius * (std::cos(theta) - 1);
                        const double z     = centerZ + (radius == 0 ? s : radius * std::sin(theta));
                        path << s + halfPath << ' ' << x << " 0 " << z
                             << " 0 0 0 0 0 0 0 0 0 0 0\n";
                    }
                }
                fsCmd_m->setType("FFT2D5");
                fsCmd_m->setNX(cells);
                fsCmd_m->setNY(cells);
                fsCmd_m->setNZ(slices);
                fsCmd_m->setPipeSizeX(box);
                fsCmd_m->setPipeSizeY(box);
                fsCmd_m->setPipeMode("NONE");
                fsCmd_m->setScatterLongitudinally(true);
                rebuildBunch();
                std::vector<Point_t> probes, momenta;
                for (int j = -12; j <= 12; ++j) {
                    probes.emplace_back(j * 0.25 * sigma, 0.0, centerZ);
                    momenta.emplace_back(0.0, 0.0, ps);
                }
                createParticles(probes, momenta);
                auto* rho          = solver_m->getRho();
                const auto density = rho->getHostMirror();
                Kokkos::deep_copy(density, 0.0);
                const int ng    = rho->getNghost();
                const double dx = box / cells;
                for (int i = 0; i < cells; ++i) {
                    const double x = -box / 2 + (i + 0.5) * dx;
                    for (int j = 0; j < cells; ++j) {
                        const double y     = -box / 2 + (j + 0.5) * dx;
                        const double value = lambda / (2 * Physics::pi * sigma * sigma)
                                             * std::exp(-(x * x + y * y) / (2 * sigma * sigma));
                        for (int k = 0; k < slices; ++k)
                            density(ng + i, ng + j, ng + k) = value;
                    }
                }
                Kokkos::deep_copy(rho->getView(), density);
                pc_m->markMomentsDirty();
                pc_m->updateMoments();
                solver_m->solvePoissons();
                solver_m->calculateLineDensity();
                solver_m->gatherFromGrid(context());
                const auto e =
                        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc_m->E.getView());
                const auto b =
                        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc_m->B.getView());
                double error2 = 0, reference2 = 0;
                for (std::size_t i = 0; i < probes.size(); ++i) {
                    const double x = probes[i][0] / sigma;
                    const double reference =
                            x == 0 ? 0 : 2 * fieldUnit * -std::expm1(-x * x / 2) / x;
                    const double increment =
                            radius == 0 ? 0 : fieldUnit * sigma / radius * curvatureCoefficient(x);
                    const double force = e(i)[0] - beta * Physics::c * b(i)[1];
                    if (std::abs(x) <= 2) {
                        error2 += std::pow(e(i)[0] - reference, 2);
                        reference2 += reference * reference;
                    }
                    // The two transverse kernels and probe locations have the
                    // same normalization; no fitted scale is applied.
                    EXPECT_NEAR(force, e(i)[0] / (1 + ps * ps), 1e-7 * fieldUnit);
                    EXPECT_NEAR(e(i)[2], 0.0, 1e-3 * fieldUnit);
                    if (csv)
                        csv << cells << ',' << radius << ',' << x << ',' << e(i)[0] << ','
                            << b(i)[1] << ',' << force << ',' << reference << ',' << increment
                            << '\n';
                }
                const double error = std::sqrt(error2 / reference2);
                std::cout << "SMOOTH_GAUSSIAN mesh=" << cells << " radius=" << radius
                          << " relative_rms_error=" << error << std::endl;
                // The logarithmic Green function, midpoint source sampling and
                // linear gather give approximately second-order convergence.
                // Measured errors at 64/128/256 cells are 6.92/1.86/0.50%;
                // require a reduction below 0.35 per refinement and <0.2% at
                // 512 cells. These are new-test convergence bounds; no existing
                // snapshot or physics tolerance is relaxed.
                if (cells > 64) EXPECT_LT(error, 0.35 * previousError);
                if (cells == 512) EXPECT_LT(error, 0.002);
                previousError = error;
            }
        }
    }

    TEST_F(TestSolve2d5, OpenReferencePathRejectsEndpointPileup) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setClosedRing(false);
        for (const double z : {-0.1, 6.1}) {
            SCOPED_TRACE(z);
            rebuildBunch();
            createParticles({{0, 0, z}}, {{0, 0, 0.4}});
            EXPECT_THROW(solver_m->scatterToGrid(context()), OpalException);
        }
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
        auto* solver = solver_m.get();
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 0}, {0, 0, 0}});
        solver->scatterToGrid(context());
        solver->solvePoissons();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::Deboosted);
        solver->gatherFromGrid<Info>(context(), *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
#if OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        expectParticleFields(0, e, b, {4493908875, 0, 0}, {0, 0, 0}, 1e3);
        expectParticleFields(1, e, b, {-4493908875, 0, 0}, {0, 0, 0}, 1e3);
#else
        expectParticleFields(0, e, b, {referenceTwoParticleForce(), 0, 0}, {0, 0, 0}, 1e3);
        expectParticleFields(1, e, b, {-referenceTwoParticleForce(), 0, 0}, {0, 0, 0}, 1e3);
#endif
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
        auto* solver = solver_m.get();
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(context());
        solver->solvePoissons();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::Deboosted);
        solver->gatherFromGrid<Info>(context(), *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
#if OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        expectParticleFields(
                0, e, b, {4493908875, 0, 0}, {0, (4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0},
                1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-4493908875, 0, 0},
                {0, -(4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
#else
        expectParticleFields(
                0, e, b, {referenceTwoParticleForce(), 0, 0},
                {0, (referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-referenceTwoParticleForce(), 0, 0},
                {0, -(referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
#endif
    }

    TEST_F(TestSolve2d5, ParabolicLineDensityHasCenteredLongitudinalForce) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(4);
        fsCmd_m->setNY(4);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(2);
        fsCmd_m->setPipeSizeY(2);
        fsCmd_m->setPipeMode("OPEN");
        constexpr double lambda0 = 1e-9, halfLength = 3.0, dz = 0.5;
        for (const bool scatter : {false, true}) {
            fsCmd_m->setScatterLongitudinally(scatter);
            for (const double ps : {0.0, 0.4, 1.0}) {
                SCOPED_TRACE(testing::Message() << "CIC=" << scatter << ", ps=" << ps);
                rebuildBunch();
                const std::vector<Point_t> probes =
                        scatter ? std::vector<Point_t>{{0, 0, 2.1}, {0, 0, 3}, {0, 0, 3.9}}
                                : std::vector<Point_t>{
                                          {0, 0, 2.25}, {0, 0, 2.75}, {0, 0, 3.25}, {0, 0, 3.75}};
                createParticles(probes, std::vector<Point_t>(probes.size(), Point_t(0, 0, ps)));
                auto* rho       = solver_m->getRho();
                const auto host = rho->getHostMirror();
                Kokkos::deep_copy(host, 0.0);
                const int ng = rho->getNghost();
                for (int k = 0; k < 12; ++k) {
                    const double zeta   = (k + 0.5) * dz - halfLength;
                    const double lambda = lambda0 * (1 - zeta * zeta / (halfLength * halfLength));
                    // One transverse cell, dx=dy=0.5 m; integral rho dx dy=lambda.
                    host(ng, ng, ng + k) = lambda / 0.25;
                }
                Kokkos::deep_copy(rho->getView(), host);
                Kokkos::deep_copy(solver_m->getEField()->getView(), Point_t(0, 0, 0));
                pc_m->markMomentsDirty();
                pc_m->updateMoments();
                solver_m->calculateLineDensity();
                solver_m->gatherFromGrid(context());
                const auto e =
                        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc_m->E.getView());
                const double coefficient = 6.36 * lambda0
                                           / (2 * Physics::pi * Physics::epsilon_0 * (1 + ps * ps)
                                              * halfLength * halfLength);
                double total = 0;
                for (std::size_t i = 0; i < probes.size(); ++i) {
                    EXPECT_NEAR(e(i)[2], coefficient * (probes[i][2] - halfLength), 1e-9);
                    total += e(i)[2];
                }
                EXPECT_NEAR(total, 0.0, 1e-9);
            }
        }
    }

    TEST_F(TestSolve2d5, ClosedLongitudinalGatherWrapsCellCenters) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNZ(12);
        fsCmd_m->setClosedRing(true);
        fsCmd_m->setScatterLongitudinally(true);
        rebuildBunch();
        createParticles({{0, 0, 0}, {0, 0, 6}}, {{0, 0, 0}, {0, 0, 0}});
        const auto gradient = solver_m->getLineDensityGradient();
        const auto host     = Kokkos::create_mirror_view(gradient);
        for (std::size_t k = 0; k < host.extent(0); ++k)
            host(k) = (k + 1) * 1e-9;
        Kokkos::deep_copy(gradient, host);
        Kokkos::deep_copy(solver_m->getEField()->getView(), Point_t(0, 0, 0));
        solver_m->gatherFromGrid(context());
        const auto e = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc_m->E.getView());
        const double expected = -6.36 / (4 * Physics::pi * Physics::epsilon_0) * 6.5e-9;
        EXPECT_NEAR(e(0)[2], expected, 1e-9);
        EXPECT_NEAR(e(1)[2], expected, 1e-9);
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
        auto* solver = solver_m.get();
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(context());
        solver->solvePoissons();
        const auto lineInfo = solver->createDiagnostic<Info>(Info::Kind::LineDensityGradient);
        solver->calculateLineDensity<Info>(*lineInfo);
        expectLineDensity(lineInfo->lineDensityGradientView_m, {});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LongitudinalField);
        solver->gatherFromGrid<Info>(context(), *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
#if OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        expectParticleFields(
                0, e, b, {4493908875, 0, 0.0},
                {0, (4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-4493908875, 0, 0.0},
                {0, -(4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
#else
        expectParticleFields(
                0, e, b, {referenceTwoParticleForce(), 0, 0.0},
                {0, (referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-referenceTwoParticleForce(), 0, 0.0},
                {0, -(referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
#endif
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
        auto* solver = solver_m.get();
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(context());
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(context(), *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
#if OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        expectParticleFields(
                0, e, b, {4493908875, 0, 0.0},
                {0, (4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-4493908875, 0, 0.0},
                {0, -(4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
#else
        expectParticleFields(
                0, e, b, {referenceTwoParticleForce(), 0, 0.0},
                {0, (referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-referenceTwoParticleForce(), 0, 0.0},
                {0, -(referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
#endif
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
        auto* solver = solver_m.get();
        createParticles(
                {{2, 0, 3}, {-2, 0, 3}, {-2, 0, 3}}, {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}},
                {false, false, true});
        solver->scatterToGrid(context());
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(context(), *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 3);
        ASSERT_EQ(b.size(), 3);
#if OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        expectParticleFields(
                0, e, b, {4493908875, 0, 0.0},
                {0, (4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-4493908875, 0, 0.0},
                {0, -(4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
        expectParticleFields(2, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
#else
        expectParticleFields(
                0, e, b, {referenceTwoParticleForce(), 0, 0.0},
                {0, (referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-referenceTwoParticleForce(), 0, 0.0},
                {0, -(referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
        expectParticleFields(2, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
#endif
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
        auto* solver = solver_m.get();
        createParticles({{4, 0, 3}, {-4, 0, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(context());
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(context(), *info);
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
        auto* solver = solver_m.get();
        createParticles({{0, 4, 3}, {0, -4, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(context());
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(context(), *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
        expectParticleFields(0, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
        expectParticleFields(1, e, b, {0, 0, 0}, {0, 0, 0}, 1e3, 1e-4);
    }

    TEST_F(TestSolve2d5, LabFrameFields_SimpleOutOfBoundsZ) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        rebuildBunch();
        createParticles({{0, 0, -1}, {0, 0, 7}}, {{0, 0, 1}, {0, 0, 1}});
        EXPECT_THROW(solver_m->scatterToGrid(context()), OpalException);
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
        auto* solver = solver_m.get();
        createParticles({{2, 0, 5}, {-2, 0, 5}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(context());
        const auto efieldInfo = solver->createDiagnostic<Info>(Info::Kind::EField);
        solver->solvePoissons<Info>(*efieldInfo);
        expectEField(efieldInfo->eFieldView_m, {});
        const auto lineInfo = solver->createDiagnostic<Info>(Info::Kind::LineDensityGradient);
        solver->calculateLineDensity<Info>(*lineInfo);
        expectLineDensity(lineInfo->lineDensityGradientView_m, {});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(context(), *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 3);
        ASSERT_EQ(b.size(), 3);
#if OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        expectParticleFields(
                0, e, b, {4493908875, 0, -6.36 / (16 * Physics::pi * Physics::epsilon_0)},
                {0, (4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e6, 1e-4);
        expectParticleFields(
                1, e, b, {-4493908875, 0, -6.36 / (16 * Physics::pi * Physics::epsilon_0)},
                {0, -(4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e6, 1e-4);
        expectParticleFields(2, e, b, {0, 0, 28.580515699e9}, {0, 0, 0}, 1e6, 1e-4);
#else
        expectParticleFields(
                0, e, b,
                {referenceTwoParticleForce(), 0, -6.36 / (16 * Physics::pi * Physics::epsilon_0)},
                {0, (referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e6, 1e-4);
        expectParticleFields(
                1, e, b,
                {-referenceTwoParticleForce(), 0, -6.36 / (16 * Physics::pi * Physics::epsilon_0)},
                {0, -(referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e6, 1e-4);
        expectParticleFields(2, e, b, {0, 0, 28.580515699e9}, {0, 0, 0}, 1e6, 1e-4);
#endif
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
        auto* solver = solver_m.get();
        createParticles({{2, 0, 5}, {-2, 0, 5}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(context());
        const auto efieldInfo = solver->createDiagnostic<Info>(Info::Kind::EField);
        solver->solvePoissons<Info>(*efieldInfo);
        expectEField(efieldInfo->eFieldView_m, {});
        const auto lineInfo = solver->createDiagnostic<Info>(Info::Kind::LineDensityGradient);
        solver->calculateLineDensity<Info>(*lineInfo);
        expectLineDensity(lineInfo->lineDensityGradientView_m, {});
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(context(), *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 3);
        ASSERT_EQ(b.size(), 3);
#if OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        expectParticleFields(
                0, e, b, {4493908875, 0, 0.0},
                {0, (4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e6, 1e-4);
        expectParticleFields(
                1, e, b, {-4493908875, 0, 0.0},
                {0, -(4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e6, 1e-4);
        expectParticleFields(2, e, b, {0, 0, 0}, {0, 0, 0}, 1e6, 1e-4);
#else
        expectParticleFields(
                0, e, b, {referenceTwoParticleForce(), 0, 0.0},
                {0, (referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e6, 1e-4);
        expectParticleFields(
                1, e, b, {-referenceTwoParticleForce(), 0, 0.0},
                {0, -(referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e6, 1e-4);
        expectParticleFields(2, e, b, {0, 0, 0}, {0, 0, 0}, 1e6, 1e-4);
#endif
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
        auto* solver = solver_m.get();
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(context());
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(context(), *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
#if OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        expectParticleFields(
                0, e, b, {4493908875, 0, 0.0},
                {0, (4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-4493908875, 0, 0.0},
                {0, -(4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
#else
        expectParticleFields(
                0, e, b, {referenceTwoParticleForce(), 0, 0.0},
                {0, (referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-referenceTwoParticleForce(), 0, 0.0},
                {0, -(referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
#endif
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
        auto* solver = solver_m.get();
        createParticles({{2, 0, 3}, {-2, 0, 3}}, {{0, 0, 1}, {0, 0, 1}});
        solver->scatterToGrid(context());
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(context(), *info);
        auto [e, b] = info->getParticleFields();
        ASSERT_EQ(e.size(), 2);
        ASSERT_EQ(b.size(), 2);
#if OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        expectParticleFields(
                0, e, b, {4493908875, 0, 0.0},
                {0, (4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-4493908875, 0, 0.0},
                {0, -(4493908875.0 / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
#else
        expectParticleFields(
                0, e, b, {referenceTwoParticleForce(), 0, 0.0},
                {0, (referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
        expectParticleFields(
                1, e, b, {-referenceTwoParticleForce(), 0, 0.0},
                {0, -(referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)), 0}, 1e3, 1e-4);
#endif
    }

    TEST_F(TestSolve2d5, LoadReferencePath_Empty) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setClosedRing(false);
        EXPECT_THROW(rebuildBunch(), OpalException);
    }

    TEST_F(TestSolve2d5, LoadReferencePath_EmptyClosedRing) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setClosedRing(true);
        EXPECT_THROW(rebuildBunch(), OpalException);
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
        auto* solver = solver_m.get();
        solver->scatterToGrid(context());
        solver->solvePoissons();
        solver->calculateLineDensity();
        const auto info = solver->createDiagnostic<Info>(Info::Kind::LabFrameFields);
        solver->gatherFromGrid<Info>(context(), *info);
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
        auto* solver = solver_m.get();
        createParticles({{2, 0, 5}, {-2, 0, 5}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}});
        (void)solver->solve(context());
        // Check the fields through the original particle bunch object
        auto& pc         = *bunch_m->getParticleContainers().front();
        const auto eHost = pc.E.getHostMirror();
        const auto bHost = pc.B.getHostMirror();
        Kokkos::deep_copy(eHost, pc.E.getView());
        Kokkos::deep_copy(bHost, pc.B.getView());
#if OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        EXPECT_NEAR(eHost(0).data_m[0], 4493908875, 1e6);
        EXPECT_NEAR(eHost(0).data_m[1], 0, 1e6);
        EXPECT_NEAR(eHost(0).data_m[2], -6.36 / (16 * Physics::pi * Physics::epsilon_0), 1e6);
        EXPECT_NEAR(eHost(1).data_m[0], -4493908875, 1e6);
        EXPECT_NEAR(eHost(1).data_m[1], 0, 1e6);
        EXPECT_NEAR(eHost(1).data_m[2], -6.36 / (16 * Physics::pi * Physics::epsilon_0), 1e6);
        EXPECT_NEAR(eHost(2).data_m[0], 0, 1e6);
        EXPECT_NEAR(eHost(2).data_m[1], 0, 1e6);
        EXPECT_NEAR(eHost(2).data_m[2], 28.580515699e9, 1e6);
        EXPECT_NEAR(bHost(0).data_m[0], 0, 1e-4);
        EXPECT_NEAR(bHost(0).data_m[1], (4493908875.0 / (std::sqrt(2.0) * Physics::c)), 1e-4);
        EXPECT_NEAR(bHost(0).data_m[2], 0, 1e-4);
        EXPECT_NEAR(bHost(1).data_m[0], 0, 1e-4);
        EXPECT_NEAR(bHost(1).data_m[1], -(4493908875.0 / (std::sqrt(2.0) * Physics::c)), 1e-4);
        EXPECT_NEAR(bHost(1).data_m[2], 0, 1e-4);
        EXPECT_NEAR(bHost(2).data_m[0], 0, 1e-4);
        EXPECT_NEAR(bHost(2).data_m[1], 0, 1e-4);
        EXPECT_NEAR(bHost(2).data_m[2], 0, 1e-4);
#else
        EXPECT_NEAR(eHost(0).data_m[0], referenceTwoParticleForce(), 1e6);
        EXPECT_NEAR(eHost(0).data_m[1], 0, 1e6);
        EXPECT_NEAR(eHost(0).data_m[2], -6.36 / (16 * Physics::pi * Physics::epsilon_0), 1e6);
        EXPECT_NEAR(eHost(1).data_m[0], -referenceTwoParticleForce(), 1e6);
        EXPECT_NEAR(eHost(1).data_m[1], 0, 1e6);
        EXPECT_NEAR(eHost(1).data_m[2], -6.36 / (16 * Physics::pi * Physics::epsilon_0), 1e6);
        EXPECT_NEAR(eHost(2).data_m[0], 0, 1e6);
        EXPECT_NEAR(eHost(2).data_m[1], 0, 1e6);
        EXPECT_NEAR(eHost(2).data_m[2], 28.580515699e9, 1e6);
        EXPECT_NEAR(bHost(0).data_m[0], 0, 1e-4);
        EXPECT_NEAR(
                bHost(0).data_m[1], (referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)),
                1e-4);
        EXPECT_NEAR(bHost(0).data_m[2], 0, 1e-4);
        EXPECT_NEAR(bHost(1).data_m[0], 0, 1e-4);
        EXPECT_NEAR(
                bHost(1).data_m[1], -(referenceTwoParticleForce() / (std::sqrt(2.0) * Physics::c)),
                1e-4);
        EXPECT_NEAR(bHost(1).data_m[2], 0, 1e-4);
        EXPECT_NEAR(bHost(2).data_m[0], 0, 1e-4);
        EXPECT_NEAR(bHost(2).data_m[1], 0, 1e-4);
        EXPECT_NEAR(bHost(2).data_m[2], 0, 1e-4);
#endif
    }

    TEST_F(TestSolve2d5, MainApiMatchesDirectStagesForScatterAndBoundaryModes) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 6}});
        fsCmd_m->setType("FFT2D5");
        fsCmd_m->setNX(12);
        fsCmd_m->setNY(12);
        fsCmd_m->setNZ(12);
        fsCmd_m->setPipeSizeX(6);
        fsCmd_m->setPipeSizeY(6);

        for (const bool scatterLongitudinally : {false, true}) {
            for (const bool closedRing : {false, true}) {
                SCOPED_TRACE(
                        testing::Message() << "scatterLongitudinally=" << scatterLongitudinally
                                           << ", closedRing=" << closedRing);
                fsCmd_m->setScatterLongitudinally(scatterLongitudinally);
                fsCmd_m->setClosedRing(closedRing);
                rebuildBunch();
                createParticles(
                        {{2, 0, 5.9}, {-2, 0, 0.1}, {0, 0, 3}}, {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}});

                // Exercise the library entry before the direct stage calls in an isolated run.
                const auto result = solver_m->solve(context());
                EXPECT_EQ(result.backendSolves, scatterLongitudinally ? 4u : 2u);
                const auto apiE = Kokkos::create_mirror(pc_m->E.getView());
                const auto apiB = Kokkos::create_mirror(pc_m->B.getView());
                Kokkos::deep_copy(apiE, pc_m->E.getView());
                Kokkos::deep_copy(apiB, pc_m->B.getView());

                clearSelfFields(*pc_m);
                solver_m->scatterToGrid(context());
                solver_m->solvePoissons();
                solver_m->calculateLineDensity();
                solver_m->gatherFromGrid(context());
                const auto directE = pc_m->E.getHostMirror();
                const auto directB = pc_m->B.getHostMirror();
                Kokkos::deep_copy(directE, pc_m->E.getView());
                Kokkos::deep_copy(directB, pc_m->B.getView());
                for (size_t particle = 0; particle < 3; ++particle) {
                    for (unsigned dimension = 0; dimension < 3; ++dimension) {
                        const auto expectedE = directE(particle)[dimension];
                        const auto expectedB = directB(particle)[dimension];
                        EXPECT_NEAR(
                                apiE(particle)[dimension], expectedE,
                                1e-10 * std::max(1.0, std::abs(expectedE)));
                        EXPECT_NEAR(
                                apiB(particle)[dimension], expectedB,
                                1e-10 * std::max(1.0, std::abs(expectedB)));
                    }
                }
            }
        }
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
        auto* solver = solver_m.get();
        createParticles(kvR, kvP);
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterCharge);
        solver->doRunSolver<Info>(context(), *info);
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
        auto* solver = solver_m.get();
        createParticles(kvR, kvP);
        const auto info = solver->createDiagnostic<Info>(Info::Kind::ScatterChargeDensity);
        solver->doRunSolver<Info>(context(), *info);
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
        auto* solver = solver_m.get();
        createParticles(kvR, kvP);
#if !OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        solver->scatterToGrid(context());
        const auto referenceSource = Kokkos::create_mirror_view_and_copy(
                Kokkos::HostSpace(), solver->getRho()->getView());
#endif
        const auto info = solver->createDiagnostic<Info>(Info::Kind::Potential);
        solver->doRunSolver<Info>(context(), *info);
        // Check the potential
#if OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        expectPotential(               //i, j, k, phi
                info->rhoView_m, {{5, 5, 2, 963538234764917},
                                     {5, 6, 2, 947474778602442},
                                     {4, 6, 2, 954504662037936},
                                     {6, 7, 2, 948706293379785}}, 10.0);
#else
        const auto potential =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), info->rhoView_m);
        for (const auto& point :
             std::array<std::array<size_t, 2>, 4>{{{5, 5}, {5, 6}, {4, 6}, {6, 7}}}) {
            const auto expected =
                    referenceGridPotential(referenceSource, point[0], point[1], 2, 0.002, 0.002);
            // Roundoff in a direct quadrature sum and FFT have different order.
            EXPECT_NEAR(potential(point[0], point[1], 2), expected, std::abs(expected) * 1e-12);
        }
#endif
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
        auto* solver = solver_m.get();
        createParticles(kvR, kvP);
#if !OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        solver->scatterToGrid(context());
        const auto referenceSource = Kokkos::create_mirror_view_and_copy(
                Kokkos::HostSpace(), solver->getRho()->getView());
#endif
        const auto info = solver->createDiagnostic<Info>(Info::Kind::EField);
        solver->doRunSolver<Info>(context(), *info);
        // Check the electric field
#if OPALX_FFT2D5_ALGORITHM_USE_IPPL_2D_POISSON_SOLVER
        expectEField(               //i, j, k, Ex,                  Ey,                   Ez
                info->eFieldView_m, {{5, 5, 2, 22636807214983296.0, 3856068046813640.0,   0.0},
                                     {5, 6, 2, -927805928809393.0, 6103969523576953.0,   0.0},
                                     {4, 6, 2, -5967249946832036.0, 12198160561896936.0,  0.0},
                                     {6, 7, 2, 12256750350007952.0, 20338717682993576.0,  0.0}});
#else
        const auto field =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), info->eFieldView_m);
        for (const auto& point :
             std::array<std::array<size_t, 2>, 4>{{{5, 5}, {5, 6}, {4, 6}, {6, 7}}}) {
            const auto expected =
                    referenceGridField(referenceSource, point[0], point[1], 2, 0.002, 0.002);
            for (unsigned axis = 0; axis < 2; ++axis) {
                EXPECT_NEAR(
                        field(point[0], point[1], 2)[axis], expected[axis],
                        std::max(1.0, std::abs(expected[axis])) * 1e-12);
            }
        }
#endif
    }

    void makeChargeDensity(FFT2D5Poisson::ScalarField2_t& chargeDensity) {
        // Fill with a delta function
        auto hostView = chargeDensity.getHostMirror();
        auto nGhost   = chargeDensity.getNghost();
        Kokkos::deep_copy(hostView, 0.0);
        hostView(nGhost, nGhost) = 1.0;
        Kokkos::deep_copy(chargeDensity.getView(), hostView);
    }

    // Independent host quadrature for the rectangular-cell potential. Unlike the
    // production corner primitive, this directly integrates the smooth kernel
    // for nonself grid displacements. The self integral is known in closed form.
    double referenceCellPotential(double x, double y, double hx, double hy) {
        if (x == 0.0 && y == 0.0) {
            const double a          = hx / 2.0;
            const double b          = hy / 2.0;
            const double averageLog = std::log(std::hypot(a, b)) - 1.5
                                      + a / (2 * b) * std::atan(b / a)
                                      + b / (2 * a) * std::atan(a / b);
            return -averageLog / (2 * Physics::pi);
        }
        static const auto quadrature = [] {
            std::array<std::array<double, 2>, 32> values{};
            for (unsigned i = 0; i < values.size(); ++i) {
                double node       = std::cos(Physics::pi * (i + 0.75) / (values.size() + 0.5));
                double derivative = 0.0;
                for (unsigned iteration = 0; iteration < 20; ++iteration) {
                    double p        = 1.0;
                    double previous = 0.0;
                    for (unsigned degree = 1; degree <= values.size(); ++degree) {
                        const double old = p;
                        p        = ((2 * degree - 1) * node * p - (degree - 1) * previous) / degree;
                        previous = old;
                    }
                    derivative        = values.size() * (node * p - previous) / (node * node - 1);
                    const double step = p / derivative;
                    node -= step;
                    if (std::abs(step) < 1e-15) {
                        break;
                    }
                }
                values[i] = {node, 2.0 / ((1 - node * node) * derivative * derivative)};
            }
            return values;
        }();
        double averageLog = 0.0;
        for (const auto& u : quadrature) {
            for (const auto& v : quadrature) {
                averageLog += 0.25 * u[1] * v[1]
                              * std::log(std::hypot(x + 0.5 * hx * u[0], y + 0.5 * hy * v[0]));
            }
        }
        return -averageLog / (2 * Physics::pi);
    }

    void verifyPoissonResults(
            const FFT2D5Poisson::ScalarField2_t& greensFn2,
            const FFT2D5Poisson::ComplexField2_t& greensTr2,
            const FFT2D5Poisson::ScalarField2_t& potential2,
            const FFT2D5Poisson::ScalarField2_t& potential,
            const FFT2D5Poisson::VectorField2_t& electric) {
        const auto& mesh          = greensFn2.get_mesh();
        const auto nx2            = static_cast<size_t>(mesh.getGridsize(0));
        const auto ny2            = static_cast<size_t>(mesh.getGridsize(1));
        const auto greens2Host    = greensFn2.getHostMirror();
        const auto greensTr2Host  = greensTr2.getHostMirror();
        const auto potential2Host = potential2.getHostMirror();
        const auto potentialHost  = potential.getHostMirror();
        const auto eHost          = electric.getHostMirror();
        Kokkos::deep_copy(greens2Host, greensFn2.getView());
        Kokkos::deep_copy(greensTr2Host, greensTr2.getView());
        Kokkos::deep_copy(potential2Host, potential2.getView());
        Kokkos::deep_copy(potentialHost, potential.getView());
        Kokkos::deep_copy(eHost, electric.getView());
        const auto nGhost2          = greensFn2.getNghost();
        const auto nGhostTr2        = greensTr2.getNghost();
        const auto nGhostPotential2 = potential2.getNghost();
        const auto nGhostPotential  = potential.getNghost();
        const auto nGhostE          = electric.getNghost();
        // Print
        if (VerboseTest) {
            for (size_t j = 0; j < ny2; ++j) {
                std::cout << j << ": ";
                for (size_t i = 0; i < nx2; ++i) {
                    std::cout << greensFn2(i + nGhost2, j + nGhost2) << " ";
                }
                std::cout << std::endl;
            }
        }
        // Check grid and contents
        ASSERT_EQ(nx2, 12);
        ASSERT_EQ(ny2, 12);
        for (size_t i = 0; i < 12; ++i) {
            for (size_t j = 0; j < 12; ++j) {
                SCOPED_TRACE(
                        ::testing::Message()
                        << "i = " << i << ", j = " << j << ", index=" << (j * 12 + i));
                ASSERT_NEAR(
                        greens2Host(i + nGhost2, j + nGhost2),
                        referenceCellPotential(std::min(i, 12 - i), std::min(j, 12 - j), 1.0, 1.0),
                        2e-13);
            }
        }
        // Check transform
        double totalGreens = 0.0;
        for (size_t j = 0; j < ny2; ++j) {
            for (size_t i = 0; i < nx2; ++i) {
                totalGreens += greens2Host(i + nGhost2, j + nGhost2);
            }
        }
        auto imagG0 = greensTr2Host(nGhostTr2, nGhostTr2).imag();
        auto realG0 = greensTr2Host(nGhostTr2, nGhostTr2).real();
        ASSERT_NEAR(totalGreens, realG0 * 144, 1e-4);
        ASSERT_NEAR(imagG0, 0.0, 1e-4);
        // Print the potential
        if (VerboseTest) {
            for (size_t j = 0; j < ny2; ++j) {
                std::cout << j << ": ";
                for (size_t i = 0; i < nx2; ++i) {
                    std::cout << potential2Host(i + nGhostPotential2, j + nGhostPotential2) << " ";
                }
                std::cout << std::endl;
            }
        }
        // Find the maximum error between the potential and the original Greens function
        double maxError = 0.0;
        for (size_t j = 0; j < ny2; ++j) {
            for (size_t i = 0; i < nx2; ++i) {
                auto diff = std::abs(
                        potential2Host(i + nGhostPotential2, j + nGhostPotential2)
                        - greens2Host(i + nGhost2, j + nGhost2));
                maxError = std::max(maxError, diff);
            }
        }
        ASSERT_NEAR(maxError, 0.0, 1e-12);
        // Now verify the potential result
        for (size_t i = 0; i < 6; ++i) {
            for (size_t j = 0; j < 6; ++j) {
                SCOPED_TRACE(
                        ::testing::Message()
                        << "i = " << i << ", j = " << j << ", index=" << (j * 6 + i));
                ASSERT_NEAR(
                        potentialHost(i + nGhostPotential, j + nGhostPotential),
                        referenceCellPotential(i, j, 1.0, 1.0), 2e-13);
            }
        }
        // Now verify the E result
        for (size_t i = 0; i < 6; ++i) {
            for (size_t j = 0; j < 6; ++j) {
                SCOPED_TRACE(
                        ::testing::Message()
                        << "i = " << i << ", j = " << j << ", index=" << (j * 6 + i));
                const auto expectedDerivative = [i, j](unsigned axis) {
                    const auto potentialAt = [axis, i, j](int cell) {
                        const double x =
                                axis == 0 ? static_cast<double>(cell) : static_cast<double>(i);
                        const double y =
                                axis == 1 ? static_cast<double>(cell) : static_cast<double>(j);
                        return referenceCellPotential(x, y, 1.0, 1.0);
                    };
                    const auto cell = static_cast<int>(axis == 0 ? i : j);
                    return (potentialAt(cell - 1) - potentialAt(cell + 1)) / 2.0;
                };
                ASSERT_NEAR(
                        eHost(i + nGhostE, j + nGhostE).data_m[0], expectedDerivative(0), 2e-13);
                ASSERT_NEAR(
                        eHost(i + nGhostE, j + nGhostE).data_m[1], expectedDerivative(1), 2e-13);
            }
        }
    }

    TEST_F(TestSolve2d5, BatchedSlicesMatchSerialPotentialHaloAndFields) {
        FFT2D5Poisson batchSolver;
        for (const auto shape : {std::array<int, 2>{6, 8}, std::array<int, 2>{8, 6}}) {
            const double dx = shape[0] == 6 ? 0.3 : 0.17;
            const double dy = shape[0] == 6 ? 0.1 : 0.29;
            const ippl::NDIndex<3> domain3{shape[0], shape[1], 73};
            FFT2D5Poisson::Mesh3_t mesh3(domain3, {dx, dy, 0.25}, {0., 0., 0.});
            ippl::FieldLayout<3> layout3(MPI_COMM_WORLD, domain3, {false, false, false});
            FFT2D5Poisson::ScalarField3_t density(mesh3, layout3), potential(mesh3, layout3);
            FFT2D5Poisson::VectorField3_t field(mesh3, layout3);
            const FFT2D5Poisson::NDIndex2_t domain2{shape[0], shape[1]};
            FFT2D5Poisson::Mesh2_t mesh2(domain2, {dx, dy}, {0., 0.});
            FFT2D5Poisson::Layout2_t layout2(MPI_COMM_WORLD, domain2, {false, false});
            FFT2D5Poisson::ScalarField2_t slice(mesh2, layout2);
            FFT2D5Poisson::VectorField2_t sliceField(mesh2, layout2);
            FFT2D5Poisson serial;
            batchSolver.initialiseViews(slice);
            serial.initialiseViews(slice);
            const auto bg = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), batchSolver.getGreensFnTr().getView());
            const auto sg = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), serial.getGreensFnTr().getView());
            ASSERT_EQ(bg.extent(0), sg.extent(0));
            ASSERT_EQ(bg.extent(1), sg.extent(1));
            for (size_t i = 0; i < bg.extent(0); ++i)
                for (size_t j = 0; j < bg.extent(1); ++j)
                    ASSERT_NEAR(Kokkos::abs(bg(i, j) - sg(i, j)), 0., 1e-13)
                            << shape[0] << " " << i << " " << j;
            const auto ghost = density.getNghost();
            // 66 active slices exercise full batches and remainders for 8/32/64.
            for (unsigned sourceCase = 0; sourceCase < 3; ++sourceCase) {
                std::vector<size_t> active;
                if (sourceCase == 0) {
                    for (size_t z = 0; z < 73; ++z)
                        if (z % 11 != 0) active.push_back(z);
                } else if (sourceCase == 1)
                    active = {2, 7, 72};
                auto source = density.getHostMirror();
                Kokkos::deep_copy(source, 0.);
                auto expectedPhi = Kokkos::create_mirror(potential.getView());
                Kokkos::deep_copy(expectedPhi, 0.);
                auto expectedE = Kokkos::create_mirror(field.getView());
                Kokkos::deep_copy(expectedE, FFT2D5Poisson::Vector3_t(0.));
                for (size_t z : active) {
                    auto rho = slice.getHostMirror();
                    Kokkos::deep_copy(rho, 0.);
                    for (int i = 0; i < shape[0]; ++i)
                        for (int j = 0; j < shape[1]; ++j) {
                            double value = (1. + 0.01 * z) * (0.3 + std::sin(0.7 * i + 0.9 * j));
                            if (z == 2)
                                value = (i == 0 && j == 0)                         ? 1.
                                        : (i == shape[0] - 1 && j == shape[1] - 1) ? -1.
                                                                                   : 0.;
                            source(i + ghost, j + ghost, z + ghost) = value;
                            rho(i + ghost, j + ghost)               = value * 1.7;
                        }
                    Kokkos::deep_copy(slice.getView(), rho);
                    serial.solve(slice, sliceField);
                    auto phi = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), slice.getView());
                    auto electric = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), sliceField.getView());
                    for (size_t i = 0; i < phi.extent(0); ++i)
                        for (size_t j = 0; j < phi.extent(1); ++j) {
                            expectedPhi(i, j, z + ghost) = phi(i, j);
                            expectedE(i, j, z + ghost)   = FFT2D5Poisson::Vector3_t(
                                    electric(i, j)[0], electric(i, j)[1], 0.);
                        }
                }
                Kokkos::deep_copy(density.getView(), source);
                auto originalSource = Kokkos::create_mirror(density.getView());
                Kokkos::deep_copy(originalSource, density.getView());
                // Begin and end with capacity 64 so the next source reuses the
                // full-batch plan while replacing its changing remainder plan.
                for (int size : {64, 1, 8, 32, 64}) {
                    SCOPED_TRACE(
                            ::testing::Message() << shape[0] << "x" << shape[1] << " batch=" << size
                                                 << " source=" << sourceCase);
                    Kokkos::deep_copy(field.getView(), FFT2D5Poisson::Vector3_t(123.));
                    Kokkos::deep_copy(potential.getView(), 123.);
                    batchSolver.solveBatched(density, field, active, size, 1.7, &potential);
                    auto actualPhi = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), potential.getView());
                    auto actualE = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), field.getView());
                    auto actualRho = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), density.getView());
                    double phiScale = 0., fieldScale = 0.;
                    for (size_t i = 0; i < actualPhi.extent(0); ++i)
                        for (size_t j = 0; j < actualPhi.extent(1); ++j)
                            for (size_t z = 0; z < actualPhi.extent(2); ++z) {
                                phiScale = std::max(phiScale, std::abs(expectedPhi(i, j, z)));
                                for (unsigned d = 0; d < 3; ++d)
                                    fieldScale =
                                            std::max(fieldScale, std::abs(expectedE(i, j, z)[d]));
                            }
                    for (size_t i = 0; i < actualPhi.extent(0); ++i)
                        for (size_t j = 0; j < actualPhi.extent(1); ++j)
                            for (size_t z = 0; z < actualPhi.extent(2); ++z) {
                                EXPECT_EQ(actualRho(i, j, z), originalSource(i, j, z));
                                EXPECT_NEAR(
                                        actualPhi(i, j, z), expectedPhi(i, j, z),
                                        phiScale * 1e-12 + 1e-30);
                                for (unsigned d = 0; d < 3; ++d)
                                    EXPECT_NEAR(
                                            actualE(i, j, z)[d], expectedE(i, j, z)[d],
                                            fieldScale * 1e-12 + 1e-30);
                            }
                }
            }
            EXPECT_THROW(
                    batchSolver.solveBatched(density, field, std::array<size_t, 1>{73}, 8),
                    OpalException);
            EXPECT_THROW(
                    batchSolver.solveBatched(density, field, std::array<size_t, 1>{1}, 0),
                    OpalException);
        }
    }

    TEST_F(TestSolve2d5, Poisson_Simple) {
        // Make a charge density
        const FFT2D5Poisson::NDIndex2_t domain{6, 6};
        FFT2D5Poisson::Mesh2_t mesh{domain, {1.0, 1.0}, {0, 0}};
        FFT2D5Poisson::Layout2_t layout{MPI_COMM_WORLD, domain, {false, false}};
        FFT2D5Poisson::ScalarField2_t rho{mesh, layout};
        makeChargeDensity(rho);
        // Make an electric field object
        FFT2D5Poisson::VectorField2_t e{mesh, layout};
        // Solve the poisson equation
        FFT2D5Poisson poisson;
        poisson.solve(rho, e);
        // Verify the Green's function
        verifyPoissonResults(
                poisson.getGreensFn(), poisson.getGreensFnTr(), poisson.getDoubledRho(), rho, e);
    }

    TEST_F(TestSolve2d5, IntegratedPoissonHasNoBoundaryCellSelfForce) {
        const FFT2D5Poisson::NDIndex2_t domain{6, 8};
        constexpr double dx = 0.3;
        constexpr double dy = 0.1;
        FFT2D5Poisson::Mesh2_t mesh{domain, {dx, dy}, {0, 0}};
        FFT2D5Poisson::Layout2_t layout{MPI_COMM_WORLD, domain, {false, false}};
        FFT2D5Poisson::ScalarField2_t rho{mesh, layout};
        FFT2D5Poisson::VectorField2_t electric{mesh, layout};
        FFT2D5Poisson poisson;
        const auto ghost = rho.getNghost();
        for (const auto& cell :
             std::array<std::array<int, 2>, 5>{{{0, 0}, {5, 0}, {0, 7}, {5, 7}, {3, 4}}}) {
            SCOPED_TRACE(::testing::Message() << "source cell=" << cell[0] << "," << cell[1]);
            auto source = rho.getHostMirror();
            Kokkos::deep_copy(source, 0.0);
            source(cell[0] + ghost, cell[1] + ghost) = 1.0;
            Kokkos::deep_copy(rho.getView(), source);
            poisson.solve(rho, electric);
            const auto field =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), electric.getView());
            const auto self = field(cell[0] + ghost, cell[1] + ghost);
            EXPECT_NEAR(self[0], 0.0, 1e-13);
            EXPECT_NEAR(self[1], 0.0, 1e-13);

            // The nearest interior probe must see the same centered free-space
            // field whether the source is in a corner or in the interior.
            const int direction   = cell[0] == 5 ? -1 : 1;
            const double expected = direction * dy / 2
                                    * (referenceCellPotential(0, 0, dx, dy)
                                       - referenceCellPotential(2 * dx, 0, dx, dy));
            EXPECT_NEAR(field(cell[0] + direction + ghost, cell[1] + ghost)[0], expected, 1e-13);
        }
    }

    TEST_F(TestSolve2d5, PoissonGradientUsesProvidedPotentialHalo) {
        const FFT2D5Poisson::NDIndex2_t domain{6, 8};
        constexpr double dx = 0.3;
        constexpr double dy = 0.1;
        FFT2D5Poisson::Mesh2_t mesh{domain, {dx, dy}, {0, 0}};
        FFT2D5Poisson::Layout2_t layout{MPI_COMM_WORLD, domain, {false, false}};
        FFT2D5Poisson::ScalarField2_t potential{mesh, layout};
        FFT2D5Poisson::VectorField2_t electric{mesh, layout};
        FFT2D5Poisson poisson;
        poisson.initialiseViews(potential);
        auto source      = potential.getHostMirror();
        const auto ghost = potential.getNghost();
        for (int i = 0; i < static_cast<int>(source.extent(0)); ++i) {
            for (int j = 0; j < static_cast<int>(source.extent(1)); ++j) {
                const double x = (i - ghost + 0.5) * dx;
                const double y = (j - ghost + 0.5) * dy;
                source(i, j)   = x * x + 3 * x * y + 2 * y * y;
            }
        }
        Kokkos::deep_copy(potential.getView(), source);
        poisson.determineField(potential, electric);
        const auto field =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), electric.getView());
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 8; ++j) {
                const double x = (i + 0.5) * dx;
                const double y = (j + 0.5) * dy;
                EXPECT_NEAR(field(i + ghost, j + ghost)[0], -(2 * x + 3 * y), 1e-12);
                EXPECT_NEAR(field(i + ghost, j + ghost)[1], -(3 * x + 4 * y), 1e-12);
            }
        }
    }

    TEST_F(TestSolve2d5, IntegratedPoissonSupportsSmallTransverseGrids) {
        for (const auto& dimensions :
             std::array<std::array<unsigned, 2>, 4>{{{1, 1}, {1, 2}, {2, 1}, {2, 2}}}) {
            SCOPED_TRACE(::testing::Message() << "grid=" << dimensions[0] << "x" << dimensions[1]);
            const FFT2D5Poisson::NDIndex2_t domain{dimensions[0], dimensions[1]};
            FFT2D5Poisson::Mesh2_t mesh{domain, {0.3, 0.1}, {0, 0}};
            FFT2D5Poisson::Layout2_t layout{MPI_COMM_WORLD, domain, {false, false}};
            FFT2D5Poisson::ScalarField2_t rho{mesh, layout};
            FFT2D5Poisson::VectorField2_t electric{mesh, layout};
            auto source      = rho.getHostMirror();
            const auto ghost = rho.getNghost();
            Kokkos::deep_copy(source, 0.0);
            source(ghost, ghost) = 1.0;
            Kokkos::deep_copy(rho.getView(), source);
            FFT2D5Poisson poisson;
            poisson.solve(rho, electric);
            const auto field =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), electric.getView());
            EXPECT_NEAR(field(ghost, ghost)[0], 0.0, 1e-13);
            EXPECT_NEAR(field(ghost, ghost)[1], 0.0, 1e-13);
            for (unsigned i = 0; i < dimensions[0]; ++i) {
                for (unsigned j = 0; j < dimensions[1]; ++j) {
                    EXPECT_TRUE(std::isfinite(field(i + ghost, j + ghost)[0]));
                    EXPECT_TRUE(std::isfinite(field(i + ghost, j + ghost)[1]));
                }
            }
        }
    }

    TEST_F(TestSolve2d5, IntegratedLogKernelSelfAndAnisotropicCells) {
        for (const auto& widths : std::array<std::array<double, 2>, 4>{
                     {{1.0, 1.0}, {0.000375, 0.00009375}, {0.02, 0.001}, {0.001, 0.02}}}) {
            const auto hx = widths[0];
            const auto hy = widths[1];
            EXPECT_NEAR(
                    FFT2D5Poisson::cellAveragedGreensFunction(0, 0, hx, hy),
                    referenceCellPotential(0, 0, hx, hy), 3e-15);
            for (const auto& offset : std::array<std::array<double, 2>, 4>{
                         {{2.0, 1.0}, {-2.0, 1.0}, {8.0, 0.0}, {10000.0, -3000.0}}}) {
                const auto scale = std::max(hx, hy);
                const auto x     = offset[0] * scale;
                const auto y     = offset[1] * scale;
                EXPECT_NEAR(
                        FFT2D5Poisson::cellAveragedGreensFunction(x, y, hx, hy),
                        referenceCellPotential(x, y, hx, hy), 3e-13);
                EXPECT_NEAR(
                        FFT2D5Poisson::cellAveragedGreensFunction(x, y, hx, hy),
                        FFT2D5Poisson::cellAveragedGreensFunction(y, x, hy, hx), 3e-13);
            }
        }
    }

    TEST_F(TestSolve2d5, IntegratedLogKernelScaleChangesOnlyPotentialGauge) {
        for (const auto& offset : std::array<std::array<double, 2>, 4>{
                     {{0.0, 0.0}, {1.0, 0.0}, {2.0, 3.0}, {10000.0, 0.0}}}) {
            const auto original =
                    FFT2D5Poisson::cellAveragedGreensFunction(offset[0], offset[1], 0.5, 0.125);
            for (const double scale : {1e-6, 1e-3, 10.0, 1e6}) {
                EXPECT_NEAR(
                        FFT2D5Poisson::cellAveragedGreensFunction(
                                scale * offset[0], scale * offset[1], scale * 0.5, scale * 0.125),
                        original - std::log(scale) / (2 * Physics::pi), 3e-14);
            }
        }
    }

    double gaussianPoissonForceError(unsigned nx, unsigned ny, double sigma) {
        const FFT2D5Poisson::NDIndex2_t domain{nx, ny};
        const double length = 12.0 * sigma;
        const double dx     = length / nx;
        const double dy     = length / ny;
        FFT2D5Poisson::Mesh2_t mesh{domain, {dx, dy}, {-length / 2, -length / 2}};
        FFT2D5Poisson::Layout2_t layout{MPI_COMM_WORLD, domain, {false, false}};
        FFT2D5Poisson::ScalarField2_t rho{mesh, layout};
        FFT2D5Poisson::VectorField2_t electric{mesh, layout};
        auto source      = rho.getHostMirror();
        const auto ghost = rho.getNghost();
        Kokkos::deep_copy(source, 0.0);
        for (unsigned i = 0; i < nx; ++i) {
            for (unsigned j = 0; j < ny; ++j) {
                const auto x                 = (i + 0.5) * dx - length / 2;
                const auto y                 = (j + 0.5) * dy - length / 2;
                source(i + ghost, j + ghost) = std::exp(-(x * x + y * y) / (2 * sigma * sigma))
                                               / (2 * Physics::pi * sigma * sigma);
            }
        }
        Kokkos::deep_copy(rho.getView(), source);
        FFT2D5Poisson poisson;
        poisson.solve(rho, electric);
        auto field = electric.getHostMirror();
        Kokkos::deep_copy(field, electric.getView());
        double maxRelativeError = 0.0;
        for (unsigned i = 0; i < nx; ++i) {
            for (unsigned j = 0; j < ny; ++j) {
                const auto x  = (i + 0.5) * dx - length / 2;
                const auto y  = (j + 0.5) * dy - length / 2;
                const auto r2 = x * x + y * y;
                if (r2 < 0.04 * sigma * sigma || r2 > 6.25 * sigma * sigma) {
                    continue;
                }
                // Gauss's law for a unit-charge round Gaussian, with epsilon=1.
                const auto coefficient =
                        -std::expm1(-r2 / (2 * sigma * sigma)) / (2 * Physics::pi * r2);
                const auto ex    = coefficient * x;
                const auto ey    = coefficient * y;
                const auto value = field(i + ghost, j + ghost);
                const auto error = std::hypot(value[0] - ex, value[1] - ey) / std::hypot(ex, ey);
                maxRelativeError = std::max(maxRelativeError, error);
            }
        }
        return maxRelativeError;
    }

    TEST_F(TestSolve2d5, IntegratedPoissonForceConvergesAndIsIndependentOfLengthUnit) {
        const auto coarse = gaussianPoissonForceError(64, 128, 0.001);
        const auto fine   = gaussianPoissonForceError(128, 256, 0.001);
        EXPECT_LT(fine, 0.0035);
        EXPECT_LT(fine, coarse / 3.5);
        EXPECT_NEAR(fine, gaussianPoissonForceError(128, 256, 1.0), 2e-11);
    }
}  // namespace
