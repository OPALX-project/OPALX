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
#include "PartBunch/PartBunch.h"
#include "PartBunch/Solve2d5.h"
#include "PartBunch/SubField.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

namespace {

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

    class FrenetSerretDiagnostic : public ParticleBaseDiagnostic {
    public:
        explicit FrenetSerretDiagnostic(const ParticleContainer_t* pc)
            : ParticleBaseDiagnostic(pc) {}

        KOKKOS_FUNCTION void frenetSerret(
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

    class ScatterDiagnostic : public Solve2d5_t::NullDiagnostic {
    public:
        explicit ScatterDiagnostic(const Solve2d5_t::ScalarGridView3D_t& rho) : NullDiagnostic() {
            Kokkos::resize(rhoView, rho.extent(0), rho.extent(1), rho.extent(2));
        }

        KOKKOS_FUNCTION void scatterCharge(const Solve2d5_t::ScalarGridView3D_t& rho) const {
            Kokkos::deep_copy(rhoView, rho);
        }

        Field_t<3>::view_type rhoView;
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

        void createParticles(size_t nPart, double pzMin, double pzMax) {
            pc->createParticles(nPart);

            std::mt19937_64 eng(42 + ippl::Comm->rank());

            const auto R_host  = pc->R.getHostMirror();
            const auto P_host  = pc->P.getHostMirror();
            const auto dt_host = pc->dt.getHostMirror();
            const auto E_host  = pc->E.getHostMirror();

            // Match mesh extents from the bunch' field container.
            const auto rmin = bunch->rmin_m;
            const auto rmax = bunch->rmax_m;

            std::uniform_real_distribution unifR_x(rmin[0] + 0.1, rmax[0] - 0.1);
            std::uniform_real_distribution unifR_y(rmin[1] + 0.1, rmax[1] - 0.1);
            std::uniform_real_distribution unifR_z(rmin[2] + 0.1, rmax[2] - 0.1);
            std::uniform_real_distribution unifP_z(pzMin, pzMax);

            const double dt = bunch->getdT();
            const double qi = pc->getChargePerParticle();

            for (size_t i = 0; i < nPart; ++i) {
                R_host(i)[0] = unifR_x(eng);
                R_host(i)[1] = unifR_y(eng);
                R_host(i)[2] = unifR_z(eng);

                P_host(i)[0] = 0.0;
                P_host(i)[1] = 0.0;
                P_host(i)[2] = unifP_z(eng);

                dt_host(i) = dt;

                // Initialize particle E to zero; solver should leave it zero in NONE mode.
                E_host(i)[0] = 0.0;
                E_host(i)[1] = 0.0;
                E_host(i)[2] = 0.0;
            }

            Kokkos::deep_copy(pc->R.getView(), R_host);
            Kokkos::deep_copy(pc->P.getView(), P_host);
            Kokkos::deep_copy(pc->dt.getView(), dt_host);
            Kokkos::deep_copy(pc->E.getView(), E_host);
            pc->setQ(qi);

            ippl::Comm->barrier();
            Kokkos::fence();
            pc->update();
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

        void expectAllParticleEZeroAndFinite(const double tol) {
            const auto E_host = pc->E.getHostMirror();
            Kokkos::deep_copy(E_host, pc->E.getView());

            const size_t localN = pc->getLocalNum();
            for (size_t i = 0; i < localN; ++i) {
                for (unsigned d = 0; d < 3; ++d) {
                    const double val = E_host(i)[d];
                    EXPECT_TRUE(isfinite(val)) << "E[" << i << "][" << d << "]=" << val;
                    EXPECT_NEAR(val, 0.0, tol) << "E[" << i << "][" << d << "]";
                }
            }
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
            for (size_t k = 0; k < rho.extent(2); ++k) {
                for (size_t j = 0; j < rho.extent(1); ++j) {
                    std::cout << k << "," << j << ": ";
                    for (size_t i = 0; i < rho.extent(0); ++i) {
                        std::cout << hostView(i, j, k) << " ";
                    }
                    std::cout << std::endl;
                }
            }
            // Check values
            for (const auto& e : expected) {
                SCOPED_TRACE(std::format("Index: {},{},{}", e.i, e.j, e.k));
                EXPECT_NEAR(hostView(e.i, e.j, e.k), e.value, 1e-6);
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
            for (size_t k = 0; k < hostView.extent(2); ++k) {
                for (size_t j = 0; j < hostView.extent(1); ++j) {
                    std::cout << k << "," << j << ": ";
                    for (size_t i = 0; i < hostView.extent(0); ++i) {
                        std::cout << hostView(i, j, k) << " ";
                    }
                    std::cout << std::endl;
                }
            }
            // Check values
            for (const auto& e : expected) {
                SCOPED_TRACE(std::format("Index: {},{},{}", e.i, e.j, e.k));
                EXPECT_NEAR(hostView(e.i, e.j, e.k)[0], e.valueX, 1e-6);
                EXPECT_NEAR(hostView(e.i, e.j, e.k)[1], e.valueY, 1e-6);
                EXPECT_NEAR(hostView(e.i, e.j, e.k)[2], e.valueZ, 1e-6);
            }
        }

        std::shared_ptr<TestableFieldSolverCmd> fsCmd;
        std::shared_ptr<DataSink> dataSink;
        std::shared_ptr<Beam> beam;
        std::shared_ptr<PartBunch_t> bunch;
        std::shared_ptr<ParticleContainer_t> pc;
    };

    TEST_F(TestSolve2d5, FieldSolverCmdTypes) {
        fsCmd->setType("OPEN2D5");
        EXPECT_NO_THROW(rebuildBunch());
        EXPECT_TRUE(dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver()) != nullptr);
        fsCmd->setType("CIRCULAR2D5");
        EXPECT_NO_THROW(rebuildBunch());
        EXPECT_TRUE(dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver()) != nullptr);
        fsCmd->setType("PLATES2D5");
        EXPECT_NO_THROW(rebuildBunch());
        EXPECT_TRUE(dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver()) != nullptr);
    }

    TEST_F(TestSolve2d5, SliceSolverSetup) {
        fsCmd->setType("OPEN2D5");
        rebuildBunch();
        const auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        EXPECT_EQ(solver->getNumSlices(), nz);
        ASSERT_FALSE(solver->getSliceMesh() == nullptr);
        EXPECT_EQ(solver->getSliceMesh()->getGridsize(0), 8);
        EXPECT_EQ(solver->getSliceMesh()->getGridsize(1), 9);
        EXPECT_EQ(solver->getSliceMesh()->getOrigin()[0], -3);
        EXPECT_EQ(solver->getSliceMesh()->getOrigin()[1], -3);
        EXPECT_DOUBLE_EQ(solver->getSliceMesh()->getMeshSpacing()[0], 0.75);
        EXPECT_NEAR(solver->getSliceMesh()->getMeshSpacing()[1], 0.6666667, 1e-6);
        ASSERT_FALSE(solver->getSliceLayout() == nullptr);
        EXPECT_EQ(solver->getSliceLayout()->getDomain()[0], 8);
        EXPECT_EQ(solver->getSliceLayout()->getDomain()[1], 9);
    }

    TEST_F(TestSolve2d5, LoadReferencePath_Missing) {
        fsCmd->setType("OPEN2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        EXPECT_ANY_THROW(solver->loadReferencePath());
    }

    TEST_F(TestSolve2d5, LoadReferencePath_ReadFail) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}},
                true);
        fsCmd->setType("OPEN2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        EXPECT_ANY_THROW(solver->loadReferencePath());
    }

    TEST_F(TestSolve2d5, LoadReferencePath_Success) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}});
        fsCmd->setType("OPEN2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        ASSERT_NO_THROW(solver->loadReferencePath());
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

    TEST_F(TestSolve2d5, ToFrenetSerret_NoParticles) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}});
        fsCmd->setType("OPEN2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        solver->loadReferencePath();
        EXPECT_NO_THROW(solver->scatterToGrid(*bunch));
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_NoReferencePath) {
        fsCmd->setType("OPEN2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{1, 2, 3}}, {{4, 5, 6}});
        const FrenetSerretDiagnostic info(pc.get());
        solver->scatterToGrid<FrenetSerretDiagnostic>(*bunch, info);
        auto [r, p] = getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {1, 2, 3}, {4, 5, 6});
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_ShortReferencePath) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}});
        fsCmd->setType("OPEN2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{1, 2, 3}}, {{4, 5, 6}});
        solver->loadReferencePath();
        const FrenetSerretDiagnostic info(pc.get());
        solver->scatterToGrid<FrenetSerretDiagnostic>(*bunch, info);
        auto [r, p] = getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {1, 2, 3}, {4, 5, 6});
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_TrivialReferencePath) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}});
        fsCmd->setType("OPEN2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 0}, {0, 0, 0.5}, {0, 0, 1}}, {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}});
        solver->loadReferencePath();
        const FrenetSerretDiagnostic info(pc.get());
        solver->scatterToGrid<FrenetSerretDiagnostic>(*bunch, info);
        auto [r, p] = info.getParticles();
        ASSERT_EQ(r.size(), 3);
        ASSERT_EQ(p.size(), 3);
        expectParticle(0, r, p, {0, 0, 0}, {0, 0, 1});
        expectParticle(1, r, p, {0, 0, 0.5}, {0, 0, 1});
        expectParticle(2, r, p, {0, 0, 1}, {0, 0, 1});
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}});
        fsCmd->setType("OPEN2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 2}}, {{0, 0, 1}});
        solver->loadReferencePath();
        const FrenetSerretDiagnostic info(pc.get());
        solver->scatterToGrid<FrenetSerretDiagnostic>(*bunch, info);
        auto [r, p] = info.getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {-0.7071068, 0, 1.7071068}, {-0.7071068, 0, 0.7071068});
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_AtRefCorner) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 1}, {1, 0, 2}, {0, 0, 3}});
        fsCmd->setType("OPEN2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{2, 0, 2}}, {{0, 0, 1}});
        solver->loadReferencePath();
        const FrenetSerretDiagnostic info(pc.get());
        solver->scatterToGrid<FrenetSerretDiagnostic>(*bunch, info);
        auto [r, p] = info.getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {0.7071068, 0, 2.4142136}, {-0.7071068, 0, 0.7071068});
    }

    TEST_F(TestSolve2d5, ToFrenetSerret_DegenerateRefPath) {
        makeReferencePathFile(
                "data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 0}, {0, 0, 1}, {1, 0, 2}});
        fsCmd->setType("OPEN2D5");
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 2}}, {{0, 0, 1}});
        solver->loadReferencePath();
        const FrenetSerretDiagnostic info(pc.get());
        solver->scatterToGrid<FrenetSerretDiagnostic>(*bunch, info);
        auto [r, p] = info.getParticles();
        ASSERT_EQ(r.size(), 1);
        ASSERT_EQ(p.size(), 1);
        expectParticle(0, r, p, {-0.7071068, 0, 1.7071068}, {-0.7071068, 0, 0.7071068});
    }

    TEST_F(TestSolve2d5, BoostToBeamFrame_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("OPEN2D5");
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

    TEST_F(TestSolve2d5, ScatterToGrid_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("OPEN2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 0}}, {{0, 0, 0}});
        solver->loadReferencePath();
        const ScatterDiagnostic info(solver->getRho()->getView());
        solver->scatterToGrid<ScatterDiagnostic>(*bunch, info);
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

    TEST_F(TestSolve2d5, SolvePoissons_Simple) {
        makeReferencePathFile("data/unit_test_DesignPath.dat", {{0, 0, 0}, {0, 0, 3}});
        fsCmd->setType("OPEN2D5");
        fsCmd->setNX(12);
        fsCmd->setNY(12);
        fsCmd->setNZ(12);
        rebuildBunch();
        auto* solver = dynamic_cast<Solve2d5_t*>(bunch->getFieldSolver());
        createParticles({{0, 0, 0}}, {{0, 0, 0}});
        solver->loadReferencePath();
        solver->scatterToGrid(*bunch);
        solver->solvePoissons();
        expectEField(
                solver->getEField()->getView(), {{1, 1, 6, -0.000150, -0.000150, 0},
                                                 {2, 1, 6, -0.000148, -0.000180, 0},
                                                 {4, 4, 6, -0.000333, -0.000333, 0},
                                                 {4, 8, 6, -0.000487, 0.000291, 0}});
    }

}  // namespace
