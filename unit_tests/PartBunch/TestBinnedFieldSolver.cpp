/**
 * @file TestBinnedFieldSolver.cpp
 * @brief Smoke tests for binned self-fields and open/periodic P3M configuration.
 *
 * This file validates that the self-field computation pathway is stable across the
 * "legacy" (no binning attached) and "binned" (adaptive bins attached) execution paths.
 *
 * The tests construct a small `PartBunch<double,3>` with a minimal `FieldSolverCmd`
 * configuration (type `"NONE"` and periodic FFT boundary conditions), populate a set of
 * particles with deterministic random initial conditions, and then invoke
 * `PartBunch::computeSelfFields()`.
 *
 * Key behaviors verified:
 * - `computeSelfFields()` does not throw for a non-binned bunch.
 * - `computeSelfFields()` does not throw when adaptive binning is attached.
 * - When using solver type `"NONE"`, the per-particle electric field `E` remains finite
 *   and (near) zero after the call.
 * - When binning is active, the current bin count is sane (between 1 and the configured
 *   maximum).
 * - Open and periodic P3M select the same solver wrapper with matching particle boundaries.
 * - Mixed and unsupported P3M boundary conditions are rejected during command validation.
 *
 * Notes:
 * - The fixture initializes IPPL and disables HDF5 output to keep the tests lightweight.
 * - The goal is a robustness/smoke check (no physics validation of non-trivial fields, since this
 * would require way more computational resources).
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <variant>

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "PartBunch/BinnedFieldSolver.h"
#include "PartBunch/FieldMirror.hpp"
#include "PartBunch/PartBunch.h"
#include "Structure/Beam.h"
#include "Structure/BinningCmd.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

namespace {

    class TestableFieldSolverCmd : public FieldSolverCmd {
    public:
        void setType(const std::string& t) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::TYPE], t);
        }

        void setBCX(const std::string& bc) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTX], bc);
        }
        void setBCY(const std::string& bc) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTY], bc);
        }
        void setBCZ(const std::string& bc) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTZ], bc);
        }

        void setGreensFunction(const std::string& greensFunction) {
            Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::GREENSF], greensFunction);
        }

        void setParallelDecomposition(bool value) {
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTX], value);
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTY], value);
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTZ], value);
        }

        void setBinsName(const std::string& binsName) {
            Attributes::setString(this->itsAttr[FIELDSOLVER::BINS], binsName);
        }

        void setP3MCutoff(double cutoff) {
            Attributes::setReal(this->itsAttr[FIELDSOLVER::P3MRCUT], cutoff);
        }

    };

    class TestableBinningCmd : public BinningCmd {
    public:
        void setMaxBins(int value) {
            Attributes::setReal(itsAttr[BINNING::MAXBINS], static_cast<double>(value));
        }

        void setAdaptiveBinning(bool value) {
            Attributes::setBool(itsAttr[BINNING::ADAPTIVEBINNING], value);
        }

        void setTablePrintFrequency(double value) {
            Attributes::setReal(itsAttr[BINNING::TABLEPRINTFREQ], value);
        }

        void setParameterString(const std::string& value) {
            Attributes::setPredefinedString(itsAttr[BINNING::PARAMETER], value);
        }
    };

    using PartBunch_t          = PartBunch<double, 3>;
    using ParticleContainer_t  = typename PartBunch_t::ParticleContainer_t;
    using AdaptBins_t          = typename PartBunch_t::AdaptBins_t;
    using CoordinateSelector_t = typename PartBunch_t::CoordinateSelector_t;

    constexpr size_t kDefaultNParticles = 64;

    class BinnedFieldSolverSmokeTest : public ::testing::Test {
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
            // test doesn't have an H5PartWrapper. Disable HDF5 for this smoke test.
            Options::enableHDF5 = false;
        }

        static void TearDownTestSuite() {
            delete gmsg;
            gmsg = nullptr;
            ippl::finalize();
        }

        void SetUp() override {
            const auto* testInfo = ::testing::UnitTest::GetInstance()->current_test_info();
            if (testInfo != nullptr
                && std::string(testInfo->name()) == "MirrorFieldZHandlesNonSlab3DDecomposition") {
                return;
            }

            // Keep mesh small so scatter/solve/gather are quick.
            constexpr double nx = 8;
            constexpr double ny = 8;
            constexpr double nz = 8;

            fsCmd = std::make_shared<TestableFieldSolverCmd>();
            fsCmd->setType("NONE");
            fsCmd->setNX(nx);
            fsCmd->setNY(ny);
            fsCmd->setNZ(nz);
            fsCmd->setBCX("PERIODIC");
            fsCmd->setBCY("PERIODIC");
            fsCmd->setBCZ("PERIODIC");
            fsCmd->setGreensFunction("STANDARD");
            fsCmd->setParallelDecomposition(true);
            fsCmd->execute();

            // Keep the concrete solver command alive; PartBunch borrows it.
            fsCmdBase = fsCmd;

            dataSink       = std::make_shared<DataSink>();
            beam           = std::make_shared<Beam>();
            Beam* testBeam = Beam::find("UNNAMED_BEAM");

            // qi/mi/lbt are used by rho scaling; but with NullSolver we mostly validate
            // "no-throw" and deterministic zero E behavior.
            bunch = std::make_shared<PartBunch_t>(
                    /*qi=*/std::vector<double>{1.0},
                    /*mi=*/std::vector<double>{1.0},
                    /*beams=*/std::vector<Beam*>{testBeam},
                    /*totalParticlesPerBeam=*/std::vector<size_t>{kDefaultNParticles},
                    /*lbt=*/1.0,
                    /*integration_method=*/"LF2", fsCmdBase.get(), dataSink.get());
            pc = bunch->getParticleContainer();
        }

        void TearDown() override {
            // Ensure device allocations can be freed between tests.
            bunch.reset();
            dataSink.reset();
            fsCmd.reset();
            fsCmdBase.reset();
            pc.reset();
        }

        void createParticles(size_t nPart, double pzMin, double pzMax) {
            pc->createParticles(nPart);

            std::mt19937_64 eng(42 + ippl::Comm->rank());

            auto R_host  = pc->R.getHostMirror();
            auto P_host  = pc->P.getHostMirror();
            auto dt_host = pc->dt.getHostMirror();
            auto E_host  = pc->E.getHostMirror();

            // Match mesh extents from the bunch' field container.
            const auto rmin = bunch->rmin_m;
            const auto rmax = bunch->rmax_m;

            std::uniform_real_distribution<double> unifR_x(rmin[0] + 0.1, rmax[0] - 0.1);
            std::uniform_real_distribution<double> unifR_y(rmin[1] + 0.1, rmax[1] - 0.1);
            std::uniform_real_distribution<double> unifR_z(rmin[2] + 0.1, rmax[2] - 0.1);
            std::uniform_real_distribution<double> unifP_z(pzMin, pzMax);

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

        void createZMirrorSymmetricParticles(size_t nPart, double pz) {
            ASSERT_EQ(nPart % 2, 0u);
            pc->createParticles(nPart);

            auto R_host  = pc->R.getHostMirror();
            auto P_host  = pc->P.getHostMirror();
            auto dt_host = pc->dt.getHostMirror();
            auto E_host  = pc->E.getHostMirror();
            auto B_host  = pc->B.getHostMirror();

            const auto rmin                  = bunch->rmin_m;
            const auto rmax                  = bunch->rmax_m;
            const Vector_t<double, 3> center = 0.5 * (rmin + rmax);
            const Vector_t<double, 3> width  = rmax - rmin;
            const double dt                  = bunch->getdT();
            const double qi                  = pc->getChargePerParticle();

            for (size_t pair = 0; pair < nPart / 2; ++pair) {
                const double phase = static_cast<double>(pair + 1);
                const double x     = center[0] + 0.18 * width[0] * std::sin(0.7 * phase);
                const double y     = center[1] + 0.16 * width[1] * std::cos(1.1 * phase);
                const double z     = 0.08 * width[2] * (0.25 + static_cast<double>(pair % 7) / 7.0);

                for (size_t side = 0; side < 2; ++side) {
                    const size_t i = 2 * pair + side;
                    R_host(i)[0]   = x;
                    R_host(i)[1]   = y;
                    R_host(i)[2]   = center[2] + (side == 0 ? -z : z);
                    P_host(i)[0]   = 0.0;
                    P_host(i)[1]   = 0.0;
                    P_host(i)[2]   = pz;
                    dt_host(i)     = dt;
                    E_host(i)      = Vector_t<double, 3>(0.0);
                    B_host(i)      = Vector_t<double, 3>(0.0);
                }
            }

            Kokkos::deep_copy(pc->R.getView(), R_host);
            Kokkos::deep_copy(pc->P.getView(), P_host);
            Kokkos::deep_copy(pc->dt.getView(), dt_host);
            Kokkos::deep_copy(pc->E.getView(), E_host);
            Kokkos::deep_copy(pc->B.getView(), B_host);
            pc->setQ(qi);

            ippl::Comm->barrier();
            Kokkos::fence();
            pc->update();
        }

        void rebuildBunch() {
            fsCmd->execute();
            fsCmdBase      = fsCmd;
            Beam* testBeam = Beam::find("UNNAMED_BEAM");
            bunch          = std::make_shared<PartBunch_t>(
                    /*qi=*/std::vector<double>{1.0},
                    /*mi=*/std::vector<double>{1.0},
                    /*beams=*/std::vector<Beam*>{testBeam},
                    /*totalParticlesPerBeam=*/std::vector<size_t>{kDefaultNParticles},
                    /*lbt=*/1.0,
                    /*integration_method=*/"LF2", fsCmdBase.get(), dataSink.get());
            pc = bunch->getParticleContainer();
        }

        void rebuildOpenBunchWithGreensFunction(const std::string& greensFunction) {
            fsCmd->setType("OPEN");
            fsCmd->setBCX("OPEN");
            fsCmd->setBCY("OPEN");
            fsCmd->setBCZ("OPEN");
            fsCmd->setGreensFunction(greensFunction);
            rebuildBunch();
        }

        void rebuildP3MBunch(const std::string& boundary) {
            fsCmd->setType("P3M");
            fsCmd->setP3MCutoff(0.25);
            fsCmd->setBCX(boundary);
            fsCmd->setBCY(boundary);
            fsCmd->setBCZ(boundary);
            fsCmd->setParallelDecomposition(true);
            fsCmd->execute();
            rebuildBunch();
        }

        std::shared_ptr<AdaptBins_t> attachBins(
                typename AdaptBins_t::bin_index_type maxBins, double alpha, double beta,
                double desiredWidth) {
            using ConcreteBins_t =
                    ParticleBinning::AdaptBins<ParticleContainer_t, CoordinateSelector_t>;

            CoordinateSelector_t selector(/*axis=*/2);
            auto bins = std::make_shared<ConcreteBins_t>(
                    *pc, selector, maxBins, alpha, beta, desiredWidth,
                    /*binningCmdName=*/"TEST_BINNING_CMD");
            bunch->setBins(bins);
            return bins;
        }

        void defineBinningCommand(
                const std::string& name, typename AdaptBins_t::bin_index_type maxBins,
                bool adaptiveBinning) {
            auto* binsCmd = new TestableBinningCmd();
            binsCmd->setOpalName(name);
            binsCmd->setMaxBins(maxBins);
            binsCmd->setAdaptiveBinning(adaptiveBinning);
            binsCmd->setTablePrintFrequency(0.0);
            binsCmd->setParameterString("VELOCITYZ");
            binsCmd->execute();
            OpalData::getInstance()->define(binsCmd);
        }

        void expectAllParticleEZeroAndFinite(double tol) {
            auto E_host = pc->E.getHostMirror();
            Kokkos::deep_copy(E_host, pc->E.getView());

            const size_t localN = pc->getLocalNum();
            for (size_t i = 0; i < localN; ++i) {
                for (unsigned d = 0; d < 3; ++d) {
                    const double val = E_host(i)[d];
                    EXPECT_TRUE(std::isfinite(val)) << "E[" << i << "][" << d << "]=" << val;
                    EXPECT_NEAR(val, 0.0, tol) << "E[" << i << "][" << d << "]";
                }
            }
        }

        std::shared_ptr<TestableFieldSolverCmd> fsCmd;
        std::shared_ptr<FieldSolverCmd> fsCmdBase;
        std::shared_ptr<DataSink> dataSink;
        std::shared_ptr<Beam> beam;
        std::shared_ptr<PartBunch_t> bunch;
        std::shared_ptr<ParticleContainer_t> pc;
    };

    TEST_F(BinnedFieldSolverSmokeTest, LegacyPath_NoBins_NoThrowAndEZero) {
        ASSERT_FALSE(bunch->hasBinning());
        createParticles(kDefaultNParticles, /*pzMin=*/0.1, /*pzMax=*/0.9);

        EXPECT_NO_THROW(bunch->computeSelfFields());
        expectAllParticleEZeroAndFinite(/*tol=*/1e-8);
    }

    TEST_F(BinnedFieldSolverSmokeTest, BinnedPath_WithBins_NoThrowAndEZero) {
        createParticles(kDefaultNParticles, /*pzMin=*/0.1, /*pzMax=*/2.0);

        constexpr AdaptBins_t::bin_index_type maxBins = 6;
        auto bins = attachBins(maxBins, /*alpha=*/1.0, /*beta=*/1.0, /*desiredWidth=*/0.3);
        ASSERT_TRUE(bunch->hasBinning());

        EXPECT_NO_THROW(bunch->computeSelfFields());

        const auto currentBins = bins->getCurrentBinCount();
        EXPECT_GE(currentBins, 1);
        EXPECT_LE(currentBins, maxBins);

        expectAllParticleEZeroAndFinite(/*tol=*/1e-8);
    }

    TEST_F(BinnedFieldSolverSmokeTest, BinnedPath_AdaptiveBinningFalseKeepsUniformMaxBins) {
        constexpr AdaptBins_t::bin_index_type maxBins = 6;
        const std::string binsName                    = "UNIT_TEST_UNIFORM_BINNING";
        defineBinningCommand(binsName, maxBins, /*adaptiveBinning=*/false);

        fsCmd->setType("OPEN");
        fsCmd->setBinsName(binsName);
        fsCmd->setBCX("OPEN");
        fsCmd->setBCY("OPEN");
        fsCmd->setBCZ("OPEN");
        rebuildBunch();

        createParticles(kDefaultNParticles, /*pzMin=*/0.1, /*pzMax=*/2.0);
        ASSERT_TRUE(bunch->hasBinning());

        EXPECT_NO_THROW(bunch->computeSelfFields());

        auto bins = bunch->getBins();
        ASSERT_NE(bins, nullptr);
        EXPECT_EQ(bins->getCurrentBinCount(), maxBins);
    }

    TEST_F(BinnedFieldSolverSmokeTest, OpenSolver_UsesStandardGreensFunctionFromFieldSolverCmd) {
        ASSERT_NO_THROW(rebuildOpenBunchWithGreensFunction("STANDARD"));

        ASSERT_NE(bunch->getFieldSolver(), nullptr);
        EXPECT_EQ(bunch->getFieldSolver()->getGreensFunction(), "STANDARD");
    }

    TEST_F(BinnedFieldSolverSmokeTest, OpenSolver_UsesIntegratedGreensFunctionFromFieldSolverCmd) {
        ASSERT_NO_THROW(rebuildOpenBunchWithGreensFunction("INTEGRATED"));

        ASSERT_NE(bunch->getFieldSolver(), nullptr);
        EXPECT_EQ(bunch->getFieldSolver()->getGreensFunction(), "INTEGRATED");
    }

    TEST_F(BinnedFieldSolverSmokeTest, BeamBeamCopyDoublesEAndCancelsBAtExactOverlap) {
        rebuildOpenBunchWithGreensFunction("INTEGRATED");
        createZMirrorSymmetricParticles(kDefaultNParticles, /*pz=*/2.0);
        attachBins(/*maxBins=*/1, /*alpha=*/1.0, /*beta=*/1.0, /*desiredWidth=*/1.0);

        const double interactionPointS =
                pc->get_sPos() + 0.5 * (bunch->rmin_m[2] + bunch->rmax_m[2]);
        bunch->setBeamBeamWindowConfig(
                bunch->rmax_m[2] - bunch->rmin_m[2], interactionPointS, bunch->rmin_m[2],
                bunch->rmax_m[2], /*copyModel=*/false);
        ASSERT_NO_THROW(bunch->computeSelfFields());
        ASSERT_TRUE(bunch->hasLastDepositedChargeBeforeBackground());
        const double physicalCharge = bunch->getLastDepositedChargeBeforeBackground();

        auto physicalE = pc->E.getHostMirror();
        auto physicalB = pc->B.getHostMirror();
        auto originalR = pc->R.getHostMirror();
        Kokkos::deep_copy(physicalE, pc->E.getView());
        Kokkos::deep_copy(physicalB, pc->B.getView());
        Kokkos::deep_copy(originalR, pc->R.getView());

        bunch->setBeamBeamWindowConfig(
                bunch->rmax_m[2] - bunch->rmin_m[2], interactionPointS, bunch->rmin_m[2],
                bunch->rmax_m[2], /*copyModel=*/true);
        ASSERT_NO_THROW(bunch->computeSelfFields());
        ASSERT_TRUE(bunch->hasLastDepositedChargeBeforeBackground());
        EXPECT_NEAR(
                bunch->getLastDepositedChargeBeforeBackground(), 2.0 * physicalCharge,
                std::max(1.0e-18, 1.0e-12 * std::abs(physicalCharge)));

        auto copiedE   = pc->E.getHostMirror();
        auto copiedB   = pc->B.getHostMirror();
        auto restoredR = pc->R.getHostMirror();
        Kokkos::deep_copy(copiedE, pc->E.getView());
        Kokkos::deep_copy(copiedB, pc->B.getView());
        Kokkos::deep_copy(restoredR, pc->R.getView());

        double eScale = 0.0;
        double bScale = 0.0;
        for (size_t i = 0; i < pc->getLocalNum(); ++i) {
            for (unsigned d = 0; d < 3; ++d) {
                eScale = std::max(eScale, std::abs(physicalE(i)[d]));
                bScale = std::max(bScale, std::abs(physicalB(i)[d]));
            }
        }
        EXPECT_GT(eScale, 0.0);
        EXPECT_GT(bScale, 0.0);

        const double eTolerance = std::max(1.0e-12, 1.0e-10 * eScale);
        const double bTolerance = std::max(1.0e-12, 1.0e-10 * bScale);
        for (size_t i = 0; i < pc->getLocalNum(); ++i) {
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_NEAR(copiedE(i)[d], 2.0 * physicalE(i)[d], eTolerance);
                EXPECT_NEAR(copiedB(i)[d], 0.0, bTolerance);
                EXPECT_DOUBLE_EQ(restoredR(i)[d], originalR(i)[d]);
            }
        }
    }

    TEST_F(BinnedFieldSolverSmokeTest, P3MOpenAndPeriodicUseSameSolverWrapperAndSelectedLayoutBC) {
        for (const auto& [boundary, expectedParticleBC] :
             {std::pair{"PERIODIC", ippl::BC::PERIODIC}, std::pair{"OPEN", ippl::BC::NO}}) {
            ASSERT_NO_THROW(rebuildP3MBunch(boundary));
            ASSERT_NE(bunch->getFieldSolver(), nullptr);
            EXPECT_EQ(bunch->getFieldSolver()->getStype(), "P3M");
            EXPECT_TRUE((std::holds_alternative<FFTTruncatedGreenSolver_t<double, 3>>(
                    bunch->getFieldSolver()->getSolver())));
            EXPECT_DOUBLE_EQ(bunch->getFieldSolver()->getP3MCutoff(), 0.25);
            ASSERT_TRUE(pc->hasP3MLayout());
            for (const auto bc : pc->getP3MLayout().getParticleBC()) {
                EXPECT_EQ(bc, expectedParticleBC);
            }

            createParticles(2, /*pzMin=*/0.1, /*pzMax=*/0.2);
            EXPECT_NO_THROW(bunch->computeSelfFields());
        }
    }

    TEST_F(BinnedFieldSolverSmokeTest, P3MRejectsMixedAndDirichletBoundaries) {
        fsCmd->setType("P3M");
        fsCmd->setP3MCutoff(0.25);
        fsCmd->setParallelDecomposition(true);

        fsCmd->setBCX("OPEN");
        fsCmd->setBCY("OPEN");
        fsCmd->setBCZ("PERIODIC");
        EXPECT_THROW(fsCmd->execute(), OpalException);

        fsCmd->setBCX("DIRICHLET");
        fsCmd->setBCY("DIRICHLET");
        fsCmd->setBCZ("DIRICHLET");
        EXPECT_THROW(fsCmd->execute(), OpalException);
    }

    TEST_F(BinnedFieldSolverSmokeTest, BunchUpdate_ImageChargeBoundsIncludeMirroredZ) {
        constexpr double zPlane = 0.0;
        bunch->setImageChargeConfiguration(true, zPlane);

        createParticles(kDefaultNParticles, /*pzMin=*/0.1, /*pzMax=*/0.9);
        bunch->bunchUpdate();

        pc->computeMinMaxR();
        const auto minR    = pc->getMinR();
        const auto maxR    = pc->getMaxR();
        const auto gridMin = bunch->rmin_m;
        const auto gridMax = bunch->rmax_m;

        const double mirroredMinZ = 2.0 * zPlane - maxR[2];
        const double mirroredMaxZ = 2.0 * zPlane - minR[2];

        EXPECT_LE(gridMin[2], std::min(minR[2], mirroredMinZ));
        EXPECT_GE(gridMax[2], std::max(maxR[2], mirroredMaxZ));
    }

    TEST_F(BinnedFieldSolverSmokeTest, MirrorFieldZHandlesNonSlab3DDecomposition) {
        if (ippl::Comm->size() != 4) {
            GTEST_SKIP() << "This mirror-field decomposition check is defined for 4 MPI ranks.";
        }

        constexpr unsigned Dim = 3;
        ippl::NDIndex<Dim> domain;
        domain[0] = ippl::Index(8);
        domain[1] = ippl::Index(8);
        domain[2] = ippl::Index(8);

        std::array<bool, Dim> decomp;
        decomp.fill(true);

        Vector_t<double, Dim> hx(1.0);
        Vector_t<double, Dim> origin(0.0);
        FieldLayout_t<Dim> layout(MPI_COMM_WORLD, domain, decomp, false);
        Mesh_t<Dim> mesh(domain, hx, origin);

        std::vector<ippl::NDIndex<Dim>> domains(4);
        domains[0][0] = ippl::Index(0, 3);
        domains[0][1] = ippl::Index(0, 7);
        domains[0][2] = ippl::Index(0, 3);

        domains[1][0] = ippl::Index(4, 7);
        domains[1][1] = ippl::Index(0, 7);
        domains[1][2] = ippl::Index(0, 3);

        domains[2][0] = ippl::Index(0, 7);
        domains[2][1] = ippl::Index(0, 3);
        domains[2][2] = ippl::Index(4, 7);

        domains[3][0] = ippl::Index(0, 7);
        domains[3][1] = ippl::Index(4, 7);
        domains[3][2] = ippl::Index(4, 7);

        layout.updateLayout(domains);

        VField_t<double, Dim> src;
        VField_t<double, Dim> dst;
        src.initialize(mesh, layout);
        dst.initialize(mesh, layout);

        const auto& localDomain = layout.getLocalNDIndex();
        const int nghost        = src.getNghost();
        auto srcHost            = Kokkos::create_mirror_view(src.getView());
        for (size_t i = 0; i < srcHost.extent(0); ++i) {
            for (size_t j = 0; j < srcHost.extent(1); ++j) {
                for (size_t k = 0; k < srcHost.extent(2); ++k) {
                    srcHost(i, j, k) = Vector_t<double, Dim>(0.0);
                }
            }
        }

        for (int i = localDomain[0].first(); i <= localDomain[0].last(); ++i) {
            for (int j = localDomain[1].first(); j <= localDomain[1].last(); ++j) {
                for (int k = localDomain[2].first(); k <= localDomain[2].last(); ++k) {
                    const int localI  = i - localDomain[0].first() + nghost;
                    const int localJ  = j - localDomain[1].first() + nghost;
                    const int localK  = k - localDomain[2].first() + nghost;
                    const double base = static_cast<double>(i + 10 * j + 100 * k);

                    Vector_t<double, Dim> value;
                    value[0]                        = base;
                    value[1]                        = 1000.0 + base;
                    value[2]                        = 2000.0 + base;
                    srcHost(localI, localJ, localK) = value;
                }
            }
        }
        Kokkos::deep_copy(src.getView(), srcHost);

        opalx::detail::mirrorField(src, dst, Dim - 1);

        auto dstHost = Kokkos::create_mirror_view(dst.getView());
        Kokkos::deep_copy(dstHost, dst.getView());

        const int globalLastZ = domain[2].last();
        for (int i = localDomain[0].first(); i <= localDomain[0].last(); ++i) {
            for (int j = localDomain[1].first(); j <= localDomain[1].last(); ++j) {
                for (int k = localDomain[2].first(); k <= localDomain[2].last(); ++k) {
                    const int localI    = i - localDomain[0].first() + nghost;
                    const int localJ    = j - localDomain[1].first() + nghost;
                    const int localK    = k - localDomain[2].first() + nghost;
                    const int mirroredK = globalLastZ - k;
                    const double base   = static_cast<double>(i + 10 * j + 100 * mirroredK);

                    EXPECT_DOUBLE_EQ(dstHost(localI, localJ, localK)[0], base);
                    EXPECT_DOUBLE_EQ(dstHost(localI, localJ, localK)[1], 1000.0 + base);
                    EXPECT_DOUBLE_EQ(dstHost(localI, localJ, localK)[2], 2000.0 + base);
                }
            }
        }
    }

}  // namespace
