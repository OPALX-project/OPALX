/**
 * @file TestBinnedFieldSolver.cpp
 * @brief Smoke tests for `PartBunch::computeSelfFields()` with and without particle binning
 * enabled.
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
 *
 * Notes:
 * - The fixture initializes IPPL and disables HDF5 output to keep the tests lightweight.
 * - The goal is a robustness/smoke check (no physics validation of non-trivial fields, since this
 * would require way more computational resources).
 */

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <random>
#include <string>

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "PartBunch/PartBunch.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utility/Inform.h"
#include "Utilities/Options.h"

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
};

using PartBunch_t = PartBunch<double, 3>;
using ParticleContainer_t = typename PartBunch_t::ParticleContainer_t;
using AdaptBins_t = typename PartBunch_t::AdaptBins_t;
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

        // Keep the concrete solver command alive; PartBunch borrows it.
        fsCmdBase = fsCmd;

        dataSink = std::make_shared<DataSink>();
        beam = std::make_shared<Beam>();
        Beam* testBeam = Beam::find("UNNAMED_BEAM");

        // qi/mi/lbt are used by rho scaling; but with NullSolver we mostly validate
        // "no-throw" and deterministic zero E behavior.
        bunch = std::make_shared<PartBunch_t>(
            /*qi=*/std::vector<double>{1.0},
            /*mi=*/std::vector<double>{1.0},
            /*beams=*/std::vector<Beam*>{testBeam},
            /*totalParticlesPerBeam=*/std::vector<size_t>{kDefaultNParticles},
            /*lbt=*/1.0,
            /*integration_method=*/"LF2",
            fsCmdBase.get(),
            dataSink.get());
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
        pc->create(nPart);

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

    std::shared_ptr<AdaptBins_t> attachBins(typename AdaptBins_t::bin_index_type maxBins,
                                             double alpha,
                                             double beta,
                                             double desiredWidth) {
        using ConcreteBins_t =
            ParticleBinning::AdaptBins<ParticleContainer_t, CoordinateSelector_t>;

        CoordinateSelector_t selector(/*axis=*/2);
        auto bins = std::make_shared<ConcreteBins_t>(
            *pc,
            selector,
            maxBins,
            alpha,
            beta,
            desiredWidth,
            /*binningCmdName=*/"TEST_BINNING_CMD");
        bunch->setBins(bins);
        return bins;
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

}  // namespace
