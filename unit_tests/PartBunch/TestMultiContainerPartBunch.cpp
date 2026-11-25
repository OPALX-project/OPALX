/**
 * @file TestMultiContainerPartBunch.cpp
 * @brief Unit tests for multi-container PartBunch construction, indexed access,
 *        aggregation, and DataSink diagnostic contracts (stems, fdext sizing).
 *
 * The construction tests also check the P3M overlap layouts: a mesh-tied cutoff
 * (`RCUT_CELLS`) that does not fit the rank regions is reduced when the bunch is built, as the
 * runtime domain update does, and a cutoff that fits is kept. The space-charge factory test
 * checks that the configured Cartesian algorithm is built.
 */

#include <gtest/gtest.h>

#include "../BunchBeam.h"
#include "../TemporaryDirectory.h"

#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "AbstractObjects/OpalData.h"
#include "Ippl.h"
#include "PartBunch/PartBunch.h"
#include "Processes/GlobalProcesses/GlobalProcess.h"
#include "SpaceCharge/SpaceChargeFactory.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/H5PartWrapperForPT.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

namespace {

    using PartBunch_t         = PartBunch<double, 3>;
    using ParticleContainer_t = typename PartBunch_t::ParticleContainer_t;

    constexpr size_t kParticlesPerBeam = 32;

    class MultiContainerPartBunchTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);

            gmsg                = new Inform(nullptr, -1);
            Options::enableHDF5 = false;
        }

        static void TearDownTestSuite() {
            delete gmsg;
            gmsg = nullptr;
            ippl::finalize();
        }

        void SetUp() override {
            auto* opal           = OpalData::getInstance();
            savedInput_m         = opal->getInputFn();
            const auto directory = outputDirectory_m.create(
                    "opalx-multi-container", ippl::Comm->getCommunicator());
            opal->storeInputFn((directory / "unit_test.opal").string());
            storageConfig.periodicParticleBoundary = true;

            dataSink = std::make_shared<DataSink>(
                    DataSinkOptions{}, std::vector<H5PartWrapper*>{}, 1, outputStem(),
                    "unit_test.opal");
            beam = std::make_shared<Beam>();
            Beam::defineBuiltins();
            testBeam = Beam::find("UNNAMED_BEAM");
            ASSERT_NE(testBeam, nullptr);

            const double q0 = 1.0e-9;
            const double q1 = 2.0e-9;
            const double m0 = 0.511e-3;
            const double m1 = 0.938;

            bunch = std::make_shared<PartBunch_t>(
                    std::vector{
                            opalx::test::bunchBeam(*testBeam, kParticlesPerBeam, q0, m0),
                            opalx::test::bunchBeam(*testBeam, kParticlesPerBeam, q1, m1)},
                    /*lbt=*/1.0, storageConfig);
            q0_m = q0;
            q1_m = q1;
            m0_m = m0;
            m1_m = m1;
        }

        void TearDown() override {
            bunch.reset();
            dataSink.reset();
            beam.reset();
            // Close every writer before deleting only this fixture's diagnostic files.
            ippl::Comm->barrier();
            OpalData::getInstance()->storeInputFn(savedInput_m);
            outputDirectory_m.reset();
            ippl::Comm->barrier();
        }

        void createParticlesInContainer(
                size_t containerIndex, size_t nPart, double pzMin, double pzMax) {
            auto pc           = bunch->getParticleContainer(containerIndex);
            const int mpiRank = ippl::Comm->rank();
            std::mt19937_64 eng(1000u + containerIndex * 97u + static_cast<unsigned>(mpiRank));

            pc->createParticles(nPart);

            auto R_host  = pc->R.getHostMirror();
            auto P_host  = pc->P.getHostMirror();
            auto dt_host = pc->dt.getHostMirror();
            auto E_host  = pc->E.getHostMirror();

            const auto rmin = bunch->cartesianDomain().lower();
            const auto rmax = bunch->cartesianDomain().upper();

            std::uniform_real_distribution<double> unifR_x(rmin[0] + 0.1, rmax[0] - 0.1);
            std::uniform_real_distribution<double> unifR_y(rmin[1] + 0.1, rmax[1] - 0.1);
            std::uniform_real_distribution<double> unifR_z(rmin[2] + 0.1, rmax[2] - 0.1);
            std::uniform_real_distribution<double> unifP_z(pzMin, pzMax);

            const double dt = bunch->getdT();
            const double qi = (containerIndex == 0) ? q0_m : q1_m;

            for (size_t i = 0; i < nPart; ++i) {
                R_host(i)[0] = unifR_x(eng);
                R_host(i)[1] = unifR_y(eng);
                R_host(i)[2] = unifR_z(eng);

                P_host(i)[0] = 0.0;
                P_host(i)[1] = 0.0;
                P_host(i)[2] = unifP_z(eng);

                dt_host(i) = dt;

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

        /// Snapshots for every container holding what the statistics and H5 outputs read.
        std::vector<opalx::moments::BeamMoments> outputMoments() const {
            std::vector<opalx::moments::BeamMoments> moments;
            for (const auto& pc : bunch->getParticleContainers()) {
                moments.push_back(
                        opalx::moments::compute(
                                *pc, dataSink->requiredMoments({true, true}, pc->hasSpin())));
            }
            return moments;
        }

        /// Stem of the diagnostic files the sinks of a test write, in the test directory.
        std::string outputStem() const { return (outputDirectory_m.path() / "unit_test").string(); }

        static std::array<Vector_t<double, 3>, 2> zeroFdPair() {
            std::array<Vector_t<double, 3>, 2> fd{};
            fd[0] = Vector_t<double, 3>(0.0);
            fd[1] = Vector_t<double, 3>(0.0);
            return fd;
        }

        opalx::test::TemporaryDirectory outputDirectory_m;
        std::string savedInput_m;
        opalx::spacecharge::CartesianDomainConfig3D storageConfig;
        std::shared_ptr<DataSink> dataSink;
        std::shared_ptr<Beam> beam;
        Beam* testBeam = nullptr;
        std::shared_ptr<PartBunch_t> bunch;
        double q0_m = 0.0;
        double q1_m = 0.0;
        double m0_m = 0.0;
        double m1_m = 0.0;
    };

    // --- PartBunch / container access ---

    TEST_F(MultiContainerPartBunchTest, TwoContainers_CountsAndDistinctPointers) {
        ASSERT_EQ(bunch->getNumParticleContainers(), 2u);
        ASSERT_EQ(bunch->getParticleContainers().size(), 2u);

        auto pc0 = bunch->getParticleContainer(0);
        auto pc1 = bunch->getParticleContainer(1);
        ASSERT_NE(pc0, nullptr);
        ASSERT_NE(pc1, nullptr);
        EXPECT_NE(pc0.get(), pc1.get());

        EXPECT_EQ(bunch->getParticleContainer(0).get(), bunch->getParticleContainer().get());
    }

    TEST_F(MultiContainerPartBunchTest, TwoContainers_PerContainerChargeAndMass) {
        auto pc0 = bunch->getParticleContainer(0);
        auto pc1 = bunch->getParticleContainer(1);

        EXPECT_DOUBLE_EQ(pc0->getChargePerParticle(), q0_m);
        EXPECT_DOUBLE_EQ(pc1->getChargePerParticle(), q1_m);
        EXPECT_DOUBLE_EQ(pc0->getMassPerParticle(), m0_m);
        EXPECT_DOUBLE_EQ(pc1->getMassPerParticle(), m1_m);
    }

    TEST_F(MultiContainerPartBunchTest, GetParticleContainer_IndexOutOfRange) {
        EXPECT_THROW(static_cast<void>(bunch->getParticleContainer(2)), std::out_of_range);
    }

    TEST_F(MultiContainerPartBunchTest, ActiveFlag_IndexOutOfRangeThrows) {
        EXPECT_THROW(static_cast<void>(bunch->isPcActive(2)), std::out_of_range);
        EXPECT_THROW(static_cast<void>(bunch->pcAtSStop(2)), std::out_of_range);
        EXPECT_THROW(static_cast<void>(bunch->getPcRoot(2)), std::out_of_range);
    }

    /// A global process that only counts how often it is applied.
    class CountingProcess : public GlobalProcess {
    public:
        explicit CountingProcess(int& calls) : calls_m(calls) {}
        size_t apply(ParticleContainer<double, 3>&, double, long long, size_t) override {
            ++calls_m;
            return 0;
        }

    private:
        int& calls_m;
    };

    TEST_F(MultiContainerPartBunchTest, ApplyGlobalProcesses_SkipsStoppedContainer) {
        PartBunch_t local(
                std::vector(2, opalx::test::bunchBeam(*testBeam, kParticlesPerBeam)), 1.0,
                storageConfig);
        int calls = 0;
        std::vector<std::unique_ptr<GlobalProcess>> processes;
        processes.push_back(std::make_unique<CountingProcess>(calls));
        local.getParticleContainer(0)->setGlobalProcesses(std::move(processes));
        local.getParticleContainer(0)->createParticles(4);
        local.resetPcActive();
        local.setdT(1e-12);
        ASSERT_TRUE(local.isPcActive(0));

        local.applyGlobalProcesses();
        EXPECT_EQ(calls, 1);

        // A beam that has passed its ZSTOP does not decay any more.
        local.setPcAtSStop(0);
        local.applyGlobalProcesses();
        EXPECT_EQ(calls, 1);
    }

    TEST_F(MultiContainerPartBunchTest, Constructor_TimeStepStartsAtZero) {
        EXPECT_EQ(bunch->getdT(), 0.0);
    }

    TEST_F(MultiContainerPartBunchTest, Constructor_ThrowsWithoutBeams) {
        EXPECT_THROW(
                static_cast<void>(std::make_shared<PartBunch_t>(
                        std::vector<BunchBeam>{}, 1.0, storageConfig)),
                OpalException);
    }

    // A mesh-tied cutoff follows the 6 m placeholder domain set at construction, before the
    // runtime geometry exists, so it must be reduced to fit the rank regions as at runtime.
    TEST_F(MultiContainerPartBunchTest, Constructor_ReducesMeshTiedCutoffToFitRankRegions) {
        opalx::spacecharge::CartesianDomainConfig3D overlapConfig;
        overlapConfig.layoutType         = opalx::spacecharge::ParticleLayoutType::SpatialOverlap;
        overlapConfig.overlapCutoffCells = 5.0;
        // Five 0.75 m cells exceed half of the 6 m domain on every rank count.
        const double requested = 5.0 * 6.0 / 8.0;

        std::shared_ptr<PartBunch_t> overlapBunch;
        ASSERT_NO_THROW(
                overlapBunch = std::make_shared<PartBunch_t>(
                        std::vector{
                                opalx::test::bunchBeam(*testBeam, kParticlesPerBeam, q0_m, m0_m),
                                opalx::test::bunchBeam(*testBeam, kParticlesPerBeam, q1_m, m1_m)},
                        1.0, overlapConfig));

        const auto& domain      = overlapBunch->cartesianDomain();
        const auto& localDomain = domain.layout().getLocalNDIndex();
        const auto& spacing     = domain.mesh().getMeshSpacing();
        const double firstCutoff =
                overlapBunch->getParticleContainer(0)->getP3MLayout().getCutoff();
        EXPECT_LT(firstCutoff, requested);
        for (unsigned dimension = 0; dimension < 3; ++dimension) {
            EXPECT_LE(firstCutoff, 0.5 * localDomain[dimension].length() * spacing[dimension]);
        }
        if (ippl::Comm->size() == 1) {
            EXPECT_NEAR(firstCutoff, 3.0, 1.0e-9);
        }
        EXPECT_DOUBLE_EQ(
                overlapBunch->getParticleContainer(1)->getP3MLayout().getCutoff(), firstCutoff);
    }

    TEST_F(MultiContainerPartBunchTest, Constructor_KeepsMeshTiedCutoffThatFits) {
        opalx::spacecharge::CartesianDomainConfig3D overlapConfig;
        overlapConfig.layoutType         = opalx::spacecharge::ParticleLayoutType::SpatialOverlap;
        overlapConfig.overlapCutoffCells = 1.0;

        const auto overlapBunch = std::make_shared<PartBunch_t>(
                std::vector{opalx::test::bunchBeam(*testBeam, kParticlesPerBeam, q0_m, m0_m)}, 1.0,
                overlapConfig);
        EXPECT_DOUBLE_EQ(
                overlapBunch->getParticleContainer()->getP3MLayout().getCutoff(), 6.0 / 8.0);
    }

    TEST_F(MultiContainerPartBunchTest, TotalNumAllContainers_SumsBoth) {
        const size_t n0 = 11u;
        const size_t n1 = 7u;
        createParticlesInContainer(0, n0, 0.1, 0.5);
        createParticlesInContainer(1, n1, 0.2, 0.6);

        EXPECT_EQ(bunch->getTotalNumAllContainers(), n0 + n1);
    }

    TEST_F(MultiContainerPartBunchTest, SpaceChargeFactoryConstructsConfiguredCartesianAlgorithm) {
        using namespace opalx::spacecharge;
        CartesianPIC3DConfig values;
        values.backend                         = PoissonSolverType::PeriodicFFT;
        values.grid.meshSize                   = storageConfig.meshSize;
        values.grid.decomposition              = storageConfig.decomposition;
        values.grid.boundingBoxIncreasePercent = storageConfig.boundingBoxIncreasePercent;
        values.boundaryConditions              = {
                FieldBoundaryCondition::Periodic, FieldBoundaryCondition::Periodic,
                FieldBoundaryCondition::Periodic};
        SpaceChargeConfig config = values;

        auto solver = makeSpaceChargeSolver(std::move(config), *bunch, dataSink.get());
        EXPECT_NE(solver, nullptr);
    }

    // --- DataSink stems and writers ---

    TEST_F(MultiContainerPartBunchTest, DiagnosticStemForContainer_SingleVsMulti) {
        EXPECT_EQ(DataSink::diagnosticStemForContainer("run", 1, 0), "run");
        EXPECT_EQ(DataSink::diagnosticStemForContainer("run", 1, 99), "run");

        EXPECT_EQ(DataSink::diagnosticStemForContainer("run", 2, 0), "run_c0");
        EXPECT_EQ(DataSink::diagnosticStemForContainer("run", 2, 1), "run_c1");
    }

    TEST_F(MultiContainerPartBunchTest, DataSink_MultiContainerInit_DoesNotThrow) {
        EXPECT_NO_THROW(DataSink(DataSinkOptions{}, {}, 2, outputStem(), "unit_test.opal"));
    }

    TEST_F(MultiContainerPartBunchTest, DataSink_dumpSDDS_ThrowsWhenFdextTooSmall) {
        createParticlesInContainer(0, 8u, 0.1, 0.3);
        createParticlesInContainer(1, 6u, 0.1, 0.3);

        DataSink ds(DataSinkOptions{}, {}, 2, outputStem(), "unit_test.opal");
        std::vector<std::array<Vector_t<double, 3>, 2>> fdTooSmall;
        fdTooSmall.push_back(zeroFdPair());

        EXPECT_THROW(ds.dumpSDDS(*bunch, fdTooSmall, outputMoments(), 1), OpalException);
    }

    TEST_F(MultiContainerPartBunchTest, DataSink_dumpSDDS_NoThrowWithMatchingFdextAndHdf5Off) {
        createParticlesInContainer(0, 8u, 0.1, 0.3);
        createParticlesInContainer(1, 6u, 0.1, 0.3);

        DataSink ds(DataSinkOptions{}, {}, 2, outputStem(), "unit_test.opal");
        std::vector<std::array<Vector_t<double, 3>, 2>> fd(2);
        fd[0] = zeroFdPair();
        fd[1] = zeroFdPair();

        EXPECT_NO_THROW(ds.dumpSDDS(*bunch, fd, outputMoments(), 1));
    }

    // When HDF5 is disabled, dumpH5 returns before validating fdext size.
    TEST_F(MultiContainerPartBunchTest, DataSink_dumpH5_NoThrowWhenHdf5OffEvenIfFdextTooSmall) {
        createParticlesInContainer(0, 4u, 0.1, 0.2);

        DataSink ds(DataSinkOptions{}, {}, 2, outputStem(), "unit_test.opal");
        std::vector<std::array<Vector_t<double, 3>, 2>> fdTooSmall;
        fdTooSmall.push_back(zeroFdPair());

        EXPECT_NO_THROW(ds.dumpH5(*bunch, fdTooSmall, outputMoments()));
        EXPECT_TRUE(ds.requiredMoments({true, false}, false).empty());
    }

    TEST_F(MultiContainerPartBunchTest, DataSink_dumpH5_ThrowsWhenFdextTooSmallAndHdf5On) {
        createParticlesInContainer(0, 4u, 0.1, 0.2);

        DataSinkOptions options;
        options.enableHDF5 = true;

        const std::string f0 = (outputDirectory_m.path() / "test_mc_h5_c0.h5").string();
        const std::string f1 = (outputDirectory_m.path() / "test_mc_h5_c1.h5").string();

        std::vector<H5PartWrapper*> wrappers;
        wrappers.push_back(new H5PartWrapperForPT(f0, options, H5_O_WRONLY));
        wrappers.push_back(new H5PartWrapperForPT(f1, options, H5_O_WRONLY));

        {
            DataSink ds(options, wrappers, 2, outputStem(), "unit_test.opal");
            std::vector<std::array<Vector_t<double, 3>, 2>> fdTooSmall;
            fdTooSmall.push_back(zeroFdPair());
            EXPECT_THROW(ds.dumpH5(*bunch, fdTooSmall, outputMoments()), OpalException);
            EXPECT_EQ(
                    ds.requiredMoments({true, false}, false), H5PartWrapperForPT::requiredMoments);
        }

        for (H5PartWrapper* w : wrappers) {
            delete w;
        }
    }

    TEST_F(MultiContainerPartBunchTest, DataSink_RewindToCheckpointAppendsDiagnosticH5InPlace) {
        createParticlesInContainer(0, 4u, 0.1, 0.2);
        createParticlesInContainer(1, 4u, 0.1, 0.2);

        DataSinkOptions options;
        options.enableHDF5 = true;

        const std::string base =
                (outputDirectory_m.path() / "test_diagnostic_checkpoint_rewind").string();
        const std::string f0 = base + "_c0.h5";
        const std::string f1 = base + "_c1.h5";
        auto w0              = std::make_unique<H5PartWrapperForPT>(f0, options, H5_O_WRONLY);
        auto w1              = std::make_unique<H5PartWrapperForPT>(f1, options, H5_O_WRONLY);
        std::vector<H5PartWrapper*> wrappers{w0.get(), w1.get()};

        {
            DataSink ds(options, wrappers, 2, base, "unit_test.opal");
            std::vector<std::array<Vector_t<double, 3>, 2>> fd(2);
            fd[0] = zeroFdPair();
            fd[1] = zeroFdPair();

            bunch->setdT(1.0e-12);
            bunch->setT(1.0e-12);
            ds.dumpH5(*bunch, fd, outputMoments());
            bunch->setT(2.0e-12);
            ds.dumpH5(*bunch, fd, outputMoments());

            bunch->setT(1.0e-12);
            ds.rewindToCheckpoint(*bunch);
            bunch->setT(2.0e-12);
            ds.dumpH5(*bunch, fd, outputMoments());
        }
        w0.reset();
        w1.reset();

        h5_prop_t props = H5CreateFileProp();
        MPI_Comm comm   = ippl::Comm->getCommunicator();
        ASSERT_NE(props, H5_ERR);
        ASSERT_NE(H5SetPropFileMPIOCollective(props, &comm), H5_ERR);
        h5_file_t file = H5OpenFile(f0.c_str(), H5_O_RDONLY, props);
        H5CloseProp(props);
        ASSERT_NE(file, static_cast<h5_file_t>(H5_ERR));
        EXPECT_EQ(H5GetNumSteps(file), 2);
        h5_float64_t time = 0.0;
        ASSERT_EQ(H5SetStep(file, 0), H5_SUCCESS);
        ASSERT_EQ(H5ReadStepAttribFloat64(file, "TIME", &time), H5_SUCCESS);
        EXPECT_DOUBLE_EQ(time, 1.0e-12);
        ASSERT_EQ(H5SetStep(file, 1), H5_SUCCESS);
        ASSERT_EQ(H5ReadStepAttribFloat64(file, "TIME", &time), H5_SUCCESS);
        EXPECT_DOUBLE_EQ(time, 2.0e-12);
        ASSERT_EQ(H5CloseFile(file), H5_SUCCESS);
    }

}  // namespace
