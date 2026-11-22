/**
 * @file TestMultiContainerPartBunch.cpp
 * @brief Unit tests for multi-container PartBunch construction, indexed access,
 *        aggregation, and DataSink diagnostic contracts (stems, fdext sizing).
 *
 * The construction tests also check the P3M overlap layouts: a mesh-tied cutoff
 * (`RCUT_CELLS`) that does not fit the rank regions is reduced when the bunch is built, as the
 * runtime domain update does, and a cutoff that fits is kept. The space-charge factory tests
 * check that the configured Cartesian algorithm is built and that FFT2D5 defers its setup.
 */

#include <gtest/gtest.h>

#include "../TemporaryDirectory.h"

#include <array>
#include <iterator>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "AbstractObjects/OpalData.h"
#include "Ippl.h"
#include "PartBunch/PartBunch.h"
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

            dataSink = std::make_shared<DataSink>();
            beam     = std::make_shared<Beam>();
            testBeam = Beam::find("UNNAMED_BEAM");
            ASSERT_NE(testBeam, nullptr);

            const double q0 = 1.0e-9;
            const double q1 = 2.0e-9;
            const double m0 = 0.511e-3;
            const double m1 = 0.938;

            bunch = std::make_shared<PartBunch_t>(
                    std::vector<double>{q0, q1}, std::vector<double>{m0, m1},
                    std::vector<Beam*>{testBeam, testBeam},
                    std::vector<size_t>{kParticlesPerBeam, kParticlesPerBeam},
                    /*lbt=*/1.0,
                    /*integration_method=*/"LF2", storageConfig);
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

            const auto rmin = bunch->rmin_m;
            const auto rmax = bunch->rmax_m;

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

    TEST_F(MultiContainerPartBunchTest, Constructor_ThrowsOnQiSizeMismatch) {
        EXPECT_THROW(
                static_cast<void>(std::make_shared<PartBunch_t>(
                        std::vector<double>{1.0}, std::vector<double>{1.0, 1.0},
                        std::vector<Beam*>{testBeam, testBeam},
                        std::vector<size_t>{kParticlesPerBeam, kParticlesPerBeam}, 1.0, "LF2",
                        storageConfig)),
                OpalException);
    }

    TEST_F(MultiContainerPartBunchTest, Constructor_ThrowsOnNullBeam) {
        EXPECT_THROW(
                static_cast<void>(std::make_shared<PartBunch_t>(
                        std::vector<double>{1.0, 1.0}, std::vector<double>{1.0, 1.0},
                        std::vector<Beam*>{testBeam, nullptr},
                        std::vector<size_t>{kParticlesPerBeam, kParticlesPerBeam}, 1.0, "LF2",
                        storageConfig)),
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
                        std::vector<double>{q0_m, q1_m}, std::vector<double>{m0_m, m1_m},
                        std::vector<Beam*>{testBeam, testBeam},
                        std::vector<size_t>{kParticlesPerBeam, kParticlesPerBeam}, 1.0, "LF2",
                        overlapConfig));

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
                std::vector<double>{q0_m}, std::vector<double>{m0_m}, std::vector<Beam*>{testBeam},
                std::vector<size_t>{kParticlesPerBeam}, 1.0, "LF2", overlapConfig);
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
        // Reconfiguration to a different solver must clear a previous FFT2D5 selection.
        for (const auto& pc : bunch->getParticleContainers())
            pc->setCompactSymmetricMoments(true);
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
        EXPECT_FALSE(solver->usesFFT2D5());
        for (const auto& pc : bunch->getParticleContainers())
            EXPECT_FALSE(pc->usesCompactSymmetricMoments());
    }

    TEST_F(MultiContainerPartBunchTest, SpaceChargeFactoryDefersFFT2D5Initialization) {
        using namespace opalx::spacecharge;
        FFT2D5Config values;
        values.referencePathFile = "";
        auto solver = makeSpaceChargeSolver(SpaceChargeConfig(values), *bunch, nullptr);

        ASSERT_NE(solver, nullptr);
        EXPECT_TRUE(solver->usesFFT2D5());
        for (const auto& pc : bunch->getParticleContainers())
            EXPECT_TRUE(pc->usesCompactSymmetricMoments());
        EXPECT_EQ(solver->backendSolveCount(), 0u);
        EXPECT_EQ(solver->redistributionCount(), 0u);
    }

    TEST_F(MultiContainerPartBunchTest, FFT2D5BeamParametersMatchLegacyWithEmptyAndUnevenRanks) {
        auto pc           = bunch->getParticleContainer();
        const size_t rank = ippl::Comm->rank();
        // Exercise an empty bunch, ownership on rank zero only, then uneven ownership.
        for (size_t added : {size_t(0), rank == 0 ? size_t(3) : size_t(0), rank + 1}) {
            pc->createParticles(added);
            auto r = pc->R.getHostMirror();
            auto p = pc->P.getHostMirror();
            for (size_t i = 0; i < pc->getLocalNum(); ++i) {
                r(i) = Vector_t<double, 3>(-0.3 + 0.01 * i, 0.2 + rank, -1.0 - i - rank);
                p(i) = Vector_t<double, 3>(0.001 * i, -0.002 * rank, 0.1 + 0.01 * i);
            }
            Kokkos::deep_copy(pc->R.getView(), r);
            Kokkos::deep_copy(pc->P.getView(), p);
            pc->markMomentsDirty();
            bunch->calcBeamParameters();
            const auto minimum    = bunch->rmin_m;
            const auto maximum    = bunch->rmax_m;
            const auto means      = pc->getMeans();
            const auto covariance = pc->getCovMatrix();
            const auto rmsR       = pc->getRmsR();
            const auto rmsP       = pc->getRmsP();

            pc->markMomentsDirty();
            bunch->rmin_m = 123.0;
            bunch->rmax_m = -123.0;
            bunch->calcBeamParameters(true);
            EXPECT_FALSE(pc->isMomentsDirty());
            for (size_t d = 0; d < 3; ++d) {
                EXPECT_EQ(bunch->rmin_m[d], minimum[d]);
                EXPECT_EQ(bunch->rmax_m[d], maximum[d]);
                EXPECT_EQ(pc->getRmsR()[d], rmsR[d]);
                EXPECT_EQ(pc->getRmsP()[d], rmsP[d]);
            }
            for (size_t i = 0; i < 6; ++i) {
                EXPECT_EQ(pc->getMeans()[i], means[i]);
                for (size_t j = 0; j < 6; ++j) {
                    EXPECT_EQ(pc->getCovMatrix()(i, j), covariance(i, j));
                }
            }
        }
    }

    // --- DataSink stems and writers ---

    TEST_F(MultiContainerPartBunchTest, DiagnosticStemForContainer_SingleVsMulti) {
        EXPECT_EQ(DataSink::diagnosticStemForContainer("run", 1, 0), "run");
        EXPECT_EQ(DataSink::diagnosticStemForContainer("run", 1, 99), "run");

        EXPECT_EQ(DataSink::diagnosticStemForContainer("run", 2, 0), "run_c0");
        EXPECT_EQ(DataSink::diagnosticStemForContainer("run", 2, 1), "run_c1");
    }

    TEST_F(MultiContainerPartBunchTest, DataSink_MultiContainerInit_DoesNotThrow) {
        EXPECT_NO_THROW(DataSink(std::vector<H5PartWrapper*>{}, false, 2));
    }

    TEST_F(MultiContainerPartBunchTest, DataSink_dumpSDDS_ThrowsWhenFdextTooSmall) {
        createParticlesInContainer(0, 8u, 0.1, 0.3);
        createParticlesInContainer(1, 6u, 0.1, 0.3);

        DataSink ds(std::vector<H5PartWrapper*>{}, false, 2);
        std::vector<std::array<Vector_t<double, 3>, 2>> fdTooSmall;
        fdTooSmall.push_back(zeroFdPair());

        const double azimuth = 0.0;
        EXPECT_THROW(ds.dumpSDDS(*bunch, fdTooSmall, 1, azimuth), OpalException);
    }

    TEST_F(MultiContainerPartBunchTest, DataSink_dumpSDDS_NoThrowWithMatchingFdextAndHdf5Off) {
        createParticlesInContainer(0, 8u, 0.1, 0.3);
        createParticlesInContainer(1, 6u, 0.1, 0.3);

        DataSink ds(std::vector<H5PartWrapper*>{}, false, 2);
        std::vector<std::array<Vector_t<double, 3>, 2>> fd(2);
        fd[0] = zeroFdPair();
        fd[1] = zeroFdPair();

        const double azimuth = 0.0;
        EXPECT_NO_THROW(ds.dumpSDDS(*bunch, fd, 1, azimuth));
    }

    // When HDF5 is disabled, dumpH5 returns before validating fdext size.
    TEST_F(MultiContainerPartBunchTest, DataSink_dumpH5_NoThrowWhenHdf5OffEvenIfFdextTooSmall) {
        ASSERT_FALSE(Options::enableHDF5);

        createParticlesInContainer(0, 4u, 0.1, 0.2);

        DataSink ds(std::vector<H5PartWrapper*>{}, false, 2);
        std::vector<std::array<Vector_t<double, 3>, 2>> fdTooSmall;
        fdTooSmall.push_back(zeroFdPair());

        EXPECT_NO_THROW(ds.dumpH5(*bunch, fdTooSmall));
    }

    TEST_F(MultiContainerPartBunchTest, DataSink_dumpH5_ThrowsWhenFdextTooSmallAndHdf5On) {
        createParticlesInContainer(0, 4u, 0.1, 0.2);

        const bool savedH5  = Options::enableHDF5;
        Options::enableHDF5 = true;

        const std::string f0 = (outputDirectory_m.path() / "test_mc_h5_c0.h5").string();
        const std::string f1 = (outputDirectory_m.path() / "test_mc_h5_c1.h5").string();

        std::vector<H5PartWrapper*> wrappers;
        wrappers.push_back(new H5PartWrapperForPT(f0, H5_O_WRONLY));
        wrappers.push_back(new H5PartWrapperForPT(f1, H5_O_WRONLY));

        {
            DataSink ds(wrappers, false, 2);
            std::vector<std::array<Vector_t<double, 3>, 2>> fdTooSmall;
            fdTooSmall.push_back(zeroFdPair());
            EXPECT_THROW(ds.dumpH5(*bunch, fdTooSmall), OpalException);
        }

        for (H5PartWrapper* w : wrappers) {
            delete w;
        }

        Options::enableHDF5 = savedH5;
    }

    // Regression for #612: deleting witnesses must move ID and phase space together,
    // including in each container's H5 output. No tracking kick is applied here, so
    // every surviving coordinate and normalized momentum must remain exactly unchanged.
    TEST_F(MultiContainerPartBunchTest, WitnessDeletionPreservesIdentityAndPhaseSpaceInH5) {
        constexpr int nLocal         = 8;
        constexpr size_t nContainers = 3;
        const int rank               = ippl::Comm->rank();
        const int ranks              = ippl::Comm->size();
        const auto comm              = ippl::Comm->getCommunicator();
        bunch                        = std::make_shared<PartBunch_t>(
                std::vector<double>{-1.0e-9, -1.0e-15, 1.0e-15},
                std::vector<double>(nContainers, 0.511e-3),
                std::vector<Beam*>(nContainers, testBeam), std::vector<size_t>(nContainers, nLocal),
                1.0, "LF2", storageConfig);
        bunch->setdT(1.0e-12);

        using PhaseSpace = std::array<double, 6>;
        using Snapshot   = std::map<long long, PhaseSpace>;
        std::array<std::array<Snapshot, nContainers>, 3> expected;
        // Remove interior and tail entries. On nonzero ranks the first round
        // empties the electron container, exercising collective I/O with no local data.
        const auto removed = [](size_t c, int owner, int i, int round) {
            if (c == 0 || round == 0) return false;
            const bool first = c == 1 ? (owner != 0 || i == 1 || i == 6 || i == 7)
                                      : (i == 0 || i == 3 || i == 7);
            return first || (round == 2 && i == 4);
        };
        const auto phase = [](size_t c, int owner, int i) {
            PhaseSpace result{};
            for (size_t d = 0; d < result.size(); ++d) {
                // Binary fractions give distinct, exactly representable sentinels:
                // positions in metres and normalized momenta in units of mc.
                result[d] = (1 + 100 * c + 10 * owner + i) / 1024.0 + d / 16.0;
            }
            return result;
        };
        for (size_t c = 0; c < nContainers; ++c) {
            auto pc = bunch->getParticleContainer(c);
            pc->createParticles(nLocal);
            auto r = pc->R.getHostMirror();
            auto p = pc->P.getHostMirror();
            for (int i = 0; i < nLocal; ++i) {
                const auto values = phase(c, rank, i);
                for (int d = 0; d < 3; ++d) {
                    r(i)[d] = values[d];
                    p(i)[d] = values[d + 3];
                }
            }
            Kokkos::deep_copy(pc->R.getView(), r);
            Kokkos::deep_copy(pc->P.getView(), p);
            pc->markMomentsDirty();
            const auto ids =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->ID.getView());
            std::vector<long long> localIds(nLocal), allIds(nLocal * ranks);
            for (int i = 0; i < nLocal; ++i)
                localIds[i] = ids(i);
            MPI_Allgather(
                    localIds.data(), nLocal, MPI_LONG_LONG, allIds.data(), nLocal, MPI_LONG_LONG,
                    comm);
            for (int round = 0; round < 3; ++round) {
                for (int owner = 0; owner < ranks; ++owner) {
                    for (int i = 0; i < nLocal; ++i) {
                        if (!removed(c, owner, i, round)) {
                            EXPECT_TRUE(
                                    expected[round][c]
                                            .emplace(allIds[owner * nLocal + i], phase(c, owner, i))
                                            .second);
                        }
                    }
                }
            }
        }

        const bool savedH5  = Options::enableHDF5;
        Options::enableHDF5 = true;
        struct RestoreH5Option {
            bool value;
            ~RestoreH5Option() { Options::enableHDF5 = value; }
        } restoreH5{savedH5};
        std::vector<std::unique_ptr<H5PartWrapperForPT>> owners;
        std::vector<H5PartWrapper*> wrappers;
        std::array<std::string, nContainers> files;
        for (size_t c = 0; c < nContainers; ++c) {
            files[c] = (outputDirectory_m.path() / ("deletion_c" + std::to_string(c) + ".h5"))
                               .string();
            owners.push_back(std::make_unique<H5PartWrapperForPT>(files[c], H5_O_WRONLY));
            wrappers.push_back(owners.back().get());
        }
        {
            DataSink ds(wrappers, false, nContainers);
            std::vector<std::array<Vector_t<double, 3>, 2>> fd(nContainers, zeroFdPair());
            for (int round = 0; round < 3; ++round) {
                for (size_t c = 0; c < nContainers; ++c) {
                    SCOPED_TRACE(::testing::Message() << "round=" << round << " container=" << c);
                    auto pc = bunch->getParticleContainer(c);
                    if (round > 0 && c > 0) {
                        const auto r = Kokkos::create_mirror_view_and_copy(
                                Kokkos::HostSpace(), pc->R.getView());
                        auto mask = pc->InvalidMask.getHostMirror();
                        for (size_t j = 0; j < pc->getLocalNum(); ++j) {
                            // Decode original local slot independently of ID and current index.
                            const int i =
                                    static_cast<int>(r(j)[0] * 1024) - 1 - 100 * c - 10 * rank;
                            mask(j) = removed(c, rank, i, round);
                        }
                        Kokkos::deep_copy(pc->InvalidMask.getView(), mask);
                        const auto deleted = pc->deleteInvalidParticles();
                        EXPECT_EQ(
                                deleted, expected[round - 1][c].size() - expected[round][c].size());
                    }
                    const auto ids = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), pc->ID.getView());
                    const auto r = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), pc->R.getView());
                    const auto p = Kokkos::create_mirror_view_and_copy(
                            Kokkos::HostSpace(), pc->P.getView());
                    std::vector<int> localSeen(expected[round][c].size(), 0),
                            globalSeen(localSeen.size());
                    for (size_t i = 0; i < pc->getLocalNum(); ++i) {
                        const auto found = expected[round][c].find(ids(i));
                        EXPECT_NE(found, expected[round][c].end());
                        if (found == expected[round][c].end()) continue;
                        ++localSeen[std::distance(expected[round][c].begin(), found)];
                        for (int d = 0; d < 3; ++d) {
                            EXPECT_EQ(r(i)[d], found->second[d]);
                            EXPECT_EQ(p(i)[d], found->second[d + 3]);
                        }
                    }
                    MPI_Allreduce(
                            localSeen.data(), globalSeen.data(), globalSeen.size(), MPI_INT,
                            MPI_SUM, comm);
                    for (int count : globalSeen)
                        EXPECT_EQ(count, 1);
                    EXPECT_EQ(pc->getTotalNum(), expected[round][c].size());
                }
                bunch->setT((round + 1) * bunch->getdT());
                ds.dumpH5(*bunch, fd);
            }
        }
        owners.clear();

        for (size_t c = 0; c < nContainers; ++c) {
            auto props      = H5CreateFileProp();
            MPI_Comm ioComm = comm;
            ASSERT_EQ(H5SetPropFileMPIOCollective(props, &ioComm), H5_SUCCESS);
            auto file = H5OpenFile(files[c].c_str(), H5_O_RDONLY, props);
            H5CloseProp(props);
            ASSERT_NE(file, static_cast<h5_file_t>(H5_ERR));
            EXPECT_EQ(H5GetNumSteps(file), 3);
            for (int round = 0; round < 3; ++round) {
                SCOPED_TRACE(::testing::Message() << "H5 round=" << round << " container=" << c);
                ASSERT_EQ(H5SetStep(file, round), H5_SUCCESS);
                // H5hut keeps the preceding step's read view; clear it before
                // querying the new global population (which decreases here).
                ASSERT_EQ(H5PartSetView(file, -1, -1), H5_SUCCESS);
                const auto n = H5PartGetNumParticles(file);
                ASSERT_EQ(n, expected[round][c].size());
                ASSERT_EQ(H5PartSetView(file, 0, n - 1), H5_SUCCESS);
                std::vector<h5_int64_t> ids(n);
                ASSERT_EQ(H5PartReadDataInt64(file, "id", ids.data()), H5_SUCCESS);
                std::array<std::vector<h5_float64_t>, 6> values;
                const std::array<const char*, 6> names{"x", "y", "z", "px", "py", "pz"};
                for (size_t d = 0; d < names.size(); ++d) {
                    values[d].resize(n);
                    ASSERT_EQ(H5PartReadDataFloat64(file, names[d], values[d].data()), H5_SUCCESS);
                }
                Snapshot actual;
                for (h5_ssize_t i = 0; i < n; ++i) {
                    PhaseSpace value{};
                    for (size_t d = 0; d < value.size(); ++d)
                        value[d] = values[d][i];
                    EXPECT_TRUE(actual.emplace(ids[i], value).second) << "duplicate H5 ID";
                }
                EXPECT_EQ(actual, expected[round][c]);
            }
            EXPECT_EQ(H5CloseFile(file), H5_SUCCESS);
        }
    }

    TEST_F(MultiContainerPartBunchTest, DataSink_RewindToCheckpointAppendsDiagnosticH5InPlace) {
        createParticlesInContainer(0, 4u, 0.1, 0.2);
        createParticlesInContainer(1, 4u, 0.1, 0.2);

        const bool savedH5  = Options::enableHDF5;
        Options::enableHDF5 = true;
        struct RestoreH5Option {
            bool value;
            ~RestoreH5Option() { Options::enableHDF5 = value; }
        } restoreH5{savedH5};

        const std::string base =
                (outputDirectory_m.path() / "test_diagnostic_checkpoint_rewind").string();
        const std::string f0 = base + "_c0.h5";
        const std::string f1 = base + "_c1.h5";
        auto w0              = std::make_unique<H5PartWrapperForPT>(f0, H5_O_WRONLY);
        auto w1              = std::make_unique<H5PartWrapperForPT>(f1, H5_O_WRONLY);
        std::vector<H5PartWrapper*> wrappers{w0.get(), w1.get()};

        {
            DataSink ds(wrappers, false, 2, base);
            std::vector<std::array<Vector_t<double, 3>, 2>> fd(2);
            fd[0] = zeroFdPair();
            fd[1] = zeroFdPair();

            bunch->setdT(1.0e-12);
            bunch->setT(1.0e-12);
            ds.dumpH5(*bunch, fd);
            bunch->setT(2.0e-12);
            ds.dumpH5(*bunch, fd);

            bunch->setT(1.0e-12);
            ds.rewindToCheckpoint(*bunch);
            bunch->setT(2.0e-12);
            ds.dumpH5(*bunch, fd);
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
