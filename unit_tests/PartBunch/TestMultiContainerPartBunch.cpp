/**
 * @file TestMultiContainerPartBunch.cpp
 * @brief Unit tests for multi-container PartBunch construction, indexed access,
 *        aggregation, and DataSink diagnostic contracts (stems, fdext sizing).
 */

#include <gtest/gtest.h>

#include <cstdio>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "PartBunch/PartBunch.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Structure/H5PartWrapperForPT.h"
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
    };

    using PartBunch_t         = PartBunch<double, 3>;
    using ParticleContainer_t = typename PartBunch_t::ParticleContainer_t;

    constexpr size_t kParticlesPerBeam = 32;

    class MultiContainerPartBunchTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);

            OpalData::getInstance()->storeInputFn("unit_test.opal");
            gmsg                = new Inform(nullptr, -1);
            Options::enableHDF5 = false;
        }

        static void TearDownTestSuite() {
            delete gmsg;
            gmsg = nullptr;
            ippl::finalize();
        }

        void SetUp() override {
            fsCmd = std::make_shared<TestableFieldSolverCmd>();
            fsCmd->setType("NONE");
            fsCmd->setNX(8);
            fsCmd->setNY(8);
            fsCmd->setNZ(8);
            fsCmd->setBCX("PERIODIC");
            fsCmd->setBCY("PERIODIC");
            fsCmd->setBCZ("PERIODIC");
            fsCmdBase = fsCmd;

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
                    /*integration_method=*/"LF2", fsCmdBase.get(), dataSink.get());
            q0_m = q0;
            q1_m = q1;
            m0_m = m0;
            m1_m = m1;
        }

        void TearDown() override {
            bunch.reset();
            dataSink.reset();
            fsCmd.reset();
            fsCmdBase.reset();
            beam.reset();
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

        std::shared_ptr<TestableFieldSolverCmd> fsCmd;
        std::shared_ptr<FieldSolverCmd> fsCmdBase;
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
                        fsCmdBase.get(), dataSink.get())),
                OpalException);
    }

    TEST_F(MultiContainerPartBunchTest, Constructor_ThrowsOnNullBeam) {
        EXPECT_THROW(
                static_cast<void>(std::make_shared<PartBunch_t>(
                        std::vector<double>{1.0, 1.0}, std::vector<double>{1.0, 1.0},
                        std::vector<Beam*>{testBeam, nullptr},
                        std::vector<size_t>{kParticlesPerBeam, kParticlesPerBeam}, 1.0, "LF2",
                        fsCmdBase.get(), dataSink.get())),
                OpalException);
    }

    TEST_F(MultiContainerPartBunchTest, Constructor_ThrowsOnNullDataSink) {
        EXPECT_THROW(
                static_cast<void>(std::make_shared<PartBunch_t>(
                        std::vector<double>{1.0, 1.0}, std::vector<double>{1.0, 1.0},
                        std::vector<Beam*>{testBeam, testBeam},
                        std::vector<size_t>{kParticlesPerBeam, kParticlesPerBeam}, 1.0, "LF2",
                        fsCmdBase.get(), nullptr)),
                OpalException);
    }

    TEST_F(MultiContainerPartBunchTest, TotalNumAllContainers_SumsBoth) {
        const size_t n0 = 11u;
        const size_t n1 = 7u;
        createParticlesInContainer(0, n0, 0.1, 0.5);
        createParticlesInContainer(1, n1, 0.2, 0.6);

        EXPECT_EQ(bunch->getTotalNumAllContainers(), n0 + n1);
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
        std::remove((OpalData::getInstance()->getInputBasename() + ".lbal").c_str());
    }

    TEST_F(MultiContainerPartBunchTest, DataSink_dumpSDDS_ThrowsWhenFdextTooSmall) {
        createParticlesInContainer(0, 8u, 0.1, 0.3);
        createParticlesInContainer(1, 6u, 0.1, 0.3);

        DataSink ds(std::vector<H5PartWrapper*>{}, false, 2);
        std::vector<std::array<Vector_t<double, 3>, 2>> fdTooSmall;
        fdTooSmall.push_back(zeroFdPair());

        const double azimuth = 0.0;
        EXPECT_THROW(ds.dumpSDDS(*bunch, fdTooSmall, azimuth), OpalException);

        std::remove((DataSink::diagnosticStemForContainer("unit_test", 2, 0) + ".stat").c_str());
        std::remove((DataSink::diagnosticStemForContainer("unit_test", 2, 1) + ".stat").c_str());
        std::remove((OpalData::getInstance()->getInputBasename() + ".lbal").c_str());
    }

    TEST_F(MultiContainerPartBunchTest, DataSink_dumpSDDS_NoThrowWithMatchingFdextAndHdf5Off) {
        createParticlesInContainer(0, 8u, 0.1, 0.3);
        createParticlesInContainer(1, 6u, 0.1, 0.3);

        DataSink ds(std::vector<H5PartWrapper*>{}, false, 2);
        std::vector<std::array<Vector_t<double, 3>, 2>> fd(2);
        fd[0] = zeroFdPair();
        fd[1] = zeroFdPair();

        const double azimuth = 0.0;
        EXPECT_NO_THROW(ds.dumpSDDS(*bunch, fd, azimuth));

        std::remove((DataSink::diagnosticStemForContainer("unit_test", 2, 0) + ".stat").c_str());
        std::remove((DataSink::diagnosticStemForContainer("unit_test", 2, 1) + ".stat").c_str());
        std::remove((OpalData::getInstance()->getInputBasename() + ".lbal").c_str());
    }

    // When HDF5 is disabled, dumpH5 returns before validating fdext size.
    TEST_F(MultiContainerPartBunchTest, DataSink_dumpH5_NoThrowWhenHdf5OffEvenIfFdextTooSmall) {
        ASSERT_FALSE(Options::enableHDF5);

        createParticlesInContainer(0, 4u, 0.1, 0.2);

        DataSink ds(std::vector<H5PartWrapper*>{}, false, 2);
        std::vector<std::array<Vector_t<double, 3>, 2>> fdTooSmall;
        fdTooSmall.push_back(zeroFdPair());

        EXPECT_NO_THROW(ds.dumpH5(*bunch, fdTooSmall));

        std::remove((DataSink::diagnosticStemForContainer("unit_test", 2, 0) + ".stat").c_str());
        std::remove((DataSink::diagnosticStemForContainer("unit_test", 2, 1) + ".stat").c_str());
        std::remove((OpalData::getInstance()->getInputBasename() + ".lbal").c_str());
    }

    TEST_F(MultiContainerPartBunchTest, DataSink_dumpH5_ThrowsWhenFdextTooSmallAndHdf5On) {
        createParticlesInContainer(0, 4u, 0.1, 0.2);

        const bool savedH5  = Options::enableHDF5;
        Options::enableHDF5 = true;

        static int h5FileCounter = 0;
        const int tag            = ++h5FileCounter;
        const std::string f0     = std::string("test_mc_h5_c0_") + std::to_string(tag) + ".h5";
        const std::string f1     = std::string("test_mc_h5_c1_") + std::to_string(tag) + ".h5";

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

        std::remove(f0.c_str());
        std::remove(f1.c_str());
        std::remove((DataSink::diagnosticStemForContainer("unit_test", 2, 0) + ".stat").c_str());
        std::remove((DataSink::diagnosticStemForContainer("unit_test", 2, 1) + ".stat").c_str());
        std::remove((OpalData::getInstance()->getInputBasename() + ".lbal").c_str());

        Options::enableHDF5 = savedH5;
    }

}  // namespace
