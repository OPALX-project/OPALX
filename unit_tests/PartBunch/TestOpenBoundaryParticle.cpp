/**
 * @file TestOpenBoundaryParticle.cpp
 * @brief Regression test for particle boundary handling with open field solver (non-P3M).
 *
 * Constructs a PartBunch with an OPEN field solver (type "OPEN") and no binning.
 * Inserts a single particle far outside the computational domain and checks that
 * `markParticlesOutside` correctly marks it as invalid. With the previous bug the
 * particle would have been wrapped periodically, staying inside the domain and
 * not be marked, causing the test to fail.
 */

#include <gtest/gtest.h>
#include <mpi.h>

#include "AbstractObjects/OpalData.h"
#include "Ippl.h"
#include "PartBunch/PartBunch.h"
#include "Structure/Beam.h"
#include "Structure/FieldSolverCmd.h"
#include "Structure/DataSink.h"
#include "Utilities/Options.h"

namespace {
    using PartBunch_t = PartBunch<double, 3>;
    using ParticleContainer_t = typename PartBunch_t::ParticleContainer_t;

    class TestableFieldSolverCmd : public FieldSolverCmd {
    public:
        void setType(const std::string& t) { Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::TYPE], t); }
        void setBCX(const std::string& bc) { Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTX], bc); }
        void setBCY(const std::string& bc) { Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTY], bc); }
        void setBCZ(const std::string& bc) { Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTZ], bc); }
        void setNX(double v) { Attributes::setReal(this->itsAttr[FIELDSOLVER::NX], v); }
        void setNY(double v) { Attributes::setReal(this->itsAttr[FIELDSOLVER::NY], v); }
        void setNZ(double v) { Attributes::setReal(this->itsAttr[FIELDSOLVER::NZ], v); }
    };

    class OpenBoundaryTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc = 0; char** argv = nullptr; ippl::initialize(argc, argv);
            OpalData::getInstance()->storeInputFn("unit_test.opal");
            Options::enableHDF5 = false;
        }
        static void TearDownTestSuite() { ippl::finalize(); }

        void SetUp() override {
            // Minimal field solver command: OPEN type, no binning.
            fsCmd = std::make_shared<TestableFieldSolverCmd>();
            fsCmd->setType("OPEN");
            fsCmd->setBCX("OPEN");
            fsCmd->setBCY("OPEN");
            fsCmd->setBCZ("OPEN");
            fsCmd->setNX(8.0); fsCmd->setNY(8.0); fsCmd->setNZ(8.0);
            fsCmd->execute();
            fsCmdBase = fsCmd;

            dataSink = std::make_shared<DataSink>();
            auto beam = std::make_shared<Beam>();
            Beam* testBeam = Beam::find("UNNAMED_BEAM");

            // Create PartBunch with a single particle container.
            bunch = std::make_shared<PartBunch_t>(
                std::vector<double>{1.0},
                std::vector<double>{1.0},
                std::vector<Beam*>{testBeam},
                std::vector<size_t>{1},
                1.0,
                "LF2",
                fsCmdBase.get(), dataSink.get());
            pc = bunch->getParticleContainer();
        }
        void TearDown() override {
            pc.reset(); bunch.reset(); dataSink.reset(); fsCmd.reset(); fsCmdBase.reset();
        }
        std::shared_ptr<TestableFieldSolverCmd> fsCmd;
        std::shared_ptr<FieldSolverCmd> fsCmdBase;
        std::shared_ptr<DataSink> dataSink;
        std::shared_ptr<PartBunch_t> bunch;
        std::shared_ptr<ParticleContainer_t> pc;
    };
}

TEST_F(OpenBoundaryTest, ParticleOutsideIsMarked) {
    // Create a particle far outside the domain (e.g., at x=100).
    pc->createParticles(1);
    auto R_host = pc->R.getHostMirror();
    R_host(0)[0] = 100.0; // x far outside
    R_host(0)[1] = 0.0;
    R_host(0)[2] = 0.0;
    Kokkos::deep_copy(pc->R.getView(), R_host);
    // Mark particles outside with zero sigma threshold (any deviation from mean is outside).
    size_t marked = pc->markParticlesOutside(0.0);
    EXPECT_EQ(marked, 1u);
}
