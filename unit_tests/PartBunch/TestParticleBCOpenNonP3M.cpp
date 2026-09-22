// Test to ensure open boundary conditions correctly set particle BC to NO for non-P3M solvers.

#include <gtest/gtest.h>
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "PartBunch/PartBunch.h"
#include "Structure/Beam.h"
#include "Structure/FieldSolverCmd.h"
#include "Structure/DataSink.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

namespace {
    class TestableFieldSolverCmd : public FieldSolverCmd {
    public:
        void setType(const std::string& t) { Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::TYPE], t); }
        void setBCX(const std::string& bc) { Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTX], bc); }
        void setBCY(const std::string& bc) { Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTY], bc); }
        void setBCZ(const std::string& bc) { Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTZ], bc); }
        void setNX(double nx) { Attributes::setReal(this->itsAttr[FIELDSOLVER::NX], nx); }
        void setNY(double ny) { Attributes::setReal(this->itsAttr[FIELDSOLVER::NY], ny); }
        void setNZ(double nz) { Attributes::setReal(this->itsAttr[FIELDSOLVER::NZ], nz); }
        void setParallelDecomposition(bool v) {
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTX], v);
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTY], v);
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTZ], v);
        }
    };

    class OpenBCNonP3MTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc = 0; char** argv = nullptr; ippl::initialize(argc, argv);
            OpalData::getInstance()->storeInputFn("unit_test.opal");
            gmsg = new Inform(nullptr, -1);
            Options::enableHDF5 = false;
        }
        static void TearDownTestSuite() { delete gmsg; gmsg = nullptr; ippl::finalize(); }
        void SetUp() override {
            fsCmd = std::make_shared<TestableFieldSolverCmd>();
            fsCmd->setType("NONE");
            fsCmd->setNX(8); fsCmd->setNY(8); fsCmd->setNZ(8);
            fsCmd->setBCX("OPEN"); fsCmd->setBCY("OPEN"); fsCmd->setBCZ("OPEN");
            fsCmd->setParallelDecomposition(true);
            fsCmd->execute();
            fsCmdBase = fsCmd;
            dataSink = std::make_shared<DataSink>();
            auto beam = std::make_shared<Beam>();
            Beam* testBeam = Beam::find("UNNAMED_BEAM");
            bunch = std::make_shared<PartBunch<double,3>>(std::vector<double>{1.0}, std::vector<double>{1.0},
                std::vector<Beam*>{testBeam}, std::vector<size_t>{64}, 1.0, "LF2", fsCmdBase.get(), dataSink.get());
            pc = bunch->getParticleContainer();
        }
        void TearDown() override { bunch.reset(); dataSink.reset(); fsCmd.reset(); fsCmdBase.reset(); pc.reset(); }
        std::shared_ptr<TestableFieldSolverCmd> fsCmd, fsCmdBase;
        std::shared_ptr<DataSink> dataSink;
        std::shared_ptr<PartBunch<double,3>> bunch;
        std::shared_ptr<ParticleContainer<double,3>> pc;
    };

    TEST_F(OpenBCNonP3MTest, ParticleBCIsNoForOpen) {
        // Verify that particle boundary conditions are set to NO (non-periodic)
        const auto& bc_vec = pc->getPL().getParticleBC();
        for (const auto& bc : bc_vec) {
            EXPECT_EQ(bc, ippl::BC::NO);
        }
    }
}
