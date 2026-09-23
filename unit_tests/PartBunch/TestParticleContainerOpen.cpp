/**
 * @file TestParticleContainerOpen.cpp
 * @brief Verify that non‑P3M (pure spatial) particle containers honour
 *        open boundary conditions (BC::NO) when constructed via PartBunch.
 *
 * The bug: PartBunch always used PERIODIC for particleBC unless the
 *        solver type was P3M. As a result, with open field boundary
 *        conditions particles were still wrapped periodically.
 *
 * This test constructs a minimal PartBunch with open BCs (type "FFT")
 * and checks that the underlying particle container (spatial layout)
 * reports a NO boundary condition.
 */

#include <gtest/gtest.h>

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "PartBunch/BinnedFieldSolver.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"

namespace {
    using PC_t = ParticleContainer<double, 3>;
    using PartBunch_t = PartBunch<double,3>;

    class TestableFieldSolverCmd : public FieldSolverCmd {
    public:
        void setType(const std::string& t) { Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::TYPE], t); }
        void setBCX(const std::string& bc) { Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTX], bc); }
        void setBCY(const std::string& bc) { Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTY], bc); }
        void setBCZ(const std::string& bc) { Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::BCFFTZ], bc); }
        void setNX(double nx) { Attributes::setReal(this->itsAttr[FIELDSOLVER::NX], nx); }
        void setNY(double ny) { Attributes::setReal(this->itsAttr[FIELDSOLVER::NY], ny); }
        void setNZ(double nz) { Attributes::setReal(this->itsAttr[FIELDSOLVER::NZ], nz); }
        void setGreensFunction(const std::string& gf) { Attributes::setPredefinedString(this->itsAttr[FIELDSOLVER::GREENSF], gf); }
        void setParallelDecomposition(bool v) {
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTX], v);
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTY], v);
            Attributes::setBool(this->itsAttr[FIELDSOLVER::PARFFTZ], v);
        }
    };
}

TEST(OpenBoundaryParticleContainer, SpatialLayoutUsesOpenBC) {
    // Initialise IPPL once for the test suite.
    static bool ipplInit = false;
    if (!ipplInit) {
        int argc = 0; char** argv = nullptr;
        ippl::initialize(argc, argv);
        ipplInit = true;
    }

    // Build a minimal field‑solver command with open BCs.
    auto fsCmd = std::make_shared<TestableFieldSolverCmd>();
    fsCmd->setType("FFT"); // non‑P3M solver
    fsCmd->setNX(8); fsCmd->setNY(8); fsCmd->setNZ(8);
    fsCmd->setBCX("OPEN"); fsCmd->setBCY("OPEN"); fsCmd->setBCZ("OPEN");
    fsCmd->setGreensFunction("STANDARD");
    fsCmd->setParallelDecomposition(true);
    fsCmd->execute();

    // Dummy DataSink and Beam required for PartBunch construction.
    auto dataSink = std::make_shared<DataSink>();
    auto beam = std::make_shared<Beam>();
    Beam* testBeam = Beam::find("UNNAMED_BEAM");

    // Construct the bunch with a single particle.
    auto bunch = std::make_shared<PartBunch_t>(
        std::vector<double>{1.0}, // qi
        std::vector<double>{1.0}, // mi
        std::vector<Beam*>{testBeam},
        std::vector<size_t>{1}, // particles per beam
        1.0, // lbt
        "LF2",
        fsCmd.get(),
        dataSink.get()
    );

    // The particle container should be a pure spatial layout (no P3M overlap).
    auto pc = bunch->getParticleContainer();
    EXPECT_FALSE(pc->hasP3MLayout());
    // The underlying spatial layout must report a NO boundary condition.
    // getPL() returns the layout (spatial or overlap). For spatial layout the
    // boundary condition can be queried via getParticleBC().
    auto bcArray = pc->getPL().getParticleBC();
    for (auto b : bcArray) {
        EXPECT_EQ(b, ippl::BC::NO);
    }
}
