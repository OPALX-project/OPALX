/*
 *  Copyright (c) 2025, Jon Thompson
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */
#include <memory>
#include <hdf5.h>
#include "AbsBeamline/BeamlineVisitor.h"
#include "AbsBeamline/OutputPlane.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeConfigBuilder.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"
#include "gtest/gtest.h"

class TestOutputPlane : public testing::Test, public OutputPlane, public BeamlineVisitor {
public:
    TestOutputPlane() : OutputPlane("Magnet") {}

    class TestableFieldSolverCmd final : public FieldSolverCmd {
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

    static std::shared_ptr<FieldSolverCmd> makeNoFieldSolver() {
        const auto fsCmd = std::make_shared<TestableFieldSolverCmd>();
        fsCmd->setType("NONE");
        fsCmd->setNX(8);
        fsCmd->setNY(8);
        fsCmd->setNZ(8);
        fsCmd->setBCX("PERIODIC");
        fsCmd->setBCY("PERIODIC");
        fsCmd->setBCZ("PERIODIC");
        return fsCmd;
    }

    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        OpalData::getInstance()->storeInputFn("unit_test.opal");
        if (gmsg == nullptr) {
            gmsg = new Inform(nullptr, -1);
        }
        Options::enableHDF5 = false;
    }
    static void TearDownTestSuite() {
        if (gmsg != nullptr) {
            delete gmsg;
            gmsg = nullptr;
        }
        ippl::finalize();
    }

protected:
    void SetUp() override {
        fsCmd_m = std::make_shared<TestableFieldSolverCmd>();
        fsCmd_m->setType("NONE");
        fsCmd_m->setNX(nx);
        fsCmd_m->setNY(ny);
        fsCmd_m->setNZ(nz);
        fsCmd_m->execute();
        beam_m         = std::make_shared<Beam>();
        Beam* testBeam = Beam::find("UNNAMED_BEAM");
        bunch_m = std::make_shared<PartBunch_t>(
                /*qi=*/std::vector{1.0},
                /*mi=*/std::vector{1.0},
                /*beams=*/std::vector{testBeam},
                /*totalParticlesPerBeam=*/std::vector{kDefaultNParticles},
                /*lbt=*/1.0,
                /*integration_method=*/"LF2",
                makeCartesianDomainConfig(
                        opalx::spacecharge::buildSpaceChargeConfig(*fsCmd_m, {})));
        pc_m = bunch_m->getParticleContainer();
        tempDir_m = std::filesystem::temp_directory_path() / "TestOutputPlane";
        std::filesystem::create_directories(tempDir_m);
    }

    void TearDown() override {
        std::filesystem::remove_all(tempDir_m);
    }

public:
    // Overrides of BeamlineVisitor
    void execute() override {}
    void visitBeamline(const Beamline&) override {}
    void visitElementBase(const ElementBase&) override {}
    void visitBox(const Box&) override {}
    void visitCollimator(const Collimator&) override {}
    void visitConstantEFieldCavity(const ConstantEFieldCavity&) override {}
    void visitDrift(const Drift&) override {}
    void visitFieldmapElement(const FieldmapElement&) override {}
    void visitFlaggedElmPtr(const FlaggedElmPtr&) override {}
    void visitLaser(const Laser&) override {}
    void visitMarker(const Marker&) override {}
    void visitMonitor(const Monitor&) override {}
    void visitMultipole(const Multipole&) override {}
    void visitMultipoleT(const MultipoleT&) override {}
    void visitRBend(const RBend&) override {}
    void visitRFCavity(const RFCavity&) override {}
    void visitScalingFFAMagnet(const ScalingFFAMagnet&) override {}
    void visitSBend(const SBend&) override {}
    void visitSolenoid(const Solenoid&) override {}
    void visitTravelingWave(const TravelingWave&) override {}
    void visitVerticalFFAMagnet(const VerticalFFAMagnet&) override {}
    void visitProbe(const Probe&) override {}
    void visitVariableRFCavity(const VariableRFCavity&) override {}
    void visitEnge(const endfieldmodel::Enge&) override {}
    void visitAsymmetricEnge(const endfieldmodel::AsymmetricEnge&) override {}
    void visitTanh(const endfieldmodel::Tanh&) override {}
    void visitOutputPlane(const OutputPlane&) override {}

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

    void moveParticle(const size_t n, const Vector_t<double, 3>& r) const {
        const auto rHost = Kokkos::create_mirror_view_and_copy(pc_m->R.getView());
        rHost(n) = r;
        Kokkos::deep_copy(pc_m->R.getView(), rHost);
    }

    std::vector<double> readDataset(const std::string& filename, const std::string& groupName,
            const std::string& datasetName) {
        const hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file < 0) {
            throw std::runtime_error("Failed to open HDF5 file");
        }
        const hid_t group = H5Gopen2(file, groupName.c_str(), H5P_DEFAULT);
        if (group < 0) {
            H5Fclose(file);
            throw std::runtime_error("Failed to open group");
        }
        const hid_t dataset = H5Dopen2(group, datasetName.c_str(), H5P_DEFAULT);
        if (dataset < 0) {
            H5Gclose(group);
            H5Fclose(file);
            throw std::runtime_error("Failed to open dataset");
        }
        const hid_t dataspace = H5Dget_space(dataset);
        const int rank = H5Sget_simple_extent_ndims(dataspace);
        if (rank != 1) {
            H5Sclose(dataspace);
            H5Dclose(dataset);
            H5Gclose(group);
            H5Fclose(file);
            throw std::runtime_error("Dataset is not 1-dimensional");
        }
        hsize_t dims[1];
        H5Sget_simple_extent_dims(dataspace, dims, nullptr);
        std::vector<double> data(dims[0]);
        H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Gclose(group);
        H5Fclose(file);
        return data;
    }

    static constexpr size_t kDefaultNParticles = 64;
    static constexpr double nx                 = 8;
    static constexpr double ny                 = 9;
    static constexpr double nz                 = 10;
    std::shared_ptr<TestableFieldSolverCmd> fsCmd_m;
    std::shared_ptr<Beam> beam_m;
    std::shared_ptr<PartBunch_t> bunch_m;
    std::shared_ptr<ParticleContainer_t> pc_m;
    std::filesystem::path tempDir_m;
};

// Does the clone API return an object with the same configuration
TEST_F(TestOutputPlane, Cloning) {
    // Set up the output plane
    setNormal({1.0, 2.0, 3.0});
    setCentre({4.0, 5.0, 6.0});
    setRadialExtent(7.0);
    // Does the clone have the same configuration
    std::unique_ptr<OutputPlane> theClone;
    theClone.reset(dynamic_cast<OutputPlane*>(clone()));
    EXPECT_NE(theClone.get(), nullptr);
    EXPECT_EQ(theClone->getNormal()[0], getNormal()[0]);
    EXPECT_EQ(theClone->getNormal()[1], getNormal()[1]);
    EXPECT_EQ(theClone->getNormal()[2], getNormal()[2]);
    EXPECT_EQ(theClone->getCentre()[0], getCentre()[0]);
    EXPECT_EQ(theClone->getCentre()[1], getCentre()[1]);
    EXPECT_EQ(theClone->getCentre()[2], getCentre()[2]);
    EXPECT_EQ(theClone->getRadialExtent(), getRadialExtent());
}

// Does a clone end up in the global list
TEST_F(TestOutputPlane, ListBuilding) {
    // Set up the output plane
    setNormal({1.0, 2.0, 3.0});
    setCentre({4.0, 5.0, 6.0});
    setRadialExtent(7.0);
    // Does the clone have the same configuration
    std::unique_ptr<OutputPlane> theClone(dynamic_cast<OutputPlane*>(clone()));
    EXPECT_EQ(outputPlanes_m.size(), 1);
}

TEST_F(TestOutputPlane, InitialParticleStashing) {
    auto fName = (tempDir_m / "capture.h5").string();
    // Set up the output plane
    setNormal({1.0, 0.0, 0.0});
    setCentre({1.0, 0.0, 0.0});
    setRadialExtent(1.0);
    setOutputFN(fName);
    const std::unique_ptr<OutputPlane> theClone(dynamic_cast<OutputPlane*>(clone()));
    // Create a particle that's before the plane
    createParticles({{0.9, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {false});
    // Operate the output plane
    Options::enableHDF5 = true;
    theClone->initialise(bunch_m.get());
    observe(bunch_m.get());
    theClone->goOffline();
    theClone->finalise();
    // Check the particles are stashed
    ASSERT_EQ(theClone->stash().size(), 1);
    EXPECT_EQ(theClone->stash()[0]->r().size(), 1);
    EXPECT_EQ(theClone->stash()[0]->p().size(), 1);
    EXPECT_EQ(theClone->stash()[0]->crossed().size(), 1);
    EXPECT_EQ(theClone->stash()[0]->t().size(), 1);
    const auto hostR = Kokkos::create_mirror_view_and_copy(theClone->stash()[0]->r());
    EXPECT_DOUBLE_EQ(hostR(0)[0], 0.9);
    EXPECT_DOUBLE_EQ(hostR(0)[1], 0.0);
    EXPECT_DOUBLE_EQ(hostR(0)[2], 0.0);
    // There should be no file yet
    EXPECT_FALSE(std::filesystem::exists(fName));
}

TEST_F(TestOutputPlane, BasicCrossingOneParticle) {
    auto fName = (tempDir_m / "capture.h5").string();
    // Set up the output plane
    setNormal({1.0, 0.0, 0.0});
    setCentre({1.0, 0.0, 0.0});
    setRadialExtent(1.0);
    setOutputFN(fName);
    const std::unique_ptr<OutputPlane> theClone(dynamic_cast<OutputPlane*>(clone()));
    // Create a particle that's before the plane
    createParticles({{0.9, 0.0, 0.0}}, {{1.0, 2.0, 3.0}}, {false});
    // Operate the output plane
    Options::enableHDF5 = true;
    theClone->initialise(bunch_m.get());
    observe(bunch_m.get());
    // Move the particle
    moveParticle(0, {1.1, 0.0, 0.0});
    observe(bunch_m.get());
    theClone->goOffline();
    theClone->finalise();
    // Check the particles are stashed
    ASSERT_EQ(theClone->stash().size(), 1);
    EXPECT_EQ(theClone->stash()[0]->r().size(), 1);
    EXPECT_EQ(theClone->stash()[0]->p().size(), 1);
    EXPECT_EQ(theClone->stash()[0]->crossed().size(), 1);
    EXPECT_EQ(theClone->stash()[0]->t().size(), 1);
    const auto hostR = Kokkos::create_mirror_view_and_copy(theClone->stash()[0]->r());
    EXPECT_DOUBLE_EQ(hostR(0)[0], 1.1);
    EXPECT_DOUBLE_EQ(hostR(0)[1], 0.0);
    EXPECT_DOUBLE_EQ(hostR(0)[2], 0.0);
    // There should now be a file
    EXPECT_TRUE(std::filesystem::exists(fName));
    // Check the position
    auto x = readDataset(fName, "Step#0", "x");
    auto y = readDataset(fName, "Step#0", "y");
    auto z = readDataset(fName, "Step#0", "z");
    EXPECT_EQ(x.size(), 1);
    EXPECT_EQ(y.size(), 1);
    EXPECT_EQ(z.size(), 1);
    EXPECT_DOUBLE_EQ(x[0], 1.0);
    EXPECT_DOUBLE_EQ(y[0], 0.0);
    EXPECT_DOUBLE_EQ(z[0], 0.0);
    // Check the momentum
    auto px = readDataset(fName, "Step#0", "px");
    auto py = readDataset(fName, "Step#0", "py");
    auto pz = readDataset(fName, "Step#0", "pz");
    EXPECT_EQ(px.size(), 1);
    EXPECT_EQ(py.size(), 1);
    EXPECT_EQ(pz.size(), 1);
    EXPECT_DOUBLE_EQ(px[0], 1.0);
    EXPECT_DOUBLE_EQ(py[0], 2.0);
    EXPECT_DOUBLE_EQ(pz[0], 3.0);
}

TEST_F(TestOutputPlane, BasicCrossingMultipleParticles) {
    auto fName = (tempDir_m / "capture.h5").string();
    // Set up the output plane
    setNormal({1.0, 0.0, 0.0});
    setCentre({1.0, 0.0, 0.0});
    setRadialExtent(1.0);
    setOutputFN(fName);
    const std::unique_ptr<OutputPlane> theClone(dynamic_cast<OutputPlane*>(clone()));
    // Create a particle that's before the plane
    createParticles({{0.9, 0.0, 0.0}, {1.0, 0.2, 0.3}}, {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
            {false, false});
    // Operate the output plane
    Options::enableHDF5 = true;
    theClone->initialise(bunch_m.get());
    observe(bunch_m.get());
    // Move the particle
    moveParticle(0, {1.1, 0.0, 0.0});
    moveParticle(1, {1.2, 0.2, 0.3});
    observe(bunch_m.get());
    theClone->goOffline();
    theClone->finalise();
    // Check the particles are stashed
    ASSERT_EQ(theClone->stash().size(), 1);
    EXPECT_EQ(theClone->stash()[0]->r().size(), 2);
    EXPECT_EQ(theClone->stash()[0]->p().size(), 2);
    EXPECT_EQ(theClone->stash()[0]->crossed().size(), 2);
    EXPECT_EQ(theClone->stash()[0]->t().size(), 2);
    // There should now be a file
    EXPECT_TRUE(std::filesystem::exists(fName));
    auto x = readDataset(fName, "Step#0", "x");
    EXPECT_EQ(x.size(), 2);
    EXPECT_DOUBLE_EQ(x[0], 1.0);
    EXPECT_DOUBLE_EQ(x[1], 1.0);
    auto y = readDataset(fName, "Step#0", "y");
    EXPECT_EQ(y.size(), 2);
    EXPECT_DOUBLE_EQ(y[0], 0.0);
    EXPECT_DOUBLE_EQ(y[1], 0.2);
    auto z = readDataset(fName, "Step#0", "z");
    EXPECT_EQ(z.size(), 2);
    EXPECT_DOUBLE_EQ(z[0], 0.0);
    EXPECT_DOUBLE_EQ(z[1], 0.3);
}

TEST_F(TestOutputPlane, BasicCrossingInvalidParticle) {
    auto fName = (tempDir_m / "capture.h5").string();
    // Set up the output plane
    setNormal({1.0, 0.0, 0.0});
    setCentre({1.0, 0.0, 0.0});
    setRadialExtent(1.0);
    setOutputFN(fName);
    const std::unique_ptr<OutputPlane> theClone(dynamic_cast<OutputPlane*>(clone()));
    // Create a particle that's before the plane
    createParticles({{0.9, 0.0, 0.0}}, {{1.0, 2.0, 3.0}}, {true});
    // Operate the output plane
    Options::enableHDF5 = true;
    theClone->initialise(bunch_m.get());
    observe(bunch_m.get());
    // Move the particle
    moveParticle(0, {1.1, 0.0, 0.0});
    observe(bunch_m.get());
    theClone->goOffline();
    theClone->finalise();
    // Check the particles are stashed
    ASSERT_EQ(theClone->stash().size(), 1);
    EXPECT_EQ(theClone->stash()[0]->r().size(), 1);
    EXPECT_EQ(theClone->stash()[0]->p().size(), 1);
    EXPECT_EQ(theClone->stash()[0]->crossed().size(), 1);
    EXPECT_EQ(theClone->stash()[0]->t().size(), 1);
    const auto hostR = Kokkos::create_mirror_view_and_copy(theClone->stash()[0]->r());
    EXPECT_DOUBLE_EQ(hostR(0)[0], 1.1);
    EXPECT_DOUBLE_EQ(hostR(0)[1], 0.0);
    EXPECT_DOUBLE_EQ(hostR(0)[2], 0.0);
    // There should not be a file
    EXPECT_FALSE(std::filesystem::exists(fName));
}
