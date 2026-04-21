//
// Unit tests for class VariableRFCavity
//
// Copyright (c) 2014, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
// All rights reserved.
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
#include <vector>
#include "AbsBeamline/BeamlineVisitor.h"
#include "AbsBeamline/VariableRFCavity.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "Algorithms/PolynomialTimeDependence.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "gtest/gtest.h"

class TestVariableRFCavity : public testing::Test, public VariableRFCavity, public BeamlineVisitor {
public:
    TestVariableRFCavity() : VariableRFCavity("Cavity") {}

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

    // Overrides of BeamlineVisitor
    void execute() override {}
    void visitBeamline(const Beamline&) override {}
    void visitComponent(const Component&) override {}
    void visitConstantEFieldCavity(const ConstantEFieldCavity&) override {}
    void visitDrift(const Drift&) override {}
    void visitFlaggedElmPtr(const FlaggedElmPtr&) override {}
    void visitMarker(const Marker&) override {}
    void visitMonitor(const Monitor&) override {}
    void visitMultipole(const Multipole&) override {}
    void visitMultipoleT(const MultipoleT&) override {}
    void visitRFCavity(const RFCavity&) override {}
    void visitScalingFFAMagnet(const ScalingFFAMagnet&) override {}
    void visitRing(const Ring&) override {}
    void visitSolenoid(const Solenoid&) override {}
    void visitTravelingWave(const TravelingWave&) override {}
    void visitVerticalFFAMagnet(const VerticalFFAMagnet&) override {}
    void visitProbe(const Probe&) override {}
    void visitVariableRFCavity(const VariableRFCavity&) override {}

    // Test helpers

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

    std::shared_ptr<FieldSolverCmd> fsCmdBase_m;
    std::shared_ptr<DataSink> dataSink_m;

    std::shared_ptr<PartBunch_t> makeBunch(const size_t numParticles) {
        dataSink_m       = std::make_shared<DataSink>();
        const auto fsCmd = std::make_shared<TestableFieldSolverCmd>();
        fsCmdBase_m      = fsCmd;
        fsCmd->setType("NONE");
        fsCmd->setNX(8);
        fsCmd->setNY(8);
        fsCmd->setNZ(8);
        fsCmd->setBCX("PERIODIC");
        fsCmd->setBCY("PERIODIC");
        fsCmd->setBCZ("PERIODIC");
        auto beam    = std::make_shared<Beam>();
        Beam* opBeam = Beam::find("UNNAMED_BEAM");
        auto bunch   = std::make_shared<PartBunch_t>(
                /*qi=*/std::vector{1.0}, /*mi=*/std::vector{1.0},
                /*beams=*/std::vector<Beam*>{opBeam},
                /*totalParticlesPerBeam=*/std::vector<size_t>{numParticles},
                /*lbt=*/1.0, /*integration_method=*/"LF2", fsCmdBase_m.get(), dataSink_m.get());
        bunch->getParticleContainer()->create(numParticles);
        return bunch;
    }

    // Is the obscure pointer-to-member function syntax appropriate here? I have
    // been doing too much python where this stuff is easy
    static void testGetSet(
            VariableRFCavity& cav1,
            std::shared_ptr<AbstractTimeDependence> (VariableRFCavity::*getMethod)() const,
            void (VariableRFCavity::*setMethod)(std::shared_ptr<AbstractTimeDependence>)) {
        const std::shared_ptr<AbstractTimeDependence> poly_1(
                std::make_shared<PolynomialTimeDependence>(std::vector(1, 1.)));
        const std::shared_ptr<AbstractTimeDependence> poly_2(
                std::make_shared<PolynomialTimeDependence>(std::vector(2, 2.)));

        (cav1.*setMethod)(poly_1);
        EXPECT_EQ((cav1.*getMethod)(), poly_1);  // shallow equals is okay
        (cav1.*setMethod)(poly_2);
        EXPECT_EQ((cav1.*getMethod)(), poly_2);  // shallow equals is okay
        (cav1.*setMethod)(poly_2);
        EXPECT_EQ((cav1.*getMethod)(), poly_2);  // shallow equals is okay
        (cav1.*setMethod)(nullptr);              // and this deletes the memory
    }

    static void testNull(const VariableRFCavity& cav1) {
        const std::shared_ptr<AbstractTimeDependence> null_poly(nullptr);
        EXPECT_DOUBLE_EQ(cav1.getLength(), 0.);
        EXPECT_EQ(cav1.getAmplitudeModel(), null_poly);
        EXPECT_EQ(cav1.getPhaseModel(), null_poly);
        EXPECT_EQ(cav1.getFrequencyModel(), null_poly);
    }
};

TEST_F(TestVariableRFCavity, TestConstructorEtc) {
    const VariableRFCavity cav1;
    EXPECT_EQ(cav1.getName(), "");
    testNull(cav1);
    const VariableRFCavity cav2("a_name");
    EXPECT_EQ(cav2.getName(), "a_name");
    testNull(cav1);
    // and now we implicitly check the destructor doesnt throw up on
    // case where everything is initialised to nullptr
}

TEST_F(TestVariableRFCavity, TestGetSet) {
    VariableRFCavity cav1;
    testGetSet(cav1, &VariableRFCavity::getAmplitudeModel, &VariableRFCavity::setAmplitudeModel);
    testGetSet(cav1, &VariableRFCavity::getPhaseModel, &VariableRFCavity::setPhaseModel);
    testGetSet(cav1, &VariableRFCavity::getFrequencyModel, &VariableRFCavity::setFrequencyModel);
    testNull(cav1);
    cav1.setLength(99.);
    EXPECT_DOUBLE_EQ(cav1.getLength(), 99.);
    cav1.setHeight(100.);
    EXPECT_DOUBLE_EQ(cav1.getHeight(), 100.);
    cav1.setWidth(101.);
    EXPECT_DOUBLE_EQ(cav1.getWidth(), 101.);
}

TEST_F(TestVariableRFCavity, GetTimeDependencyValues) {
    const auto amplPoly  = std::make_shared<PolynomialTimeDependence>(std::vector{1.0});
    const auto freqPoly  = std::make_shared<PolynomialTimeDependence>(std::vector{2.0});
    const auto phasePoly = std::make_shared<PolynomialTimeDependence>(std::vector{3.0});
    VariableRFCavity cav1;
    cav1.setAmplitudeModel(amplPoly);
    cav1.setFrequencyModel(freqPoly);
    cav1.setPhaseModel(phasePoly);
    EXPECT_DOUBLE_EQ(cav1.getAmplitude(0), 1.0);
    EXPECT_DOUBLE_EQ(cav1.getFrequency(0), 2.0);
    EXPECT_DOUBLE_EQ(cav1.getPhase(0), 3.0);
}

TEST_F(TestVariableRFCavity, TimeDependencyNames) {
    const auto amplPoly  = std::make_shared<PolynomialTimeDependence>(std::vector{1.0});
    const auto freqPoly  = std::make_shared<PolynomialTimeDependence>(std::vector{2.0});
    const auto phasePoly = std::make_shared<PolynomialTimeDependence>(std::vector{3.0});
    AbstractTimeDependence::setTimeDependence("AMPL", amplPoly);
    AbstractTimeDependence::setTimeDependence("FREQ", freqPoly);
    AbstractTimeDependence::setTimeDependence("PHASE", phasePoly);
    VariableRFCavity cav1;
    cav1.setAmplitudeName("AMPL");
    cav1.setFrequencyName("FREQ");
    cav1.setPhaseName("PHASE");
    cav1.setHeight(1.0);
    cav1.setWidth(1.0);
    EXPECT_NO_THROW(cav1.accept(*this));
    EXPECT_DOUBLE_EQ(cav1.getAmplitude(0), 1.0);
    EXPECT_DOUBLE_EQ(cav1.getFrequency(0), 2.0);
    EXPECT_DOUBLE_EQ(cav1.getPhase(0), 3.0);
    // Initialise failures
    cav1.setHeight(0.0);
    EXPECT_ANY_THROW(cav1.accept(*this));
    cav1.setHeight(1.0);
    cav1.setWidth(0.0);
    EXPECT_ANY_THROW(cav1.accept(*this));
}

TEST_F(TestVariableRFCavity, TestAssignmentNull) {
    const VariableRFCavity cav1;
    VariableRFCavity cav2;
    EXPECT_EQ(cav2.getLength(), 0);  // stop compiler "optimising" to copy constructor
    cav2 = cav1;
    testNull(cav2);
    const VariableRFCavity cav3(cav2);
    testNull(cav3);  // now this is really the copy constructor
}

TEST_F(TestVariableRFCavity, TestAssignmentValue) {
    const std::shared_ptr<AbstractTimeDependence> poly1(
            new PolynomialTimeDependence(std::vector(1, 1.)));
    const std::shared_ptr<AbstractTimeDependence> poly2(
            new PolynomialTimeDependence(std::vector(1, 2.)));
    const std::shared_ptr<AbstractTimeDependence> poly3(
            new PolynomialTimeDependence(std::vector(1, 3.)));
    VariableRFCavity cav1;
    cav1.setPhaseModel(poly1);
    cav1.setAmplitudeModel(poly2);
    cav1.setFrequencyModel(poly3);
    cav1.setLength(99.);
    const VariableRFCavity cav2(cav1);
    EXPECT_EQ(cav1.getPhaseModel()->getValue(1.), cav2.getPhaseModel()->getValue(1.));
    EXPECT_EQ(cav1.getAmplitudeModel()->getValue(1.), cav2.getAmplitudeModel()->getValue(1.));
    EXPECT_EQ(cav1.getFrequencyModel()->getValue(1.), cav2.getFrequencyModel()->getValue(1.));
    EXPECT_DOUBLE_EQ(cav1.getLength(), cav2.getLength());
}

TEST_F(TestVariableRFCavity, TestClone) {
    VariableRFCavity cav1;
    cav1.setLength(99.);
    const auto* cav2 = dynamic_cast<VariableRFCavity*>(cav1.clone());
    EXPECT_DOUBLE_EQ(cav1.getLength(), cav2->getLength());
    delete cav2;
}

TEST_F(TestVariableRFCavity, TestInitialiseFinalise) {
    // nothing to do here
}

TEST_F(TestVariableRFCavity, TestGetGeometry) {
    VariableRFCavity cav1;
    const VariableRFCavity& cav2(cav1);
    EXPECT_EQ(&cav1.getGeometry(), &cav2.getGeometry());
    cav1.setLength(99.);
    EXPECT_EQ(cav1.getGeometry().getElementLength(), cav1.getLength());
}

TEST_F(TestVariableRFCavity, TestApplyField) {
    VariableRFCavity cav1;
    const std::shared_ptr<AbstractTimeDependence> poly1(new PolynomialTimeDependence({1.0, 2.0}));
    const std::shared_ptr<AbstractTimeDependence> poly2(new PolynomialTimeDependence({3.0, 4.0}));
    const std::shared_ptr<AbstractTimeDependence> poly3(new PolynomialTimeDependence({5.0, 6.0}));
    cav1.setAmplitudeModel(poly1);
    cav1.setFrequencyModel(poly2);
    cav1.setPhaseModel(poly3);
    cav1.setLength(2.);
    cav1.setWidth(3.);
    cav1.setHeight(4.);
    const Vector_t<double, 3> R({1., 1., 1.});
    for (double t = 0.0; t < 10.0e-9; t += 1.0e-9) {
        Vector_t<double, 3> B({0., 0., 0.});
        Vector_t<double, 3> E({0., 0., 0.});
        const double phase     = poly3->getValue(t);
        const double amplitude = poly1->getValue(t);
        const double integralF = poly2->getIntegral(t) * Units::MHz2Hz;
        const double e_test =
                amplitude * sin(Physics::two_pi * integralF + phase) * Units::MVpm2Vpm;
        EXPECT_FALSE(cav1.apply(R, Vector_t<double, 3>(0.0), t, E, B));
        EXPECT_NEAR(0., E[0], 1.e-6);
        EXPECT_NEAR(0., E[1], 1.e-6);
        EXPECT_NEAR(e_test, E[2], 1.e-6);
        EXPECT_NEAR(0., B[0], 1.e-6);
        EXPECT_NEAR(0., B[1], 1.e-6);
        EXPECT_NEAR(0., B[2], 1.e-6);
    }
}

TEST_F(TestVariableRFCavity, TestApplyBoundingBox) {
    VariableRFCavity cav1;
    std::shared_ptr<AbstractTimeDependence> poly1(new PolynomialTimeDependence(std::vector(1, 1.)));
    std::shared_ptr<AbstractTimeDependence> poly2(new PolynomialTimeDependence(std::vector(2, 2.)));
    std::shared_ptr<AbstractTimeDependence> poly3(new PolynomialTimeDependence(std::vector(3, 3.)));
    cav1.setAmplitudeModel(poly1);
    cav1.setFrequencyModel(poly2);
    cav1.setPhaseModel(poly3);
    cav1.setLength(2.);
    cav1.setHeight(3.);
    cav1.setWidth(4.);
    Vector_t<double, 3> R({0., 0., 1.});
    Vector_t<double, 3> B({0., 0., 0.});
    Vector_t<double, 3> E({0., 0., 0.});
    double t = 0;
    EXPECT_FALSE(cav1.apply(R, Vector_t<double, 3>(0.0), t, E, B));
    R[2] = 2. - 1e-9;
    EXPECT_FALSE(cav1.apply(R, Vector_t<double, 3>(0.0), t, E, B));
    R[2] = 1.e-9;
    EXPECT_FALSE(cav1.apply(R, Vector_t<double, 3>(0.0), t, E, B));
    R[2] = -1.e-9;
    EXPECT_FALSE(cav1.apply(R, Vector_t<double, 3>(0.0), t, E, B));
    R[2] = 2. + 1.e-9;
    EXPECT_FALSE(cav1.apply(R, Vector_t<double, 3>(0.0), t, E, B));
    R[2] = 1.;
    R[1] = -1.5 - 1e-9;
    EXPECT_TRUE(cav1.apply(R, Vector_t<double, 3>(0.0), t, E, B));
    R[1] = +1.5 + 1e-9;
    EXPECT_TRUE(cav1.apply(R, Vector_t<double, 3>(0.0), t, E, B));
    R[1] = 0.;
    EXPECT_FALSE(cav1.apply(R, Vector_t<double, 3>(0.0), t, E, B));
    R[0] = -2. - 1e-9;
    EXPECT_TRUE(cav1.apply(R, Vector_t<double, 3>(0.0), t, E, B));
    R[0] = +2. + 1e-9;
    EXPECT_TRUE(cav1.apply(R, Vector_t<double, 3>(0.0), t, E, B));
    R[0] = 0.;
    EXPECT_FALSE(cav1.apply(R, Vector_t<double, 3>(0.0), t, E, B));
}

TEST_F(TestVariableRFCavity, BunchFields) {
    // Set up the cavity
    const auto amplPoly = std::make_shared<PolynomialTimeDependence>(std::vector{1.0});
    const auto freqPoly = std::make_shared<PolynomialTimeDependence>(std::vector{1.0});
    const auto phasePoly =
            std::make_shared<PolynomialTimeDependence>(std::vector{Physics::pi / 2.0});
    setAmplitudeModel(amplPoly);
    setFrequencyModel(freqPoly);
    setPhaseModel(phasePoly);
    constexpr double width  = 1.0;
    constexpr double length = 1.0;
    setLength(length);
    setHeight(2 * width);
    setWidth(2 * width);
    // Make the bunch
    std::vector<double> line(11);
    const auto bunch = makeBunch(line.size());
    const auto pc    = bunch->getParticleContainer();
    // Create the local views and data
    std::vector<Vector_t<double, 3>> localR(line.size());
    const auto hostR = Kokkos::create_mirror_view(pc->R.getView());
    const auto hostE = Kokkos::create_mirror_view(pc->E.getView());
    // Set the particle positions
    const double stepSize = width / static_cast<double>(line.size() - 1);
    for (size_t i = 0; i < line.size(); ++i) {
        localR[i] = {static_cast<double>(i) * stepSize - width / 2.0, 0.0, length / 2.0};
        hostR(i)  = {localR[i][0], localR[i][1], localR[i][2] + length / 2.0};
    }
    Kokkos::deep_copy(pc->R.getView(), hostR);
    pc->setQ(pc->getChargePerParticle());
    ippl::Comm->barrier();
    Kokkos::fence();
    // Register the bunch with the element
    bunch->setT(0.0);
    double startField, endField;
    initialise(bunch.get(), startField, endField);
    EXPECT_NE(RefPartBunch_m, nullptr);
    // Get the fields for all particles
    apply(pc);
    // Extract the fields from the GPU
    Kokkos::deep_copy(hostE, pc->E.getView());
    Kokkos::fence();
    for (size_t i = 0; i < line.size(); ++i) {
        line[i] = std::hypot(hostE(i)[0], hostE(i)[1], hostE(i)[2]);
    }
    for (size_t i = 0; i < line.size(); ++i) {
        EXPECT_DOUBLE_EQ(line[i], 1.0 * Units::MVpm2Vpm);
    }
    // Get the field for one of the particles
    Vector_t<double, 3> singleE{}, singleB{};
    EXPECT_FALSE(apply(0, 0.0, singleE, singleB));
    EXPECT_DOUBLE_EQ(singleE[2], 1.0 * Units::MVpm2Vpm);
    // Done
    finalise();
    EXPECT_EQ(RefPartBunch_m, nullptr);
}

TEST_F(TestVariableRFCavity, ReferenceParticle) {
    // Set up the cavity
    const auto amplPoly = std::make_shared<PolynomialTimeDependence>(std::vector{1.0});
    const auto freqPoly = std::make_shared<PolynomialTimeDependence>(std::vector{1.0});
    const auto phasePoly =
            std::make_shared<PolynomialTimeDependence>(std::vector{Physics::pi / 2.0});
    setAmplitudeModel(amplPoly);
    setFrequencyModel(freqPoly);
    setPhaseModel(phasePoly);
    setLength(10);
    setHeight(2.0);
    setWidth(2.0);
    const Vector_t<double, 3> R({0.0, 0.0, 5.0});
    Vector_t<double, 3> B({0.0, 0.0, 0.0});
    Vector_t<double, 3> E({0.0, 0.0, 0.0});
    EXPECT_FALSE(applyToReferenceParticle(R, {}, 0.0, E, B));
    EXPECT_DOUBLE_EQ(E[2], 1.0 * Units::MVpm2Vpm);
}

TEST_F(TestVariableRFCavity, OddApis) {
    const VariableRFCavity cav1;
    // Dimensions dummy override
    double a{}, b{};
    EXPECT_NO_THROW(cav1.getDimensions(a, b));
    // The cavity does not make a bend
    EXPECT_FALSE(cav1.bends());
    // Self assignment
    VariableRFCavity cav2;
    cav2.setLength(3.0);
    cav2 = cav1;
    EXPECT_DOUBLE_EQ(cav2.getLength(), 0.0);
    EXPECT_NO_THROW(cav2 = cav2);
    EXPECT_DOUBLE_EQ(cav2.getLength(), 0.0);
    // No implementation of field
    EXPECT_ANY_THROW(cav1.getField());
    EXPECT_ANY_THROW(cav2.getField());
}
