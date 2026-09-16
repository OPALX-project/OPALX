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
#include <cmath>
#include <limits>
#include <vector>
#include "AbsBeamline/BeamlineVisitor.h"
#include "AbsBeamline/VariableRFCavity.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "Algorithms/PolynomialTimeDependence.h"
#include "Algorithms/ParallelTracker.h"
#include "Beamlines/FlaggedBeamline.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
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
    void visitElementBase(const ElementBase&) override {}
    void visitCollimator(const Collimator&) override {}
    void visitConstantEFieldCavity(const ConstantEFieldCavity&) override {}
    void visitDrift(const Drift&) override {}
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
                /*lbt=*/1.0, /*integration_method=*/"LF2",
                opalx::spacecharge::CartesianDomainConfig3D{.periodicParticleBoundary = true});
        bunch->getParticleContainer()->createParticles(numParticles);
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
        EXPECT_FALSE(cav1.applyToReferenceParticle(R, Vector_t<double, 3>(0.0), t, E, B));
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
    EXPECT_FALSE(cav1.applyToReferenceParticle(R, Vector_t<double, 3>(0.0), t, E, B));
    R[2] = 2. - 1e-9;
    EXPECT_FALSE(cav1.applyToReferenceParticle(R, Vector_t<double, 3>(0.0), t, E, B));
    R[2] = 1.e-9;
    EXPECT_FALSE(cav1.applyToReferenceParticle(R, Vector_t<double, 3>(0.0), t, E, B));
    R[2] = -1.e-9;
    EXPECT_FALSE(cav1.applyToReferenceParticle(R, Vector_t<double, 3>(0.0), t, E, B));
    R[2] = 2. + 1.e-9;
    EXPECT_FALSE(cav1.applyToReferenceParticle(R, Vector_t<double, 3>(0.0), t, E, B));
    R[2] = 1.;
    R[1] = -1.5 - 1e-9;
    EXPECT_TRUE(cav1.applyToReferenceParticle(R, Vector_t<double, 3>(0.0), t, E, B));
    R[1] = +1.5 + 1e-9;
    EXPECT_TRUE(cav1.applyToReferenceParticle(R, Vector_t<double, 3>(0.0), t, E, B));
    R[1] = 0.;
    EXPECT_FALSE(cav1.applyToReferenceParticle(R, Vector_t<double, 3>(0.0), t, E, B));
    R[0] = -2. - 1e-9;
    EXPECT_TRUE(cav1.applyToReferenceParticle(R, Vector_t<double, 3>(0.0), t, E, B));
    R[0] = +2. + 1e-9;
    EXPECT_TRUE(cav1.applyToReferenceParticle(R, Vector_t<double, 3>(0.0), t, E, B));
    R[0] = 0.;
    EXPECT_FALSE(cav1.applyToReferenceParticle(R, Vector_t<double, 3>(0.0), t, E, B));
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
        hostR(i)  = localR[i];
    }
    Kokkos::deep_copy(pc->R.getView(), hostR);
    pc->setQ(pc->getChargePerParticle());
    ippl::Comm->barrier();
    Kokkos::fence();
    // Register the bunch with the element
    bunch->setT(0.0);
    bunch->setdT(0.0); // This existing peak-field check samples exactly t=0.
    initialise(bunch.get());
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
    // Get the field at one particle's position via the position overload
    Vector_t<double, 3> singleE{}, singleB{};
    const auto hostR0 = Kokkos::create_mirror_view(pc->R.getView());
    Kokkos::deep_copy(hostR0, pc->R.getView());
    EXPECT_FALSE(applyToReferenceParticle(hostR0(0), Vector_t<double, 3>{}, 0.0, singleE, singleB));
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
    // The field-support interval follows the body length.
    double a{}, b{};
    EXPECT_NO_THROW(cav1.getFieldExtent(a, b));
    EXPECT_DOUBLE_EQ(a, 0.0);
    EXPECT_DOUBLE_EQ(b, 0.0);
    VariableRFCavity cavWithLength;
    cavWithLength.setLength(3.0);
    EXPECT_NO_THROW(cavWithLength.getFieldExtent(a, b));
    EXPECT_DOUBLE_EQ(a, 0.0);
    EXPECT_DOUBLE_EQ(b, 3.0);
    // The cavity does not make a bend
    // Self assignment
    VariableRFCavity cav2;
    cav2.setLength(3.0);
    cav2 = cav1;
    EXPECT_DOUBLE_EQ(cav2.getLength(), 0.0);
    EXPECT_NO_THROW(cav2 = cav2);
    EXPECT_DOUBLE_EQ(cav2.getLength(), 0.0);
}

TEST_F(TestVariableRFCavity, FieldSupportMatchesBodyLength) {
    const auto amplPoly = std::make_shared<PolynomialTimeDependence>(std::vector{1.0});
    const auto freqPoly = std::make_shared<PolynomialTimeDependence>(std::vector{1.0});
    const auto phasePoly =
            std::make_shared<PolynomialTimeDependence>(std::vector{Physics::pi / 2.0});
    setAmplitudeModel(amplPoly);
    setFrequencyModel(freqPoly);
    setPhaseModel(phasePoly);
    setLength(10.0);
    setHeight(2.0);
    setWidth(2.0);

    double zBegin = -1.0;
    double zEnd   = -1.0;
    getFieldExtent(zBegin, zEnd);
    EXPECT_DOUBLE_EQ(zBegin, 0.0);
    EXPECT_DOUBLE_EQ(zEnd, 10.0);

    Vector_t<double, 3> E{}, B{};
    EXPECT_FALSE(applyToReferenceParticle({0.0, 0.0, -0.1}, {}, 0.0, E, B));
    EXPECT_DOUBLE_EQ(E[2], 0.0);
    EXPECT_FALSE(applyToReferenceParticle({0.0, 0.0, 5.0}, {}, 0.0, E, B));
    EXPECT_DOUBLE_EQ(E[2], 1.0 * Units::MVpm2Vpm);
}

TEST_F(TestVariableRFCavity, DeviceSamplesAllTimeModelsAtBorisMidpoint) {
    // Nonconstant models expose sampling at the start, at the end, or with
    // frequency*f(t) instead of the required integral of frequency.
    setAmplitudeModel(std::make_shared<PolynomialTimeDependence>(std::vector{1.25, 2.e7}));
    setFrequencyModel(std::make_shared<PolynomialTimeDependence>(std::vector{5., 1.e8}));
    setPhaseModel(std::make_shared<PolynomialTimeDependence>(std::vector{0.2, -2.e6}));
    setLength(0.1);
    setWidth(0.2);
    setHeight(0.4);
    const auto bunch = makeBunch(1);
    const auto pc = bunch->getParticleContainer();
    auto r = Kokkos::create_mirror_view(pc->R.getView());
    auto e = Kokkos::create_mirror_view(pc->E.getView());
    auto b = Kokkos::create_mirror_view(pc->B.getView());
    const Vector_t<double, 3> initialE{2., -3., 4.}, initialB{0.1, -0.2, 0.3};
    r(0) = {0.01, -0.02, 0.05};
    Kokkos::deep_copy(pc->R.getView(), r);
    initialise(bunch.get());
    for (const double dt : {4.e-9, 2.e-9}) {
        SCOPED_TRACE(dt);
        constexpr double t = 7.e-9;
        bunch->setT(t);
        bunch->setdT(dt);
        e(0) = initialE;
        b(0) = initialB;
        Kokkos::deep_copy(pc->E.getView(), e);
        Kokkos::deep_copy(pc->B.getView(), b);
        apply(pc);
        Kokkos::deep_copy(e, pc->E.getView());
        Kokkos::deep_copy(b, pc->B.getView());
        const double tm = t + 0.5 * dt;
        const double expected = (1.25 + 2.e7 * tm) * 1.e6
                * std::sin(Physics::two_pi * 1.e6 * (5. * tm + 0.5e8 * tm * tm)
                           + 0.2 - 2.e6 * tm);
        auto hostE = initialE, refE = initialE, hostB = initialB, refB = initialB;
        apply(r(0), {}, tm, hostE, hostB);
        EXPECT_FALSE(applyToReferenceParticle(r(0), {}, tm, refE, refB));
        EXPECT_NEAR(e(0)[2], initialE[2] + expected,
                    16 * std::numeric_limits<double>::epsilon() * std::abs(expected));
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_DOUBLE_EQ(e(0)[d], hostE[d]);
            EXPECT_DOUBLE_EQ(e(0)[d], refE[d]);
            EXPECT_DOUBLE_EQ(b(0)[d], initialB[d]);
            EXPECT_DOUBLE_EQ(hostB[d], initialB[d]);
            EXPECT_DOUBLE_EQ(refB[d], initialB[d]);
        }
        EXPECT_DOUBLE_EQ(e(0)[0], initialE[0]);
        EXPECT_DOUBLE_EQ(e(0)[1], initialE[1]);
    }
    finalise();
}

TEST_F(TestVariableRFCavity, DeviceAndHostHaveIdenticalRectangularSupport) {
    setAmplitudeModel(std::make_shared<PolynomialTimeDependence>(std::vector{1.}));
    setFrequencyModel(std::make_shared<PolynomialTimeDependence>(std::vector{1.}));
    setPhaseModel(std::make_shared<PolynomialTimeDependence>(std::vector{Physics::pi / 2.}));
    setLength(0.5);
    setWidth(0.25);
    setHeight(0.75);
    struct Sample { Vector_t<double, 3> r; bool field, material; };
    const std::vector<Sample> samples{
        {{0., 0., 0.}, true, false},
        {{0., 0., std::nextafter(0., -1.)}, false, false},
        {{0., 0., std::nextafter(0.5, 0.)}, true, false},
        {{0., 0., 0.5}, false, false},
        {{0.125, 0.375, 0.25}, true, false},
        {{-0.125, -0.375, 0.25}, true, false},
        {{std::nextafter(0.125, 1.), 0., 0.25}, false, true},
        {{std::nextafter(-0.125, -1.), 0., 0.25}, false, true},
        {{0., std::nextafter(0.375, 1.), 0.25}, false, true},
        {{0., std::nextafter(-0.375, -1.), 0.25}, false, true},
        // Another arm of a ring can have the same local z but be metres away.
        {{8., 0., 0.25}, false, true},
        {{8., 0., 0.5}, false, false}};
    const auto bunch = makeBunch(samples.size());
    const auto pc = bunch->getParticleContainer();
    auto r = Kokkos::create_mirror_view(pc->R.getView());
    auto e = Kokkos::create_mirror_view(pc->E.getView());
    auto b = Kokkos::create_mirror_view(pc->B.getView());
    const Vector_t<double, 3> initialE{2., -3., 4.}, initialB{0.1, -0.2, 0.3};
    for (size_t i = 0; i < samples.size(); ++i) {
        r(i) = samples[i].r;
        e(i) = initialE;
        b(i) = initialB;
    }
    Kokkos::deep_copy(pc->R.getView(), r);
    Kokkos::deep_copy(pc->E.getView(), e);
    Kokkos::deep_copy(pc->B.getView(), b);
    bunch->setT(-1.e-9);
    bunch->setdT(2.e-9); // Peak field at the physical midpoint t=0.
    initialise(bunch.get());
    apply(pc);
    Kokkos::deep_copy(e, pc->E.getView());
    Kokkos::deep_copy(b, pc->B.getView());
    for (size_t i = 0; i < samples.size(); ++i) {
        SCOPED_TRACE(i);
        EXPECT_EQ(isInside(samples[i].r), samples[i].field);
        auto hostE = initialE, refE = initialE, hostB = initialB, refB = initialB;
        apply(samples[i].r, {}, 0., hostE, hostB);
        EXPECT_EQ(applyToReferenceParticle(samples[i].r, {}, 0., refE, refB), samples[i].material);
        const Vector_t<double, 3> expectedE{2., -3., samples[i].field ? 1000004. : 4.};
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_DOUBLE_EQ(e(i)[d], expectedE[d]);
            EXPECT_DOUBLE_EQ(hostE[d], expectedE[d]);
            EXPECT_DOUBLE_EQ(refE[d], expectedE[d]);
            EXPECT_DOUBLE_EQ(b(i)[d], initialB[d]);
            EXPECT_DOUBLE_EQ(hostB[d], initialB[d]);
            EXPECT_DOUBLE_EQ(refB[d], initialB[d]);
        }
    }
    finalise();
}

TEST_F(TestVariableRFCavity, ParallelTrackerVisitorRegistersLiveCavity) {
    // Exercise real visitor dispatch without executing an OPALX tracking run.
    // Previously the inherited DefaultVisitor handler silently dropped this
    // element, despite the cavity's direct host/device field tests passing.
    struct Observation {
        unsigned clones = 0;
        PartBunch_t* attached = nullptr;
        VariableRFCavity* runtime = nullptr;
    } observed;
    struct ObservedCavity final : VariableRFCavity {
        explicit ObservedCavity(Observation& state) : VariableRFCavity("registeredRF"), state(state) {}
        ElementBase* clone() const override {
            ++state.clones;
            auto* copy = new ObservedCavity(*this);
            state.runtime = copy;
            return copy;
        }
        void initialise(PartBunch_t* bunch) override {
            VariableRFCavity::initialise(bunch);
            state.attached = bunch;
        }
        Observation& state;
    } cavity(observed);
    struct BorrowedBunchTracker final : ParallelTracker {
        BorrowedBunchTracker(const Beamline& line, PartBunch_t& bunch)
            : ParallelTracker(line, false) { itsBunch_m = &bunch; }
    };
    AbstractTimeDependence::setTimeDependence("REGISTER_RF_A",
            std::make_shared<PolynomialTimeDependence>(std::vector{0.05}));
    AbstractTimeDependence::setTimeDependence("REGISTER_RF_F",
            std::make_shared<PolynomialTimeDependence>(std::vector{3.}));
    AbstractTimeDependence::setTimeDependence("REGISTER_RF_P",
            std::make_shared<PolynomialTimeDependence>(std::vector{Physics::pi / 2.}));
    cavity.setAmplitudeName("REGISTER_RF_A");
    cavity.setFrequencyName("REGISTER_RF_F");
    cavity.setPhaseName("REGISTER_RF_P");
    cavity.setLength(0.1);
    cavity.setWidth(1.);
    cavity.setHeight(0.2);
    const auto bunch = makeBunch(1);
    const auto pc = bunch->getParticleContainer();
    FlaggedBeamline line;
    BorrowedBunchTracker tracker(line, *bunch);
    cavity.accept(tracker);
    ASSERT_EQ(observed.clones, 1u);
    ASSERT_EQ(observed.attached, bunch.get());
    ASSERT_NE(observed.runtime, nullptr);
    EXPECT_NE(observed.runtime, &cavity);
    EXPECT_DOUBLE_EQ(observed.runtime->getAmplitude(0.), 0.05);
    Vector_t<double, 3> hostE(0), hostB(0);
    EXPECT_FALSE(observed.runtime->applyToReferenceParticle({0., 0., 0.05}, {}, 0., hostE, hostB));
    EXPECT_DOUBLE_EQ(hostE[2], 50000.);
    Kokkos::deep_copy(pc->R.getView(), Vector_t<double, 3>{0., 0., 0.05});
    Kokkos::deep_copy(pc->E.getView(), Vector_t<double, 3>(0));
    bunch->setT(-1.e-9);
    bunch->setdT(2.e-9);
    observed.runtime->apply(pc);
    const auto fields = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->E.getView());
    for (unsigned d = 0; d < 3; ++d) EXPECT_DOUBLE_EQ(fields(0)[d], hostE[d]);
}
