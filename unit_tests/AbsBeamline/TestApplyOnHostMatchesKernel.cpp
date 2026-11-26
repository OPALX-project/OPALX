/**
 * \file TestApplyOnHostMatchesKernel.cpp
 * \brief For every element with a field: the bunch kernel apply(pc) and the one-point
 *        applyOnHost() give the same field at the same points.
 *
 * Each element is sampled on points inside, at the edges of and outside its field window and
 * aperture. The kernel runs at StepTime{t, dt}; applyOnHost() gets time.fieldTime(). Error
 * terms (DVOLT, DLAG, DKS) are set where an element has them, so both callers must include
 * them. Both start from the same nonzero E and B, so both must add rather than overwrite.
 */

#include "AbsBeamline/ConstantEFieldCavity.h"
#include "AbsBeamline/FieldmapElement.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/VariableRFCavity.h"
#include "Algorithms/PolynomialTimeDependence.h"
#include "BeamlineGeometry/Geometry.h"
#include "Fields/Fieldmap.h"
#include "Ippl.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"
#include "Structure/Beam.h"

#include "gtest/gtest.h"

#include "../BunchBeam.h"
#include "../TemporaryDirectory.h"

#include <cmath>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

extern Inform* gmsg;

namespace {
    using Vector = Vector_t<double, 3>;

    /// Points around the field window [zBegin, zEnd] of @p element, on and off axis.
    std::vector<Vector> samplePoints(const ElementBase& element) {
        double zBegin = 0.0, zEnd = 0.0;
        element.getFieldExtent(zBegin, zEnd);
        const double margin = 0.1 * (zEnd - zBegin);
        std::vector<Vector> points;
        constexpr int nz = 23;
        for (int iz = 0; iz < nz; ++iz) {
            const double z = zBegin - margin + (zEnd - zBegin + 2 * margin) * iz / (nz - 1);
            for (double x : {0.0, 0.004, -0.013, 0.03})
                for (double y : {0.0, -0.006, 0.02})
                    points.push_back(Vector(x, y, z));
        }
        points.push_back(Vector(0.0, 0.0, zBegin));  // window edges exactly
        points.push_back(Vector(0.0, 0.0, zEnd));
        return points;
    }

    /// A one-dimensional Astra RF map over [0, length] with Ez = sin(pi z / length).
    void writeAstraMap(const std::string& path, double length) {
        std::ofstream f(path);
        f << "AstraDynamic 8 TRUE\n100.0\n";
        constexpr int n = 41;
        for (int i = 0; i < n; ++i) {
            const double z = length * i / (n - 1);
            f << z << " " << std::sin(Physics::pi * z / length) + 0.1 << "\n";
        }
    }

    /// A G4beamline cylinder map with Bz = 0.25 T, Br = 0.05 T, z in [0.1, 0.6] m, r < 40 mm.
    void writeCylinderMap(const std::string& path) {
        std::ofstream f(path);
        f << "param normB=1. current=1.\n";
        f << "cylinder Z0=100 nZ=11 dZ=50 nR=5 dR=10\n";
        for (const double value : {0.25, 0.05}) {
            f << (value == 0.25 ? "Bz\n" : "Br\n");
            for (int iz = 0; iz < 11; ++iz)
                for (int ir = 0; ir < 5; ++ir)
                    f << value * (1.0 + 0.1 * iz + 0.05 * ir) << (ir + 1 == 5 ? '\n' : '\t');
        }
    }
}  // namespace

class ApplyOnHostMatchesKernelTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        gmsg = new Inform(nullptr, -1);
    }
    static void TearDownTestSuite() {
        delete gmsg;
        gmsg = nullptr;
        ippl::finalize();
    }
    void SetUp() override {
        directory_m.create("opalx-apply-on-host");
        Beam::defineBuiltins();
        bunch_m = makeBunch(1);
    }

    static std::shared_ptr<PartBunch_t> makeBunch(size_t count) {
        auto bunch = std::make_shared<PartBunch_t>(
                std::vector{opalx::test::bunchBeam(*Beam::find("UNNAMED_BEAM"), count)}, 1.0,
                opalx::spacecharge::CartesianDomainConfig3D{.periodicParticleBoundary = true});
        bunch->getParticleContainer()->createParticles(count);
        return bunch;
    }
    void TearDown() override {
        bunch_m.reset();
        Fieldmap::clearDictionary();
        directory_m.reset();
    }

    std::string file(const std::string& name) const { return (directory_m.path() / name).string(); }

    /// Run the kernel and applyOnHost() on samplePoints(element) and compare the fields.
    void expectSameField(ElementBase& element) {
        SCOPED_TRACE(element.getName());
        const std::vector<Vector> points = samplePoints(element);
        const auto bunch                 = makeBunch(points.size());
        auto pc                          = bunch->getParticleContainer();
        auto R                           = Kokkos::create_mirror_view(pc->R.getView());
        for (size_t i = 0; i < points.size(); ++i)
            R(i) = points[i];
        const Vector initialE(1.0, 2.0, 3.0), initialB(0.1, 0.2, 0.3);
        Kokkos::deep_copy(pc->R.getView(), R);
        Kokkos::deep_copy(pc->E.getView(), initialE);
        Kokkos::deep_copy(pc->B.getView(), initialB);

        const StepTime time{1.0e-9, 2.0e-10, 0};
        element.apply(pc, time);
        const auto E = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->E.getView());
        const auto B = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->B.getView());

        size_t withField = 0;
        for (size_t i = 0; i < points.size(); ++i) {
            SCOPED_TRACE(i);
            Vector hostE(initialE), hostB(initialB);
            element.applyOnHost(points[i], time.fieldTime(), hostE, hostB);
            bool changed = false;
            for (unsigned d = 0; d < 3; ++d)
                changed = changed || hostE(d) != initialE(d) || hostB(d) != initialB(d);
            if (changed) ++withField;
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_NEAR(E(i)(d), hostE(d), 1e-12 * std::max(1.0, std::abs(hostE(d))));
                EXPECT_NEAR(B(i)(d), hostB(d), 1e-12 * std::max(1.0, std::abs(hostB(d))));
            }
        }
        // Not vacuous: the element has a field at a good part of the points.
        EXPECT_GT(withField, points.size() / 10);
    }

    opalx::test::TemporaryDirectory directory_m;
    std::shared_ptr<PartBunch_t> bunch_m;
};

TEST_F(ApplyOnHostMatchesKernelTest, Multipole) {
    Multipole multipole("MULTIPOLE");
    multipole.getGeometry().setElementLength(1.0);
    multipole.setNormalComponent(1, 2.0);
    multipole.setNormalComponent(3, 40.0);
    multipole.setSkewComponent(2, 5.0);
    multipole.setAperture(ApertureType::RECTANGULAR, {0.05, 0.05});
    expectSameField(multipole);
}

TEST_F(ApplyOnHostMatchesKernelTest, MultipoleTStraightAndCurved) {
    for (const double angle : {0.0, Physics::pi / 8.0}) {
        MultipoleT magnet(angle == 0.0 ? "MULTIPOLET_STRAIGHT" : "MULTIPOLET_CURVED");
        magnet.setBendAngle(angle, false);
        magnet.setElementLength(1.0);
        magnet.setAperture(0.1, 0.1);
        magnet.setFringeField(0.5, 0.05, 0.05);
        magnet.setRotation(0.0);
        magnet.setEntranceAngle(0.0);
        magnet.setMaxOrder(3, 5);
        magnet.setTransProfile({0.1, 1.5, 3.0});
        magnet.initialise(bunch_m.get());
        expectSameField(magnet);
    }
}

TEST_F(ApplyOnHostMatchesKernelTest, SBendAndRBend) {
    SBend sector("SBEND");
    sector.getGeometry() = Geometry::makeSBend(1.0, 0.2);
    sector.getGeometry().setElementLength(1.0);
    sector.getGeometry().setBendAngle(0.2);
    sector.setB(-1.0);
    sector.setFullGap(0.04);
    sector.setAperture(ApertureType::RECTANGULAR, {0.05, 0.05});
    expectSameField(sector);

    RBend rectangular("RBEND");
    rectangular.getGeometry() = Geometry::makeRBend(1.0, 0.2);
    rectangular.getGeometry().setElementLength(1.0);
    rectangular.setB(-1.0);
    rectangular.setFullGap(0.04);
    rectangular.setAperture(ApertureType::RECTANGULAR, {0.05, 0.05});
    expectSameField(rectangular);
}

TEST_F(ApplyOnHostMatchesKernelTest, SolenoidWithError) {
    const std::string map = file("solenoid.g4blmap");
    writeCylinderMap(map);
    Solenoid solenoid("SOLENOID");
    solenoid.setFieldMapFN(map);
    solenoid.setKS(2.0);
    solenoid.setDKS(0.1);
    solenoid.initialise(bunch_m.get());
    expectSameField(solenoid);
}

TEST_F(ApplyOnHostMatchesKernelTest, FieldmapElement) {
    const std::string map = file("fieldmap.g4blmap");
    writeCylinderMap(map);
    FieldmapElement element("FIELDMAP");
    element.setFieldMapFN(map);
    element.setBScale(1.5);
    element.initialise(nullptr);
    Fieldmap::readMap(map);
    expectSameField(element);
}

TEST_F(ApplyOnHostMatchesKernelTest, RFCavityIdealAndMapWithErrors) {
    RFCavity ideal("RFCAVITY_IDEAL");
    ideal.setCavityType("IDEAL");
    ideal.getGeometry().setElementLength(0.2);
    ideal.setAmplitudem(5.0);
    ideal.setAmplitudeError(0.5);
    ideal.setFrequencym(Physics::two_pi * 100.e6);
    ideal.setPhasem(0.3);
    ideal.setPhaseError(0.1);
    ideal.setAperture(ApertureType::RECTANGULAR, {0.02, 0.02});
    ideal.initialise(bunch_m.get());
    expectSameField(ideal);

    const std::string map = file("cavity.astra");
    writeAstraMap(map, 0.3);
    RFCavity cavity("RFCAVITY_MAP");
    cavity.setFieldMapFN(map);
    cavity.getGeometry().setElementLength(0.3);
    cavity.setAmplitudem(10.0);
    cavity.setAmplitudeError(0.7);
    cavity.setFrequencym(Physics::two_pi * 100.e6);
    cavity.setPhasem(0.4);
    cavity.setPhaseError(-0.2);
    cavity.initialise(bunch_m.get());
    expectSameField(cavity);
}

TEST_F(ApplyOnHostMatchesKernelTest, TravelingWaveAllThreeCellsWithErrors) {
    const std::string map = file("traveling.astra");
    writeAstraMap(map, 0.1);
    TravelingWave cavity("TRAVELINGWAVE");
    cavity.setFieldMapFN(map);
    cavity.getGeometry().setElementLength(0.4);
    cavity.setNumCells(9);
    cavity.setMode(1.0 / 3.0);
    cavity.setAmplitudem(10.0);
    cavity.setAmplitudeError(0.7);
    cavity.setFrequencym(Physics::two_pi * 100.e6);
    cavity.setPhasem(0.4);
    cavity.setPhaseError(-0.2);
    cavity.initialise(bunch_m.get());
    expectSameField(cavity);
}

TEST_F(ApplyOnHostMatchesKernelTest, VariableRFCavity) {
    VariableRFCavity cavity("VARIABLE_RF_CAVITY");
    cavity.setAmplitudeModel(std::make_shared<PolynomialTimeDependence>(std::vector{1.0, 1.e8}));
    cavity.setFrequencyModel(std::make_shared<PolynomialTimeDependence>(std::vector{100.0}));
    cavity.setPhaseModel(std::make_shared<PolynomialTimeDependence>(std::vector{0.3, 1.e9}));
    cavity.setLength(0.5);
    cavity.setWidth(0.05);
    cavity.setHeight(0.03);
    cavity.initialise(bunch_m.get());
    expectSameField(cavity);
}

TEST_F(ApplyOnHostMatchesKernelTest, ConstantEFieldCavity) {
    ConstantEFieldCavity cavity("CONSTANTEFIELDCAVITY");
    cavity.getGeometry().setElementLength(0.5);
    cavity.setEx(1.0e5);
    cavity.setEy(-2.0e5);
    cavity.setEz(3.0e6);
    cavity.setAperture(ApertureType::ELLIPTICAL, {0.02, 0.01});
    expectSameField(cavity);
}
