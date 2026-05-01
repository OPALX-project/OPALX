#include "AbsBeamline/Solenoid.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/IndexMap.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/MultipoleRep.h"
#include "BeamlineCore/RFCavityRep.h"
#include "BeamlineCore/SolenoidRep.h"
#include "BeamlineCore/TravelingWaveRep.h"
#include "BeamlineGeometry/NullGeometry.h"
#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"
#include "Fields/Fieldmap.h"
#include "Physics/Units.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
#include "Structure/MeshGenerator.h"

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class SolenoidPlacementTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        gmsg                = new Inform(nullptr, -1);
        Options::enableHDF5 = false;
        std::filesystem::create_directories("data");
    }

    static void TearDownTestSuite() {
        delete gmsg;
        gmsg = nullptr;
        ippl::finalize();
    }

    void SetUp() override {
        const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        testStem_        = std::string("solenoid_") + info->name();
        OpalData::getInstance()->storeInputFn(testStem_ + ".opal");
        OpalData::getInstance()->setOpenMode(OpalData::OpenMode::WRITE);
        cleanupOutputs();
    }

    void TearDown() override {
        Fieldmap::freeMap(fieldmapFile_.string());
        cleanupOutputs();
    }

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
        EXPECT_NE(opBeam, nullptr);

        auto bunch = std::make_shared<PartBunch_t>(
                /*qi=*/std::vector{1.0}, /*mi=*/std::vector{1.0},
                /*beams=*/std::vector<Beam*>{opBeam},
                /*totalParticlesPerBeam=*/std::vector<size_t>{numParticles},
                /*lbt=*/1.0, /*integration_method=*/"LF2", fsCmdBase_m.get(), dataSink_m.get());
        bunch->getParticleContainer()->createParticles(numParticles);
        return bunch;
    }

    std::filesystem::path writeXZFieldmap(
            const std::string& filename, double zBeginCm, double zEndCm, int nz, double rBeginCm,
            double rEndCm, int nr, double bzValue = 1.0, double brValue = 0.0) {
        fieldmapFile_ = std::filesystem::path("data") / filename;
        std::ofstream out(fieldmapFile_);
        out << "2DMagnetoStatic XZ\n";
        out << zBeginCm << " " << zEndCm << " " << nz - 1 << "\n";
        out << rBeginCm << " " << rEndCm << " " << nr - 1 << "\n";
        for (int ir = 0; ir < nr; ++ir) {
            for (int iz = 0; iz < nz; ++iz) {
                out << bzValue << " " << brValue << "\n";
            }
        }
        return fieldmapFile_;
    }

    std::filesystem::path outputPath(const std::string& suffix) const {
        return std::filesystem::path("data") / (testStem_ + suffix);
    }

    void cleanupOutputs() const {
        std::filesystem::remove(outputPath("_ElementPositions.txt"));
        std::filesystem::remove(outputPath("_ElementPositions.py"));
        std::filesystem::remove(outputPath("_ElementPositions.sdds"));
        std::filesystem::remove(outputPath("_ElementPositions.vtk"));
        std::filesystem::remove(outputPath("_ElementPositions.html"));
        std::filesystem::remove(outputPath("_ElementPositions.gpl"));
        if (!fieldmapFile_.empty()) {
            std::filesystem::remove(fieldmapFile_);
        }
    }

    static std::tuple<double, double, double> parsePositionLine(const std::string& line) {
        const auto quote = line.rfind('"');
        EXPECT_NE(quote, std::string::npos);
        std::istringstream values(line.substr(quote + 1));
        double z = 0.0, x = 0.0, y = 0.0;
        values >> z >> x >> y;
        return {z, x, y};
    }

    static std::vector<std::string> splitWhitespace(const std::string& line) {
        std::vector<std::string> tokens;
        std::istringstream in(line);
        std::string token;
        while (in >> token) {
            tokens.push_back(token);
        }
        return tokens;
    }

    std::string testStem_;
    std::filesystem::path fieldmapFile_;
    std::shared_ptr<FieldSolverCmd> fsCmdBase_m;
    std::shared_ptr<DataSink> dataSink_m;
};

class DummyBeamline final : public Beamline {
public:
    DummyBeamline() : Beamline("dummy") {}

    ElementType getType() const override { return ElementType::BEAMLINE; }
    BGeometryBase& getGeometry() override { return geometry_; }
    const BGeometryBase& getGeometry() const override { return geometry_; }
    void accept(BeamlineVisitor& visitor) const override { visitor.visitBeamline(*this); }
    ElementBase* clone() const override { return new DummyBeamline(*this); }
    void iterate(BeamlineVisitor&, bool) const override {}

private:
    NullGeometry geometry_;
};

TEST_F(SolenoidPlacementTest, FieldMapEdgesAndSupportEnvelopeFollowFieldMap) {
    const auto mapFile = writeXZFieldmap("solenoid_edges.map", -5.0, 15.0, 4, 0.0, 3.0, 3);

    SolenoidRep solenoid("SOL1");
    solenoid.setElementLength(1.0);
    solenoid.setFieldMapFN(mapFile.string());
    solenoid.setElementPosition(1.0);

    double startField = 1.0;
    double endField   = 0.0;
    try {
        solenoid.initialise(nullptr, startField, endField);
    } catch (const OpalException& ex) {
        FAIL() << ex.where() << ": " << ex.what();
    } catch (...) {
        FAIL() << "non-OPAL exception";
    }

    double fieldBegin = 0.0, fieldEnd = 0.0;
    solenoid.getFieldExtend(fieldBegin, fieldEnd);

    double bodyBegin = 0.0, bodyEnd = 0.0;
    solenoid.getElementDimensions(bodyBegin, bodyEnd);

    EXPECT_NEAR(startField, 0.95, 1e-12);
    EXPECT_NEAR(endField, 1.15, 1e-12);
    EXPECT_NEAR(solenoid.getElementLength(), 1.0, 1e-12);
    EXPECT_NEAR(fieldBegin, -0.05, 1e-12);
    EXPECT_NEAR(fieldEnd, 0.15, 1e-12);
    EXPECT_NEAR(bodyBegin, 0.0, 1e-12);
    EXPECT_NEAR(bodyEnd, 1.0, 1e-12);
    EXPECT_NEAR(solenoid.getEdgeToBegin().getOrigin()(2), 0.0, 1e-12);
    EXPECT_NEAR(solenoid.getEdgeToEnd().getOrigin()(2), 1.0, 1e-12);
    EXPECT_TRUE(solenoid.isInside({0.0, 0.0, 0.00}));
    EXPECT_TRUE(solenoid.isInside({0.0, 0.0, 0.14}));
    EXPECT_FALSE(solenoid.isInside({0.0, 0.0, 0.20}));

    double horizontalRadius = 0.0;
    double verticalRadius   = 0.0;
    ASSERT_NO_THROW({
        const bool hasSupport = solenoid.getSupportEnvelope(horizontalRadius, verticalRadius);
        ASSERT_TRUE(hasSupport);
    });
    EXPECT_NEAR(horizontalRadius, 0.03, 1e-12);
    EXPECT_NEAR(verticalRadius, 0.03, 1e-12);
}

TEST_F(SolenoidPlacementTest, LatticeExportsUseFieldMapEdgesAndSolenoidMeshType) {
    const auto mapFile = writeXZFieldmap("solenoid_lattice.map", -5.0, 15.0, 4, 0.0, 3.0, 3);

    SolenoidRep solenoid("SOL1");
    solenoid.setElementLength(1.0);
    solenoid.setFieldMapFN(mapFile.string());
    solenoid.setElementPosition(1.0);

    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);
    OpalBeamline beamline;
    try {
        beamline.visit(solenoid, visitor, *bunch);
    } catch (const OpalException& ex) {
        FAIL() << ex.where() << ": " << ex.what();
    } catch (...) {
        FAIL() << "non-OPAL exception";
    }
    ASSERT_NO_THROW(beamline.compute3DLattice());

    const auto elements = beamline.getElements();
    ASSERT_EQ(elements.size(), 1u);
    const auto& placedComponent = *elements.begin();
    const PlacedElement placed  = beamline.getPlacedElement(placedComponent);
    EXPECT_NEAR(placed.getNominalEntryTransform().getOrigin()(2), 1.0, 1e-12);
    EXPECT_NEAR(placed.getNominalExitTransform().getOrigin()(2), 2.0, 1e-12);
    EXPECT_EQ(beamline.getElements(Vector_t<double, 3>(0.0, 0.0, 1.10)).size(), 1u);
    EXPECT_TRUE(beamline.getElements(Vector_t<double, 3>(0.0, 0.0, 1.30)).empty());

    ASSERT_NO_THROW(beamline.save3DLattice());

    std::ifstream txt(outputPath("_ElementPositions.txt"));
    ASSERT_TRUE(txt.is_open());

    std::string textLine;
    bool foundBegin = false;
    bool foundEnd   = false;
    while (std::getline(txt, textLine)) {
        if (textLine.find("\"BEGIN: SOL1\"") != std::string::npos) {
            const auto [z, x, y] = parsePositionLine(textLine);
            const auto entry     = placed.getNominalEntryTransform().getOrigin();
            EXPECT_NEAR(z, entry(2), 1e-12);
            EXPECT_NEAR(x, entry(0), 1e-12);
            EXPECT_NEAR(y, entry(1), 1e-12);
            foundBegin = true;
        } else if (textLine.find("\"END: SOL1\"") != std::string::npos) {
            const auto [z, x, y] = parsePositionLine(textLine);
            const auto exit      = placed.getNominalExitTransform().getOrigin();
            EXPECT_NEAR(z, exit(2), 1e-12);
            EXPECT_NEAR(x, exit(0), 1e-12);
            EXPECT_NEAR(y, exit(1), 1e-12);
            foundEnd = true;
        }
    }
    EXPECT_TRUE(foundBegin);
    EXPECT_TRUE(foundEnd);

    std::ifstream py(outputPath("_ElementPositions.py"));
    ASSERT_TRUE(py.is_open());
    const std::string script(
            (std::istreambuf_iterator<char>(py)), std::istreambuf_iterator<char>());
    EXPECT_NE(script.find("color = [5]"), std::string::npos);
    EXPECT_NE(script.find("numVertices = [432]"), std::string::npos);
    EXPECT_NE(
            script.find("os.path.getmtime(script_file) > os.path.getmtime(vtk_file)"),
            std::string::npos);
}

TEST_F(SolenoidPlacementTest, ElementPositionsSDDSMarksSolenoidColumn) {
    auto solenoid = std::make_shared<SolenoidRep>("SOL1");

    IndexMap map;
    IndexMap::value_t elements;
    elements.insert(solenoid);
    map.add(1.0, 2.0, elements);
    map.saveSDDS(0.0);

    std::ifstream sdds(outputPath("_ElementPositions.sdds"));
    ASSERT_TRUE(sdds.is_open());

    std::string line;
    bool foundSolenoidRow = false;
    while (std::getline(sdds, line)) {
        if (line.find("\"SOL1\"") == std::string::npos) {
            continue;
        }

        const auto tokens = splitWhitespace(line);
        ASSERT_GE(tokens.size(), 12u);
        EXPECT_EQ(tokens[7], "1");
        EXPECT_EQ(tokens[10], "0");
        foundSolenoidRow = true;
    }
    EXPECT_TRUE(foundSolenoidRow);
}

TEST_F(SolenoidPlacementTest, DriftMeshesAsBlueCylinderUsingFirstNonDriftReference) {
    const auto mapFile = writeXZFieldmap("solenoid_drift_mesh.map", -5.0, 15.0, 4, 0.0, 3.0, 3);

    DriftRep drift("DR1");
    drift.setElementLength(0.5);
    drift.setElementPosition(0.0);

    SolenoidRep solenoid("SOL1");
    solenoid.setElementLength(1.0);
    solenoid.setFieldMapFN(mapFile.string());
    solenoid.setElementPosition(0.5);

    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);
    OpalBeamline beamline;
    beamline.visit(drift, visitor, *bunch);
    beamline.visit(solenoid, visitor, *bunch);
    beamline.compute3DLattice();
    beamline.save3DLattice();

    std::ifstream py(outputPath("_ElementPositions.py"));
    ASSERT_TRUE(py.is_open());
    const std::string script(
            (std::istreambuf_iterator<char>(py)), std::istreambuf_iterator<char>());
    EXPECT_NE(script.find("color = [8, 5]"), std::string::npos);
    EXPECT_NE(script.find("numVertices = [144, 432]"), std::string::npos);
}

TEST_F(SolenoidPlacementTest, QuadrupoleMeshesAsPoleBodyRatherThanGenericCylinder) {
    MultipoleRep quadrupole("Q1");
    quadrupole.setElementLength(0.4);
    quadrupole.setElementPosition(0.0);
    quadrupole.setAperture(ApertureType::ELLIPTICAL, std::vector<double>{0.02, 0.03, 1.0});
    quadrupole.setNormalComponent(1, 1.0);

    auto bunch = makeBunch(0);
    DummyBeamline beamlineForVisitor;
    DefaultVisitor visitor(beamlineForVisitor, false, false);
    OpalBeamline beamline;
    beamline.visit(quadrupole, visitor, *bunch);
    beamline.compute3DLattice();
    beamline.save3DLattice();

    std::ifstream py(outputPath("_ElementPositions.py"));
    ASSERT_TRUE(py.is_open());
    const std::string script(
            (std::istreambuf_iterator<char>(py)), std::istreambuf_iterator<char>());
    EXPECT_NE(script.find("color = [2]"), std::string::npos);
    EXPECT_NE(script.find("numVertices = [32]"), std::string::npos);
}

TEST_F(SolenoidPlacementTest, RFCavityMeshesAsBulgedCellStructure) {
    RFCavityRep cavity("RFC1");
    cavity.setElementLength(0.7);
    cavity.setElementPosition(0.0);
    cavity.setAperture(ApertureType::ELLIPTICAL, std::vector<double>{0.02, 0.03, 1.0});

    MeshGenerator mesh;
    mesh.add(cavity);
    mesh.write(testStem_);

    std::ifstream py(outputPath("_ElementPositions.py"));
    ASSERT_TRUE(py.is_open());
    const std::string script(
            (std::istreambuf_iterator<char>(py)), std::istreambuf_iterator<char>());
    EXPECT_NE(script.find("color = [6]"), std::string::npos);
    EXPECT_NE(script.find("numVertices = [1008]"), std::string::npos);
}

TEST_F(SolenoidPlacementTest, TravelingWaveMeshesAsPeriodicStructure) {
    TravelingWaveRep travelingWave("TW1");
    travelingWave.setElementLength(0.8);
    travelingWave.setElementPosition(0.0);
    travelingWave.setAperture(ApertureType::ELLIPTICAL, std::vector<double>{0.02, 0.03, 1.0});

    MeshGenerator mesh;
    mesh.add(travelingWave);
    mesh.write(testStem_);

    std::ifstream py(outputPath("_ElementPositions.py"));
    ASSERT_TRUE(py.is_open());
    const std::string script(
            (std::istreambuf_iterator<char>(py)), std::istreambuf_iterator<char>());
    EXPECT_NE(script.find("color = [7]"), std::string::npos);
    EXPECT_NE(script.find("numVertices = [1152]"), std::string::npos);
}
