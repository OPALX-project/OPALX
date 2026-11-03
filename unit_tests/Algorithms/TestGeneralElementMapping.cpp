/**
 * \file TestGeneralElementMapping.cpp
 * \brief Unit tests for the 3D grid lookup of active elements and its reference-orbit checks.
 *
 * A grid with a single cell holds every element in that cell, so its answers are those of
 * the exact box test alone; finer grids must return the same sets. Point sampling of the
 * bunch box gives an independent check that no true overlap is missed.
 */

#include <gtest/gtest.h>

#include <unistd.h>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <random>
#include <set>
#include <string>

#include "Ippl.h"

#include "AbsBeamline/CyclotronSector.h"
#include "Algorithms/GeneralElementMapping.h"
#include "Attributes/Attributes.h"
#include "BasicActions/Option.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/MultipoleRep.h"
#include "BeamlineCore/MultipoleTRep.h"
#include "BeamlineCore/SBendRep.h"
#include "Physics/Physics.h"
#include "Utilities/ElementList.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <limits>
#include <vector>

extern Inform* gmsg;

namespace {
    using Vector = Vector_t<double, 3>;

    constexpr double reach = 0.5;

    /// Rotation about y by theta, as OpalElement composes THETA.
    Quaternion rotationAboutY(double theta) {
        return Quaternion(std::cos(0.5 * theta), 0, std::sin(0.5 * theta), 0);
    }

    Quaternion rotation(double theta, double phi, double psi) {
        const Quaternion y(std::cos(0.5 * theta), 0, std::sin(0.5 * theta), 0);
        const Quaternion x(std::cos(0.5 * phi), std::sin(0.5 * phi), 0, 0);
        const Quaternion z(std::cos(0.5 * psi), 0, 0, std::sin(0.5 * psi));
        return y * (x * z);
    }

    std::shared_ptr<DriftRep> makeDrift(
            const std::string& name, double length, const Vector& origin,
            const Quaternion& orientation = Quaternion(1, 0, 0, 0)) {
        auto drift = std::make_shared<DriftRep>(name);
        drift->getGeometry().setElementLength(length);
        drift->setCSTrafoGlobal2Local(CoordinateSystemTrafo(origin, orientation.conjugate()));
        return drift;
    }

    /// The bunch-to-lab transform of a bunch frame at @p origin, turned by @p orientation.
    CoordinateSystemTrafo bunchToLab(
            const Vector& origin, const Quaternion& orientation = Quaternion(1, 0, 0, 0)) {
        return CoordinateSystemTrafo(origin, orientation.conjugate()).inverted();
    }

    std::set<std::string> names(const GeneralElementMapping::ElementSet& elements) {
        std::set<std::string> result;
        for (const auto& element : elements)
            result.insert(element->getName());
        return result;
    }
}  // namespace

class GeneralElementMappingTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        if (gmsg == nullptr) gmsg = new Inform(nullptr, -1);
    }

    static void TearDownTestSuite() {
        delete gmsg;
        gmsg = nullptr;
        ippl::finalize();
    }
};

TEST_F(GeneralElementMappingTest, EmptyLatticeSelectsNothing) {
    const GeneralElementMapping mapping(ElementList{}, {8, 8, 8}, reach);
    EXPECT_TRUE(mapping.query(Vector(-1.0), Vector(1.0), bunchToLab(Vector(0.0))).empty());
}

TEST_F(GeneralElementMappingTest, DriftIsActiveOnlyWhereTheBunchBoxReachesIt) {
    const ElementList lattice = {makeDrift("D", 1.0, Vector(0.0))};
    const GeneralElementMapping mapping(lattice, {8, 8, 8}, reach);
    const auto identity = bunchToLab(Vector(0.0));

    EXPECT_EQ(
            names(mapping.query(Vector({-0.1, -0.1, 0.4}), Vector({0.1, 0.1, 0.6}), identity)),
            std::set<std::string>{"D"});
    EXPECT_TRUE(
            mapping.query(Vector({-0.1, -0.1, 2.0}), Vector({0.1, 0.1, 3.0}), identity).empty());
    EXPECT_TRUE(mapping.query(Vector({0.6, -0.1, 0.4}), Vector({0.7, 0.1, 0.6}), identity).empty());
    EXPECT_EQ(
            names(mapping.query(Vector({0.45, -0.1, 0.4}), Vector({0.7, 0.1, 0.6}), identity)),
            std::set<std::string>{"D"});
    // A zero-size box at a point inside the element.
    EXPECT_EQ(
            names(mapping.query(Vector({0.0, 0.0, 0.5}), Vector({0.0, 0.0, 0.5}), identity)),
            std::set<std::string>{"D"});
}

TEST_F(GeneralElementMappingTest, BunchFrameIsTransformedToTheLab) {
    // Bunch frame placed at lab z = 10 m: a bunch-frame box around its origin reaches a drift
    // at lab z = 10 m but not one at the lab origin.
    const ElementList lattice = {
            makeDrift("NEAR_ORIGIN", 1.0, Vector(0.0)),
            makeDrift("FAR", 1.0, Vector({0.0, 0.0, 10.0}))};
    const GeneralElementMapping mapping(lattice, {8, 8, 8}, reach);

    EXPECT_EQ(
            names(mapping.query(Vector(-0.1), Vector(0.1), bunchToLab(Vector({0, 0, 10.2})))),
            std::set<std::string>{"FAR"});
}

TEST_F(GeneralElementMappingTest, ParallelBranchesAreSeparated) {
    const ElementList lattice = {
            makeDrift("A", 2.0, Vector(0.0)), makeDrift("B", 2.0, Vector({2.0, 0.0, 0.0}))};
    const GeneralElementMapping mapping(lattice, {8, 8, 8}, reach);

    EXPECT_EQ(
            names(mapping.query(Vector(-0.2), Vector(0.2), bunchToLab(Vector({0, 0, 1})))),
            std::set<std::string>{"A"});
    EXPECT_EQ(
            names(mapping.query(Vector(-0.2), Vector(0.2), bunchToLab(Vector({2, 0, 1})))),
            std::set<std::string>{"B"});
}

TEST_F(GeneralElementMappingTest, RotatedElementIsTestedAsAnOrientedBox) {
    // A 4 m drift turned by 45 degrees: its lab bounding box contains the point (2.5, 0, 0.5),
    // but the oriented element box does not.
    const ElementList lattice = {
            makeDrift("TURNED", 4.0, Vector(0.0), rotationAboutY(0.25 * Physics::pi))};
    const GeneralElementMapping mapping(lattice, {8, 8, 8}, reach);

    EXPECT_TRUE(mapping.query(Vector(-0.05), Vector(0.05), bunchToLab(Vector({2.5, 0.0, 0.5})))
                        .empty());
    EXPECT_EQ(
            names(mapping.query(Vector(-0.05), Vector(0.05), bunchToLab(Vector({1.4, 0.0, 1.4})))),
            std::set<std::string>{"TURNED"});
    // A bunch frame turned with the element, with its long axis along the element.
    EXPECT_EQ(
            names(mapping.query(
                    Vector({-0.1, -0.1, -0.5}), Vector({0.1, 0.1, 0.5}),
                    bunchToLab(Vector({1.4, 0.0, 1.4}), rotationAboutY(0.25 * Physics::pi)))),
            std::set<std::string>{"TURNED"});
}

TEST_F(GeneralElementMappingTest, ElementReachOverridesTheGlobalReach) {
    auto narrow = makeDrift("NARROW", 1.0, Vector(0.0));
    narrow->setReach(0.1);
    const ElementList lattice = {narrow, makeDrift("WIDE", 1.0, Vector({0.0, 0.0, 2.0}))};
    const GeneralElementMapping mapping(lattice, {8, 8, 8}, reach);
    const auto identity = bunchToLab(Vector(0.0));

    EXPECT_TRUE(mapping.query(Vector({0.2, 0, 0.4}), Vector({0.3, 0, 0.6}), identity).empty());
    EXPECT_EQ(
            names(mapping.query(Vector({0.2, 0, 2.4}), Vector({0.3, 0, 2.6}), identity)),
            std::set<std::string>{"WIDE"});
}

TEST_F(GeneralElementMappingTest, BunchOutsideAFiniteApertureIsActiveWithinReach) {
    auto collimator = makeDrift("COLL", 1.0, Vector(0.0));
    collimator->setAperture(ApertureType::RECTANGULAR, {0.01, 0.01});
    const GeneralElementMapping mapping(ElementList{collimator}, {8, 8, 8}, reach);
    const auto identity = bunchToLab(Vector(0.0));

    EXPECT_EQ(
            names(mapping.query(Vector({0.2, -0.01, 0.4}), Vector({0.3, 0.01, 0.6}), identity)),
            std::set<std::string>{"COLL"});
    EXPECT_TRUE(
            mapping.query(Vector({0.6, -0.01, 0.4}), Vector({0.7, 0.01, 0.6}), identity).empty());
}

TEST_F(GeneralElementMappingTest, MisalignmentMovesTheElementBox) {
    auto shifted = makeDrift("SHIFTED", 1.0, Vector(0.0));
    shifted->setMisalignment(
            CoordinateSystemTrafo(Vector({1.0, 0.0, 0.0}), Quaternion(1, 0, 0, 0)));
    const GeneralElementMapping mapping(ElementList{shifted}, {8, 8, 8}, reach);
    const auto identity = bunchToLab(Vector(0.0));

    EXPECT_EQ(
            names(mapping.query(Vector({0.9, -0.1, 0.4}), Vector({1.1, 0.1, 0.6}), identity)),
            std::set<std::string>{"SHIFTED"});
    EXPECT_TRUE(
            mapping.query(Vector({-0.1, -0.1, 0.4}), Vector({0.1, 0.1, 0.6}), identity).empty());
}

TEST_F(GeneralElementMappingTest, GridAnswersEqualTheSingleCellAnswersAndMissNoOverlap) {
    std::mt19937 random(20261103);
    std::uniform_real_distribution<double> position(-10.0, 10.0);
    std::uniform_real_distribution<double> angle(-Physics::pi, Physics::pi);
    std::uniform_real_distribution<double> length(0.1, 3.0);
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    ElementList lattice;
    for (int i = 0; i < 40; ++i) {
        auto drift = makeDrift(
                "E" + std::to_string(i), length(random),
                Vector({position(random), position(random), position(random)}),
                rotation(angle(random), angle(random), angle(random)));
        if (unit(random) < 0.5)
            drift->setAperture(ApertureType::RECTANGULAR, {length(random), 0.2 * length(random)});
        lattice.push_back(drift);
    }

    const GeneralElementMapping oracle(lattice, {1, 1, 1}, reach);
    const std::array<std::array<unsigned, 3>, 3> grids = {{{8, 8, 8}, {32, 4, 2}, {3, 17, 5}}};
    std::vector<GeneralElementMapping> mappings;
    for (const auto& grid : grids)
        mappings.emplace_back(lattice, grid, reach);

    int nonEmpty = 0;
    for (int q = 0; q < 300; ++q) {
        const Vector half({0.5 * length(random), 0.5 * length(random), length(random)});
        const auto toLab = bunchToLab(
                Vector({position(random), position(random), position(random)}),
                rotation(angle(random), angle(random), angle(random)));
        const auto expected = oracle.query(-half, half, toLab);
        nonEmpty += expected.empty() ? 0 : 1;
        for (std::size_t g = 0; g < mappings.size(); ++g)
            EXPECT_EQ(names(mappings[g].query(-half, half, toLab)), names(expected))
                    << "query " << q << " grid " << g;

        // Every element that contains a sampled point of the bunch box must be selected.
        for (int i = 0; i <= 4; ++i) {
            for (int j = 0; j <= 4; ++j) {
                for (int k = 0; k <= 4; ++k) {
                    const Vector local(
                            {half(0) * (0.5 * i - 1.0), half(1) * (0.5 * j - 1.0),
                             half(2) * (0.5 * k - 1.0)});
                    const Vector lab = toLab.transformTo(local);
                    for (const auto& element : lattice) {
                        const Vector p =
                                (element->getMisalignment() * element->getCSTrafoGlobal2Local())
                                        .transformTo(lab);
                        const auto [lower, upper] = element->getLocalSupportBox(reach).getCorners();
                        const bool inside = p(0) > lower(0) && p(0) < upper(0) && p(1) > lower(1)
                                            && p(1) < upper(1) && p(2) > lower(2)
                                            && p(2) < upper(2);
                        if (inside) EXPECT_EQ(expected.count(element), 1u) << element->getName();
                    }
                }
            }
        }
    }
    EXPECT_GT(nonEmpty, 20);  // The random queries must exercise overlaps.
}

TEST_F(GeneralElementMappingTest, ReferenceFarOffAnUnaperturedElementIsReported) {
    const ElementList lattice = {makeDrift("D", 1.0, Vector(0.0))};
    const GeneralElementMapping mapping(lattice, {8, 8, 8}, reach);

    EXPECT_TRUE(mapping.checkReferencePoint(Vector({0.3, 0.0, 0.5}), 0.5).empty());
    const auto findings = mapping.checkReferencePoint(Vector({0.8, 0.0, 0.5}), 0.5);
    ASSERT_EQ(findings.size(), 1u);
    EXPECT_EQ(findings.front().element->getName(), "D");
    EXPECT_NEAR(findings.front().offset, 0.8, 1e-12);
}

TEST_F(GeneralElementMappingTest, ReferenceThatWouldBeDeletedStopsTheRun) {
    auto quad = std::make_shared<MultipoleRep>("QB");
    quad->getGeometry().setElementLength(1.0);
    quad->setAperture(ApertureType::RECTANGULAR, {0.05, 0.05});
    quad->setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector(0.0), Quaternion(1, 0, 0, 0)));

    {
        const GeneralElementMapping mapping(ElementList{quad}, {8, 8, 8}, reach);
        EXPECT_THROW(mapping.checkReferencePoint(Vector({0.3, 0.0, 0.5}), 12.0), OpalException);
        // On the axis, or outside the body, nothing is deleted.
        EXPECT_NO_THROW(mapping.checkReferencePoint(Vector({0.0, 0.0, 0.5}), 12.0));
        EXPECT_NO_THROW(mapping.checkReferencePoint(Vector({0.3, 0.0, 1.5}), 12.0));
    }
    {
        quad->setReach(0.2);  // the other branch is now outside this quad's box
        const GeneralElementMapping mapping(ElementList{quad}, {8, 8, 8}, reach);
        EXPECT_NO_THROW(mapping.checkReferencePoint(Vector({0.3, 0.0, 0.5}), 12.0));
    }
    {
        quad->setReach(reach);
        quad->setFlagDeleteOnTransverseExit(false);
        const GeneralElementMapping mapping(ElementList{quad}, {8, 8, 8}, reach);
        EXPECT_NO_THROW(mapping.checkReferencePoint(Vector({0.3, 0.0, 0.5}), 12.0));
    }
}

TEST_F(GeneralElementMappingTest, StraightLatticeOnAxisPassesTheChecks) {
    ElementList lattice;
    for (int cell = 0; cell < 4; ++cell) {
        auto quad = std::make_shared<MultipoleRep>("Q" + std::to_string(cell));
        quad->getGeometry().setElementLength(0.2);
        quad->setAperture(ApertureType::ELLIPTICAL, {0.03, 0.03});
        quad->setCSTrafoGlobal2Local(
                CoordinateSystemTrafo(Vector({0.0, 0.0, 1.0 * cell}), Quaternion(1, 0, 0, 0)));
        lattice.push_back(quad);
        lattice.push_back(
                makeDrift("D" + std::to_string(cell), 0.8, Vector({0.0, 0.0, cell + 0.2})));
    }
    const GeneralElementMapping mapping(lattice, {8, 8, 8}, reach);
    for (double z = -0.5; z < 4.5; z += 0.01) {
        EXPECT_NO_THROW(
                EXPECT_TRUE(mapping.checkReferencePoint(Vector({0.001, -0.002, z}), z).empty()));
    }
}

TEST_F(GeneralElementMappingTest, MappingOptionsAreValidatedAndPersist) {
    Option exemplar;
    std::unique_ptr<Option> command(exemplar.clone("MAPPING_OPTIONS"));
    EXPECT_EQ(
            Attributes::getRealArray(*command->findAttribute("ELEMENTMAPPINGGRID")),
            (std::vector<double>{8, 8, 8}));
    EXPECT_DOUBLE_EQ(Attributes::getReal(*command->findAttribute("ELEMENTMAPPINGREACH")), 0.5);

    Attributes::setRealArray(*command->findAttribute("ELEMENTMAPPINGGRID"), {4, 16, 32});
    Attributes::setReal(*command->findAttribute("ELEMENTMAPPINGREACH"), 0.25);
    command->execute();
    const std::array<unsigned, 3> grid{4, 16, 32};
    EXPECT_EQ(Options::elementMappingGrid, grid);
    EXPECT_DOUBLE_EQ(Options::elementMappingReach, 0.25);

    std::unique_ptr<Option> next(exemplar.clone("NEXT_MAPPING_OPTIONS"));
    next->execute();  // A later OPTION without these attributes keeps the values.
    EXPECT_EQ(Options::elementMappingGrid, grid);
    EXPECT_DOUBLE_EQ(Options::elementMappingReach, 0.25);

    for (const auto& invalid :
         {std::vector<double>{}, std::vector<double>{8, 8}, std::vector<double>{8, 8, 8, 8},
          std::vector<double>{0, 8, 8}, std::vector<double>{8, 1.5, 8},
          std::vector<double>{8, -1, 8}, std::vector<double>{4096, 4096, 4096}}) {
        Attributes::setRealArray(*next->findAttribute("ELEMENTMAPPINGGRID"), invalid);
        EXPECT_THROW(next->execute(), OpalException);
        EXPECT_EQ(Options::elementMappingGrid, grid);
    }
    Attributes::setRealArray(*next->findAttribute("ELEMENTMAPPINGGRID"), {4, 16, 32});
    for (double invalid :
         {0.0, -1.0, std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::quiet_NaN()}) {
        Attributes::setReal(*next->findAttribute("ELEMENTMAPPINGREACH"), invalid);
        EXPECT_THROW(next->execute(), OpalException);
        EXPECT_DOUBLE_EQ(Options::elementMappingReach, 0.25);
    }

    Options::elementMappingGrid  = {8, 8, 8};
    Options::elementMappingReach = 0.5;
}

namespace {
    /// Branch B: a 90-degree sector bend of radius 1 m with a 0.02 m half-aperture, entrance at
    /// the lab origin. Branch A continues straight along the bend's entrance tangent.
    std::shared_ptr<SBendRep> makeBranchBend() {
        auto bend           = std::make_shared<SBendRep>("B_BEND");
        bend->getGeometry() = Geometry::makeSBend(0.5 * Physics::pi, 1.0);
        bend->setAperture(ApertureType::RECTANGULAR, {0.02, 0.02});
        bend->setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector(0.0), Quaternion(1, 0, 0, 0)));
        return bend;
    }
}  // namespace

TEST_F(GeneralElementMappingTest, StraightBranchAlongABendTangentStaysOutsideTheArcBoxes) {
    auto bend = makeBranchBend();
    bend->setReach(0.05);
    const GeneralElementMapping mapping(ElementList{bend}, {8, 8, 8}, reach);
    const auto identity = bunchToLab(Vector(0.0));

    // (0, 0, 0.5) lies 0.118 m off the arc: inside the arc's bounding box, outside its band.
    EXPECT_TRUE(mapping.query(Vector({-0.01, -0.01, 0.45}), Vector({0.01, 0.01, 0.55}), identity)
                        .empty());
    EXPECT_NO_THROW(mapping.checkReferencePoint(Vector({0.0, 0.0, 0.5}), 0.5));

    // On the arc at s = 0.5 m the bend is active.
    const Vector onArc({std::cos(0.5) - 1.0, 0.0, std::sin(0.5)});
    EXPECT_EQ(
            names(mapping.query(onArc - Vector(0.01), onArc + Vector(0.01), identity)),
            std::set<std::string>{"B_BEND"});
}

TEST_F(GeneralElementMappingTest, ReachErrorSaysWhetherReachCanClearTheOrbit) {
    // Within the default reach of the bend: lowering REACH clears it.
    {
        const GeneralElementMapping mapping(ElementList{makeBranchBend()}, {8, 8, 8}, reach);
        try {
            mapping.checkReferencePoint(Vector({0.0, 0.0, 0.5}), 0.5);
            FAIL() << "expected the reach check to stop";
        } catch (const OpalException& e) {
            EXPECT_NE(e.what().find("set its REACH below"), std::string::npos) << e.what();
        }
    }
    // A MULTIPOLET whose field region (HAPERT = 2 m) covers the orbit outside its APERTURE:
    // no REACH can clear it.
    auto magnet = std::make_shared<MultipoleTRep>("MT");
    magnet->setBendAngle(0.0, false);
    magnet->setElementLength(1.0);
    magnet->setAperture(0.4, 2.0);
    magnet->setFringeField(0.5, 0.05, 0.05);
    magnet->setMaxOrder(5, 10);
    magnet->setTransProfile({1});
    magnet->ElementBase::setAperture(ApertureType::RECTANGULAR, {0.05, 0.05});
    magnet->setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector(0.0), Quaternion(1, 0, 0, 0)));
    const GeneralElementMapping mapping(ElementList{magnet}, {8, 8, 8}, reach);
    try {
        mapping.checkReferencePoint(Vector({0.7, 0.0, 0.5}), 0.5);
        FAIL() << "expected the reach check to stop";
    } catch (const OpalException& e) {
        EXPECT_NE(e.what().find("REACH cannot exclude"), std::string::npos) << e.what();
    }
}

TEST_F(GeneralElementMappingTest, ReportsWhetherAPointLiesWithinAnElementsTransverseReach) {
    const auto near = makeDrift("NEAR", 1.0, Vector(0.0));
    const auto far  = makeDrift("FAR", 1.0, Vector({0.6, 0.0, 0.0}));
    const GeneralElementMapping mapping(ElementList{near, far}, {8, 8, 8}, reach);

    // Longitudinal position does not matter, only the transverse reach.
    EXPECT_TRUE(mapping.coversTransversely(*near, Vector({0.1, 0.0, 5.0})));
    EXPECT_FALSE(mapping.coversTransversely(*far, Vector({0.0, 0.0, 0.5})));
    EXPECT_TRUE(mapping.coversTransversely(*far, Vector({0.2, 0.0, 0.5})));
}

TEST_F(GeneralElementMappingTest, SectorElementsHaveNoTransverseReachToMiss) {
    // Synthetic PSI layout from TestCyclotronSector: r = 1.0 .. 1.4 m, one 45-degree sector.
    const std::string file = (std::filesystem::temp_directory_path()
                              / ("opalx-mapping-sector-" + std::to_string(getpid()) + ".map"))
                                     .string();
    {
        std::ofstream out(file);
        out << "1000 100 0 5.625\n";
        for (int i = 0; i < 13; ++i)
            out << "header ";
        out << "5 8 a b c d e 2 a b c d 0.100000000E+01-0.100000000E+01 a b c d e f LREC= a b "
               "c d e\n";
        for (int i = 0; i < 5; ++i) {
            if (i) out << "a b c d e f\n";
            for (int channel = 0; channel < 4; ++channel)
                for (int j = 0; j < 8; ++j)
                    out << (channel == 0 ? 10 * (1 + 0.1 * i) : 0) << ' ';
            out << '\n';
        }
    }
    auto sector = std::make_shared<CyclotronSector>("SECTOR");
    sector->configure(CyclotronSectorFieldMap::read(file), 8, -0.05, 0.05, 1, {});
    sector->setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector(0.0), Quaternion(1, 0, 0, 0)));
    const GeneralElementMapping mapping(ElementList{sector}, {8, 8, 8}, reach);

    // Just past the sector's far edge at 60 degrees: outside its local x range, which is an
    // angular coordinate here, not a transverse one.
    const double theta = Physics::pi / 3.0;
    EXPECT_TRUE(mapping.coversTransversely(
            *sector, Vector({1.2 * std::cos(theta) - 1.2, 0.0, 1.2 * std::sin(theta)})));
    std::filesystem::remove(file);
}
