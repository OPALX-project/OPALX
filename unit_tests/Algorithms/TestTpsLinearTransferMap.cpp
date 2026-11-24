#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/TpsLinearTransferMap.h"
#include "Ippl.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

namespace {
    using Matrix = matrix6x6_t;
    using Point  = std::array<double, 6>;
    using Series = tps::Tps<double>;

    Matrix identity() {
        Matrix result(0.0);
        for (unsigned i = 0; i < 6; ++i)
            result(i, i) = 1;
        return result;
    }
    // Independent dense arithmetic oracle, not OPALX prod() or TPS composition.
    Matrix multiply(const Matrix& a, const Matrix& b) {
        Matrix result(0.0);
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                for (unsigned k = 0; k < 6; ++k)
                    result(i, j) += a(i, k) * b(k, j);
        return result;
    }
    LinearTransferMapReference boundary(double z, double angle = 0) {
        const CoordinateSystemTrafo pose(
                Vector_t<double, 3>(3, -2, z),
                Quaternion(std::cos(angle / 2), 0, std::sin(angle / 2), 0));
        LinearTransferMapReference result;
        result.position   = pose.getOrigin();
        result.xAxis      = pose.rotateFrom(Vector_t<double, 3>(1, 0, 0));
        result.yAxis      = pose.rotateFrom(Vector_t<double, 3>(0, 1, 0));
        result.sAxis      = pose.rotateFrom(Vector_t<double, 3>(0, 0, 1));
        result.momentum   = 2.3 * result.sAxis;
        result.time       = 0.1 + z * 1e-9;
        result.pathLength = z;
        return result;
    }
    LinearTransferMapDraft segment(
            const Matrix& matrix, const LinearTransferMapReference& in,
            const LinearTransferMapReference& out) {
        LinearTransferMapDraft result;
        result.matrix   = matrix;
        result.entrance = in;
        result.exit     = out;
        return result;
    }
    void expectMatrix(const tps::Map<double>& map, const Matrix& expected) {
        const auto actual = map.jacobian(Point{});
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                EXPECT_NEAR(actual(i, j), expected(i, j), 2e-13) << i << ',' << j;
    }
    class TpsLinearTransferMapTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }
        static void TearDownTestSuite() { ippl::finalize(); }
        void SetUp() override { previousOrder = Series::getGlobalTruncOrder(); }
        void TearDown() override { Series::setGlobalTruncOrder(previousOrder); }
        int previousOrder;
    };
}  // namespace

TEST_F(TpsLinearTransferMapTest, CoefficientsChartsEvaluationAndSourcePreservation) {
    auto matrix = identity();
    for (unsigned i = 0; i < 6; ++i)
        for (unsigned j = 0; j < 6; ++j)
            matrix(i, j) += 0.01 * (i + 1) * (j + 2);
    auto source                          = segment(matrix, boundary(0, 0.4), boundary(1, -0.2));
    source.entrance.positionCorrection   = Vector_t<double, 3>(1e-6, -2e-6, 3e-6);
    source.entrance.timeCorrection       = 1e-12;
    source.entrance.pathLengthCorrection = 2e-13;
    source.finiteDifferenceSteps.fill(1e-4);
    source.richardsonLevels = 2;
    const auto saved        = source;
    const auto map          = opalx::maps::toTpsLinearTransferMap(source, "pass0:segment0");
    EXPECT_EQ(map.context().order(), 1);
    EXPECT_EQ(map.context().segments(), (std::vector<std::string>{"pass0:segment0"}));
    const auto& input = map.context().input();
    EXPECT_EQ(input.convention, tps::CoordinateConvention::mechanical);
    EXPECT_TRUE(input.canonicalPairs.empty());
    EXPECT_EQ(input.frame.name, map.context().output().frame.name);
    EXPECT_EQ(input.referenceTime, source.entrance.time - source.entrance.timeCorrection);
    ASSERT_TRUE(input.referenceMomentum);
    EXPECT_NEAR(*input.referenceMomentum, 2.3, 2e-13);
    EXPECT_EQ(input.momentumUnit, "p/(mc)");
    const std::array<std::string, 6> names{"x", "xprime", "y", "yprime", "zeta", "delta"};
    for (unsigned i = 0; i < 6; ++i) {
        EXPECT_EQ(input.coordinates[i].name, names[i]);
        EXPECT_EQ(input.coordinates[i].unit, i % 2 ? "1" : "m");
        EXPECT_EQ(input.coordinates[i].scale, 1);
        EXPECT_EQ(map.polynomial()[i][0], 0);
    }
    for (unsigned i = 0; i < 3; ++i) {
        EXPECT_EQ(
                input.frame.origin[i],
                source.entrance.position(i) - source.entrance.positionCorrection(i));
        EXPECT_EQ(input.frame.axes[i], source.entrance.xAxis(i));
        EXPECT_EQ(input.frame.axes[3 + i], source.entrance.yAxis(i));
        EXPECT_EQ(input.frame.axes[6 + i], source.entrance.sAxis(i));
        EXPECT_EQ(source.entrance.position(i), saved.entrance.position(i));
        EXPECT_EQ(source.entrance.positionCorrection(i), saved.entrance.positionCorrection(i));
    }
    EXPECT_EQ(source.entrance.timeCorrection, saved.entrance.timeCorrection);
    EXPECT_EQ(source.entrance.pathLengthCorrection, saved.entrance.pathLengthCorrection);
    EXPECT_EQ(source.finiteDifferenceSteps, saved.finiteDifferenceSteps);
    EXPECT_EQ(source.richardsonLevels, saved.richardsonLevels);
    expectMatrix(map, matrix);
    const Point point{0.01, -0.02, 0.03, -0.04, 0.05, -0.06};
    const auto actual = map.evaluate(point);
    for (unsigned i = 0; i < 6; ++i) {
        double expected = 0;
        for (unsigned j = 0; j < 6; ++j) {
            expected += matrix(i, j) * point[j];
            EXPECT_EQ(source.matrix(i, j), saved.matrix(i, j));
        }
        EXPECT_NEAR(actual[i], expected, 2e-13);
    }
}

TEST_F(TpsLinearTransferMapTest, NoncommutingCoupledCompositionAndAmbientOrder) {
    auto drift = identity(), lens = identity(), coupled = identity();
    drift(0, 1) = 0.3;
    drift(2, 3) = 0.3;
    lens(1, 0)  = -0.7;
    lens(3, 2)  = 0.7;
    for (unsigned i = 0; i < 6; ++i)
        for (unsigned j = 0; j < 6; ++j)
            coupled(i, j) += 0.005 * (i + 1) * (6 - j);
    const auto a        = segment(drift, boundary(0), boundary(1, 0.1));
    const auto b        = segment(lens, a.exit, boundary(2, 0.3));
    const auto c        = segment(coupled, b.exit, boundary(3, -0.2));
    const auto expected = multiply(coupled, multiply(lens, drift));
    EXPECT_GT(std::abs(multiply(lens, drift)(0, 0) - multiply(drift, lens)(0, 0)), 0.1);
    for (int order : {0, 1, 4, Series::EXACT}) {
        Series::setGlobalTruncOrder(order);
        const std::vector<tps::Map<double>> maps{
                opalx::maps::toTpsLinearTransferMap(a, "a"),
                opalx::maps::toTpsLinearTransferMap(b, "b"),
                opalx::maps::toTpsLinearTransferMap(c, "c")};
        const auto combined = tps::composeSequence<double>(maps, 1);
        expectMatrix(combined, expected);
        EXPECT_EQ(combined.context().input(), maps.front().context().input());
        EXPECT_EQ(combined.context().output(), maps.back().context().output());
        EXPECT_EQ(combined.context().segments(), (std::vector<std::string>{"a", "b", "c"}));
        EXPECT_EQ(Series::getGlobalTruncOrder(), order);
    }
}

TEST_F(TpsLinearTransferMapTest, RejectsMismatchedBoundaryMetadata) {
    const auto first = segment(identity(), boundary(0), boundary(1));
    const auto inner = opalx::maps::toTpsLinearTransferMap(first, "inner");
    auto second      = segment(identity(), first.exit, boundary(2));
    EXPECT_NO_THROW(
            (void)tps::compose(opalx::maps::toTpsLinearTransferMap(second, "outer"), inner));
    for (unsigned change = 0; change < 4; ++change) {
        auto mismatch = second;
        if (change == 0) mismatch.entrance.position(0) += 0.01;
        if (change == 1) mismatch.entrance.time += 1e-9;
        if (change == 2) mismatch.entrance.momentum *= 1.1;
        if (change == 3) {
            const auto rotated      = boundary(1, 0.3);
            mismatch.entrance.xAxis = rotated.xAxis;
            mismatch.entrance.yAxis = rotated.yAxis;
            mismatch.entrance.sAxis = rotated.sAxis;
        }
        const auto outer = opalx::maps::toTpsLinearTransferMap(mismatch, "outer");
        EXPECT_THROW((void)tps::compose(outer, inner), std::invalid_argument);
        EXPECT_EQ(Series::getGlobalTruncOrder(), previousOrder);
    }
    for (bool changeScale : {false, true}) {
        auto output = inner.context().output();
        if (changeScale)
            output.coordinates[0].scale = 1000;
        else
            output.coordinates[0].unit = "mm";
        const tps::Map<double> mismatched(
                tps::MapContext(inner.context().input(), output, 1, {"inner"}), inner.polynomial());
        EXPECT_THROW(
                (void)tps::compose(
                        opalx::maps::toTpsLinearTransferMap(second, "outer"), mismatched),
                std::invalid_argument);
    }
}

TEST_F(TpsLinearTransferMapTest, KeepsReferenceSlopesAtFixedPlaneReturn) {
    auto source          = segment(identity(), boundary(0), boundary(1));
    source.exit.momentum = Vector_t<double, 3>(0.2, -0.1, 1.8);
    const auto map       = opalx::maps::toTpsLinearTransferMap(source, "return");
    EXPECT_NEAR(map.context().output().reference[1], 0.2 / 1.8, 2e-13);
    EXPECT_NEAR(map.context().output().reference[3], -0.1 / 1.8, 2e-13);
    EXPECT_EQ(map.context().output().frame.axes[8], 1);
    for (double value : map.evaluate(Point{}))
        EXPECT_EQ(value, 0);
}

TEST_F(TpsLinearTransferMapTest, TraversalIdsAreIndependentOfOwnerNames) {
    auto a           = segment(identity(), boundary(0), boundary(1));
    auto b           = segment(identity(), a.exit, boundary(2));
    a.activeElements = b.activeElements = {"Q"};
    const auto first                    = opalx::maps::toTpsLinearTransferMap(a, "pass0:0");
    EXPECT_THROW(
            (void)tps::compose(opalx::maps::toTpsLinearTransferMap(b, "pass0:0"), first),
            std::invalid_argument);
    const auto joined = tps::compose(opalx::maps::toTpsLinearTransferMap(b, "pass0:1"), first);
    EXPECT_EQ(joined.context().segments(), (std::vector<std::string>{"pass0:0", "pass0:1"}));
}

TEST_F(TpsLinearTransferMapTest, RejectsInvalidMatrixAndReferences) {
    const auto valid = segment(identity(), boundary(0), boundary(1));
    EXPECT_THROW(opalx::maps::toTpsLinearTransferMap(valid, ""), OpalException);
    for (unsigned change = 0; change < 8; ++change) {
        auto invalid = valid;
        if (change == 0) invalid.matrix(0, 0) = std::numeric_limits<double>::quiet_NaN();
        if (change == 1) invalid.matrix(5, 5) = std::numeric_limits<double>::infinity();
        if (change == 2) invalid.entrance.momentum = Vector_t<double, 3>(0.0);
        if (change == 3) invalid.exit.momentum = Vector_t<double, 3>(1, 0, 0);
        if (change == 4) invalid.entrance.xAxis *= 2;
        if (change == 5) invalid.exit.timeCorrection = std::numeric_limits<double>::infinity();
        if (change == 6)
            invalid.exit.positionCorrection(0) = std::numeric_limits<double>::quiet_NaN();
        if (change == 7)
            invalid.exit.pathLengthCorrection = std::numeric_limits<double>::infinity();
        EXPECT_THROW(opalx::maps::toTpsLinearTransferMap(invalid, "invalid"), OpalException)
                << change;
    }
}

TEST_F(TpsLinearTransferMapTest, PersistentRecordOwnsPolynomialAndCompensatedMetadata) {
    auto matrix = identity();
    matrix(0, 3) = .27;
    matrix(5, 1) = -.19;
    auto draft = segment(matrix, boundary(0, .2), boundary(1, -.3));
    draft.entrance.positionCorrection = {1e-17, 2e-17, -3e-17};
    draft.exit.timeCorrection = 1e-20;
    draft.exit.pathLengthCorrection = -2e-17;
    draft.segment = 7;
    draft.activeElements = {"FRINGE"};
    const auto stored = std::make_shared<const LinearTransferMap>(draft, "record:7");
    draft.matrix(0, 3) = 99;
    draft.activeElements.clear();
    expectMatrix(stored->polynomial(), matrix);
    EXPECT_EQ(stored->id(), "record:7");
    EXPECT_EQ(stored->metadata().segment, 7u);
    EXPECT_EQ(stored->metadata().activeElements, (std::vector<std::string>{"FRINGE"}));
    EXPECT_EQ(stored->metadata().entrance.positionCorrection(2), -3e-17);
    EXPECT_EQ(stored->metadata().exit.timeCorrection, 1e-20);
    EXPECT_EQ(stored->metadata().exit.pathLengthCorrection, -2e-17);
    const auto jacobian = stored->matrix();
    for (unsigned i = 0; i < 6; ++i)
        for (unsigned j = 0; j < 6; ++j) EXPECT_EQ(jacobian(i, j), matrix(i, j));
}
