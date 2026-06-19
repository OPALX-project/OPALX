#include <gtest/gtest.h>

#include <array>
#include <memory>
#include <string>

#include "Ippl.h"
#include "PartBunch/FieldAccumulator.hpp"
#include "PartBunch/FieldContainer.hpp"
#include "PartBunch/FieldOps.hpp"

namespace {

    class FieldOpsAndAccumulatorTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }

        std::unique_ptr<FieldContainer<double, 3>> makeFieldContainer() {
            ippl::Vector<int, 3> nr        = 6;
            ippl::Vector<double, 3> rmin   = -3.0;
            ippl::Vector<double, 3> rmax   = 3.0;
            ippl::Vector<double, 3> origin = rmin;
            ippl::Vector<double, 3> hr     = (rmax - rmin) / ippl::Vector<double, 3>(nr);
            std::array<bool, 3> decomp     = {true, true, true};

            ippl::NDIndex<3> domain;
            for (unsigned i = 0; i < 3; ++i) {
                domain[i] = ippl::Index(nr[i]);
            }

            auto fields = std::make_unique<FieldContainer<double, 3>>(
                    hr, rmin, rmax, decomp, domain, origin, true);
            fields->initializeFields("OPEN");
            return fields;
        }

        std::shared_ptr<VField_t<double, 3>> makeVectorField(
                FieldContainer<double, 3>& fields) {
            auto field = std::make_shared<VField_t<double, 3>>();
            field->initialize(fields.getMesh(), fields.getFL());
            return field;
        }

        double maxScalarError(Field_t<3>& field, double expected) {
            auto view     = field.getView();
            double maxErr = 0.0;
            ippl::parallel_reduce(
                    "TestFieldOpsAndAccumulator::maxScalarError", field.getFieldRangePolicy(),
                    KOKKOS_LAMBDA(
                            const ippl::RangePolicy<3>::index_array_type& idx, double& local) {
                        const double diff = apply(view, idx) - expected;
                        const double err  = diff < 0.0 ? -diff : diff;
                        if (err > local) {
                            local = err;
                        }
                    },
                    Kokkos::Max<double>(maxErr));
            return maxErr;
        }

        double maxVectorError(
                VField_t<double, 3>& field, const Vector_t<double, 3>& expected) {
            auto view     = field.getView();
            double maxErr = 0.0;
            ippl::parallel_reduce(
                    "TestFieldOpsAndAccumulator::maxVectorError", field.getFieldRangePolicy(),
                    KOKKOS_LAMBDA(
                            const ippl::RangePolicy<3>::index_array_type& idx, double& local) {
                        const Vector_t<double, 3> value = apply(view, idx);
                        for (unsigned d = 0; d < 3; ++d) {
                            const double diff = value[d] - expected[d];
                            const double err  = diff < 0.0 ? -diff : diff;
                            if (err > local) {
                                local = err;
                            }
                        }
                    },
                    Kokkos::Max<double>(maxErr));
            return maxErr;
        }
    };

    TEST_F(FieldOpsAndAccumulatorTest, ScalarOpsSetScaleAndShiftField) {
        auto fields = makeFieldContainer();

        opalx::fieldops::setScalarField(fields->getRho(), 2.0);
        opalx::fieldops::scaleAndShiftScalarField(fields->getRho(), 3.0, -1.0);

        EXPECT_NEAR(maxScalarError(fields->getRho(), 5.0), 0.0, 1e-12);
    }

    TEST_F(FieldOpsAndAccumulatorTest, AccumulatorClearsAndAccumulatesUnflippedField) {
        auto fields = makeFieldContainer();
        auto eTmp   = makeVectorField(*fields);
        auto bTmp   = makeVectorField(*fields);

        const Vector_t<double, 3> ePrime{1.0, 2.0, 3.0};
        opalx::fieldops::setVectorField(fields->getE(), ePrime);
        opalx::fieldops::setVectorField(*eTmp, Vector_t<double, 3>{9.0, 9.0, 9.0});
        opalx::fieldops::setVectorField(*bTmp, Vector_t<double, 3>{9.0, 9.0, 9.0});

        FieldAccumulator<double, 3> accumulator;
        accumulator.bind(eTmp, bTmp);
        accumulator.clear();
        accumulator.accumulate(
                fields->getE(), 1.0, Vector_t<double, 3>{0.0, 0.0, 0.0});

        EXPECT_NEAR(maxVectorError(*eTmp, ePrime), 0.0, 1e-12);
        EXPECT_NEAR(maxVectorError(*bTmp, Vector_t<double, 3>{0.0, 0.0, 0.0}), 0.0, 1e-12);
    }

    TEST_F(FieldOpsAndAccumulatorTest, AccumulatorAppliesFlippedImageSignRule) {
        auto fields = makeFieldContainer();
        auto eTmp   = makeVectorField(*fields);
        auto bTmp   = makeVectorField(*fields);

        opalx::fieldops::setVectorField(fields->getE(), Vector_t<double, 3>{1.0, 2.0, 3.0});

        FieldAccumulator<double, 3> accumulator;
        accumulator.bind(eTmp, bTmp);
        accumulator.clear();
        accumulator.accumulate(
                fields->getE(), 1.0, Vector_t<double, 3>{0.0, 0.0, 0.0},
                -1.0, 2);

        EXPECT_NEAR(
                maxVectorError(*eTmp, Vector_t<double, 3>{-1.0, -2.0, 3.0}), 0.0, 1e-12);
        EXPECT_NEAR(maxVectorError(*bTmp, Vector_t<double, 3>{0.0, 0.0, 0.0}), 0.0, 1e-12);
    }

}  // namespace
