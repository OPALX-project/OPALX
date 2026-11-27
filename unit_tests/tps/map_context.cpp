#include <Kokkos_Core.hpp>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <tps/tps.hpp>
#include "Utilities/CLRangeError.h"
using S = tps::Tps<double>;
void require(bool condition, const char* msg) {
    if (!condition) throw std::runtime_error(msg);
}
void near(double a, double b) {
    require(std::isfinite(a) && std::abs(a - b) < 2e-13, "map context numerical mismatch");
}
template <class F>
void rejects(F f) {
    bool caught = false;
    try {
        f();
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    require(caught, "expected context rejection");
}
tps::PhaseSpaceChart chart(std::string name) {
    tps::PhaseSpaceChart c;
    c.name              = std::move(name);
    c.coordinates       = {{"q", "m", 1}, {"p", "normalized", 1}};
    c.reference         = {0, 0};
    c.frame.name        = "lab";
    c.referenceTime     = 0;
    c.referenceMomentum = 1;
    c.momentumUnit      = "mc";
    c.convention        = tps::CoordinateConvention::canonical;
    c.canonicalPairs    = {{0, 1}};
    return c;
}
int main() {
    Kokkos::ScopeGuard kokkos;
    try {
        auto a = chart("A"), b = chart("B"), c = chart("C");
        // Fresh executable: construction/composition works without global-order setup.
        require(S::getGlobalTruncOrder() == S::EXACT, "default scalar order");
        //! [context map composition]
        auto q = S::makeVariable(2, 0), p = S::makeVariable(2, 1);
        std::array<S, 2> drift{q + 2.0 * p, p}, kick{q, p - 3.0 * q * q};
        double& retained = drift[0][0];
        tps::Map<double> d(tps::MapContext(a, b, 4, {"drift"}), drift),
                k(tps::MapContext(b, c, 4, {"kick"}), kick);
        const std::array<double, 2> x{0.1, 0.2};
        auto combined = tps::compose(k, d, 4);
        auto y        = combined.evaluate(x);
        near(y[0], 0.5);
        near(y[1], -0.55);
        auto sequential = k.evaluate(d.evaluate(x));
        near(y[0], sequential[0]);
        near(y[1], sequential[1]);
        require(combined.context().input() == a && combined.context().output() == c,
                "composed charts");
        require(combined.context().segments() == std::vector<std::string>({"drift", "kick"}),
                "segment tracking order");
        auto jac = combined.jacobian(x);
        near(jac(0, 0), 1);
        near(jac(0, 1), 2);
        near(jac(1, 0), -3);
        near(jac(1, 1), -5);
        //! [context map composition]
        auto id   = tps::identityMap<double>(a, 4);
        auto same = tps::compose(d, id, 4).evaluate(x);
        near(same[0], 0.5);
        near(same[1], 0.2);
        auto rightId = tps::identityMap<double>(b, 4);
        same         = tps::compose(rightId, d, 4).evaluate(x);
        near(same[0], 0.5);
        std::array<tps::Map<double>, 2> sequence{d, k};
        auto total = tps::composeSequence<double>(sequence, 4);
        near(total.evaluate(x)[1], -0.55);
        rejects([&] {
            (void)tps::composeSequence<double>({}, 4);
        });
        rejects([&] {
            (void)tps::compose(d, k, 4);
        });
        rejects([&] {
            (void)tps::compose(k, d, 5);
        });
        rejects([&] {
            (void)tps::compose(k, d, -1);
        });
        rejects([&] {
            (void)tps::identityMap<double>(a, 0);
        });
        rejects([&] {
            (void)d.evaluate(std::array<double, 1>{0});
        });
        rejects([&] {
            (void)d.jacobian(std::array<double, 1>{0});
        });
        // Precision is explicit and operations restore the caller's scalar setting.
        S::setGlobalTruncOrder(0);
        auto projected = tps::compose(k, d, 1);
        require(S::getGlobalTruncOrder() == 0, "global precision restoration");
        near(projected.evaluate(x)[1], 0.2);
        near(projected.jacobian(x)(1, 1), 1);
        require(S::getGlobalTruncOrder() == 0, "Jacobian restored precision");
        auto minimum = tps::compose(k, d.truncated(2));
        require(minimum.context().order() == 2, "minimum context order");
        auto constantProjection = d.truncated(0);
        near(constantProjection.evaluate(x)[0], 0);
        S::setGlobalTruncOrder(4);
        // An affine inner map requires explicit stored-polynomial composition.
        auto shifted = drift;
        shifted[0] += 1.0;
        tps::Map<double> affine(tps::MapContext(a, b, 4), shifted);
        rejects([&] {
            (void)tps::compose(k, affine, 4);
        });
        auto affineResult = tps::compose(k, affine, 4, tps::CompositionPolicy::storedPolynomial);
        near(affineResult.evaluate(x)[1], 0.2 - 3 * 1.5 * 1.5);
        auto translatedConstant =
                tps::compose(k, affine, 0, tps::CompositionPolicy::storedPolynomial);
        near(translatedConstant.evaluate(x)[0], 1);
        near(translatedConstant.evaluate(x)[1], -3);
        std::array<tps::Map<double>, 2> affineSequence{affine, k};
        auto constantSequence = tps::composeSequence<double>(
                affineSequence, 0, tps::CompositionPolicy::storedPolynomial);
        near(constantSequence.evaluate(x)[1], -3);
        // Constant outputs retain input dimension through context and have zero derivatives.
        std::array<S, 2> constants{S(2.0), S(4.0)};
        tps::Map<double> constant(tps::MapContext(b, c, 4), constants);
        auto values = tps::compose(constant, d, 4).evaluate(x);
        near(values[0], 2);
        near(values[1], 4);
        near(constant.truncated(0).evaluate(x)[0], 2);
        auto zeros = constant.jacobian(x);
        near(zeros(0, 0), 0);
        near(zeros(1, 1), 0);
        // Each boundary field contributes to exact join compatibility.
        for (int field = 0; field < 9; ++field) {
            auto bad = b;
            switch (field) {
                case 0:
                    bad.coordinates[0].unit = "cm";
                    break;
                case 1:
                    bad.coordinates[0].scale = 2;
                    break;
                case 2:
                    bad.reference[0] = 1;
                    break;
                case 3:
                    bad.frame.origin[0] = 1;
                    break;
                case 4:
                    bad.referenceTime = 1;
                    break;
                case 5:
                    bad.referenceMomentum = 2;
                    break;
                case 6:
                    bad.coordinates[0].name = "other";
                    break;
                case 7:
                    bad.canonicalPairs = {{1, 0}};
                    break;
                case 8:
                    bad.parameters = {{"gauge", "1/m", 2}};
                    break;
            }
            tps::Map<double> wrong(tps::MapContext(bad, c, 4), kick);
            rejects([&] {
                (void)tps::compose(wrong, d, 4);
            });
        }
        // Named chart parameters are values, units and names, not display-only labels.
        auto parameterized       = b;
        parameterized.parameters = {{"gauge", "1/m", 2}};
        auto parameterInner      = tps::Map<double>(tps::MapContext(a, parameterized, 4), drift);
        auto parameterOuter      = tps::Map<double>(tps::MapContext(parameterized, c, 4), kick);
        near(tps::compose(parameterOuter, parameterInner, 4).evaluate(x)[1], -0.55);
        for (int field = 0; field < 3; ++field) {
            auto changed = parameterized;
            if (field == 0) changed.parameters[0].value = 3;
            if (field == 1) changed.parameters[0].name = "different";
            if (field == 2) changed.parameters[0].unit = "1/cm";
            auto wrong = tps::Map<double>(tps::MapContext(changed, c, 4), kick);
            rejects([&] {
                (void)tps::compose(wrong, parameterInner, 4);
            });
        }
        for (int field = 0; field < 4; ++field) {
            auto invalid = parameterized;
            if (field == 0) invalid.parameters[0].value = std::numeric_limits<double>::infinity();
            if (field == 1) invalid.parameters[0].name.clear();
            if (field == 2) invalid.parameters[0].unit.clear();
            if (field == 3) invalid.parameters.push_back(invalid.parameters.front());
            rejects([&] {
                (void)tps::MapContext(invalid, c, 4);
            });
        }
        tps::Map<double> duplicate(tps::MapContext(b, c, 4, {"drift"}), kick);
        rejects([&] {
            (void)tps::compose(duplicate, d, 4);
        });
        require(S::getGlobalTruncOrder() == 4, "rejection changed precision");
        rejects([&] {
            (void)tps::MapContext(a, b, 4, {"x", "x"});
        });
        rejects([&] {
            (void)tps::MapContext(a, b, -1);
        });
        for (int field = 0; field < 10; ++field) {
            auto invalid = a;
            switch (field) {
                case 0:
                    invalid.coordinates.clear();
                    break;
                case 1:
                    invalid.reference.clear();
                    break;
                case 2:
                    invalid.coordinates[0].scale = 0;
                    break;
                case 3:
                    invalid.reference[0] = std::numeric_limits<double>::quiet_NaN();
                    break;
                case 4:
                    invalid.frame.axes[0] = 2;
                    break;
                case 5:
                    invalid.frame.axes[0] = -1;
                    break;
                case 6:
                    invalid.coordinates[1].name = "q";
                    break;
                case 7:
                    invalid.canonicalPairs = {{0, 0}};
                    break;
                case 8:
                    invalid.referenceMomentum = -1;
                    break;
                case 9:
                    invalid.frame.name.clear();
                    break;
            }
            rejects([&] {
                (void)tps::MapContext(invalid, b, 4);
            });
        }
        rejects([&] {
            (void)tps::Map<double>(tps::MapContext(a, b, 1), kick);
        });
        rejects([&] {
            std::array<S, 1> missing{q};
            (void)tps::Map<double>(tps::MapContext(a, b, 4), missing);
        });
        auto insufficient = drift;
        S::setGlobalTruncOrder(1);
        insufficient[0] = q * q;
        rejects([&] {
            (void)tps::Map<double>(tps::MapContext(a, b, 4), insufficient);
        });
        require(S::getGlobalTruncOrder() == 1, "constructor restored precision");
        // A rectangular context map joins a two-variable and a one-variable chart.
        auto one = a;
        one.name = "one";
        one.coordinates.resize(1);
        one.reference.resize(1);
        one.convention = tps::CoordinateConvention::unspecified;
        one.canonicalPairs.clear();
        S::setGlobalTruncOrder(4);
        std::array<S, 1> sum{q + p};
        tps::Map<double> reduction(tps::MapContext(a, one, 4), sum);
        auto oneVar = S::makeVariable(1, 0);
        std::array<S, 2> expand{oneVar, 2.0 * oneVar};
        tps::Map<double> expansion(tps::MapContext(one, c, 4), expand);
        auto rectangle = tps::compose(expansion, reduction, 4).evaluate(x);
        near(rectangle[0], 0.3);
        near(rectangle[1], 0.6);
        // Composition must be safe when nonlinear products grow a fresh cache.
        auto three = a;
        three.name = "three";
        three.coordinates.push_back({"z", "m", 1});
        three.reference.push_back(0);
        three.convention = tps::CoordinateConvention::unspecified;
        three.canonicalPairs.clear();
        auto u = S::makeVariable(3, 0), v = S::makeVariable(3, 1), w = S::makeVariable(3, 2);
        std::array<S, 3> squares{u * u, v, w};
        tps::Map<double> nonlinear(tps::MapContext(three, three, 4), squares);
        auto fourth = tps::compose(nonlinear, nonlinear, 4);
        near(fourth.evaluate(std::array<double, 3>{0.5, 0.2, 0.3})[0], 0.0625);
        // Copy-on-write map ownership isolates later edits to source coefficients.
        retained = 99;
        near(d.evaluate(x)[0], 0.5);
        drift[0][0] = 88;
        near(d.evaluate(x)[0], 0.5);
        // An exception inside algebra/cache setup also restores global order.
        auto largeIdentity = tps::identityMap<double>(a, 65536);
        S::setGlobalTruncOrder(1);
        bool overflowRejected = false;
        try {
            (void)tps::compose(largeIdentity, largeIdentity, 65536);
        } catch (const ::CLRangeError&) {
            overflowRejected = true;
        }
        require(overflowRejected && S::getGlobalTruncOrder() == 1,
                "failed composition restored precision");
        S::setGlobalTruncOrder(S::EXACT);
        std::cout << "Map context and composition checks passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
