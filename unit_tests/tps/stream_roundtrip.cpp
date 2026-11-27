#include <Kokkos_Core.hpp>
#include <complex>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include "tps/Algebra/Tps.hpp"

namespace {
    void require(bool valid, const char* message) {
        if (!valid) throw std::runtime_error(message);
    }
    template <class T>
    void same(const tps::Tps<T>& actual, const tps::Tps<T>& expected) {
        require(actual.getVariables() == expected.getVariables(), "variable count");
        require(actual.getMaxOrder() == expected.getMaxOrder(), "stored degree");
        require(actual.getTruncOrder() == expected.getTruncOrder(), "precision");
        require(actual.getSize() == expected.getSize(), "coefficient count");
        for (int i = 0; i < expected.getSize(); ++i)
            require(actual[i] == expected[i], "round-trip coefficient changed");
    }
    template <class T>
    void check(int variables, int order) {
        using S = tps::Tps<T>;
        S::setGlobalTruncOrder(S::EXACT);
        S polynomial(order, variables);
        tps::TpsMonomial power(variables);
        power[0]          = order;
        polynomial[power] = T(1);
        polynomial[0]     = T(-2.5);
        if (variables > 1) {
            power[0]             = 0;
            power[variables - 1] = order;
            polynomial[power]    = T(.125);
        }
        if constexpr (std::is_same_v<T, std::complex<double>>) polynomial[1] = T(1.25, -.5);
        // Exactly representable coefficients isolate token boundaries from the
        // writer's existing 14-decimal-digit coefficient formatting policy.
        const S scalar(T(1.25));
        const auto linear = S::makeVariable(variables, 0);
        std::stringstream stream;
        stream << polynomial << scalar << linear;
        S restored, restoredScalar, restoredLinear;
        stream >> restored >> restoredScalar >> restoredLinear;
        require(!stream.fail(), "stream extraction failed");
        same(restored, polynomial);
        same(restoredScalar, scalar);
        same(restoredLinear, linear);
        stream >> std::ws;
        require(stream.eof(), "unexpected trailing records");
    }
}  // namespace
int main() {
    Kokkos::ScopeGuard runtime;
    try {
        for (int order : {2, 99, 100, 101}) {
            check<double>(1, order);
            check<std::complex<double>>(1, order);
        }
        check<double>(2, 100);
        check<std::complex<double>>(2, 100);
        check<double>(6, 2);
        // Previously emitted whitespace-delimited input remains readable.
        std::istringstream old("Tps 2 2147483647 1\n1.0  2\n0.0 -1\n");
        tps::Tps<double> restored;
        old >> restored;
        require(!old.fail() && restored[2] == 1 && restored[0] == 0, "old format input");
        std::cout << "11 coefficient-exact round trips of three consecutive objects passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
