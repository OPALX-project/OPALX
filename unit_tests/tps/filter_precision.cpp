#include <Kokkos_Core.hpp>
#include <array>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "tps/Algebra/Tps.hpp"

namespace {
    using S = tps::Tps<double>;
    void require(bool valid, const char* message) {
        if (!valid) throw std::runtime_error(message);
    }
    S polynomial(int variables, int precision) {
        // Explicit whitespace makes this fixture independent of the stream writer.
        std::ostringstream text;
        text << "Tps 2 " << precision << ' ' << variables << "\n2";
        for (int i = 0; i < variables; ++i)
            text << " 0";
        text << "\n1 2";
        for (int i = 1; i < variables; ++i)
            text << " 0";
        text << "\n0";
        for (int i = 0; i < variables; ++i)
            text << " -1";
        std::istringstream input(text.str());
        S result;
        input >> result;
        require(!input.fail(), "fixture parse");
        return result;
    }
}  // namespace

int main() {
    Kokkos::ScopeGuard runtime;
    try {
        for (int variables : {1, 2, 6})
            for (int precision : {4, S::EXACT})
                for (int cap : {0, 1, 2, 5, S::EXACT})
                    for (const auto range :
                         {std::array{0, 2}, std::array{2, 2}, std::array{3, 5}, std::array{0, 0}}) {
                        S::setGlobalTruncOrder(S::EXACT);
                        const auto source = polynomial(variables, precision);
                        tps::TpsMonomial quadratic(variables);
                        quadratic[0] = 2;
                        S::setGlobalTruncOrder(cap);
                        const auto result = source.filter(range[0], range[1]);
                        require(S::getGlobalTruncOrder() == cap, "filter changed global cap");
                        require(result.getVariables() == variables, "filter changed variables");
                        require(result.getMaxOrder() == std::min(2, range[1]), "stored degree");
                        S::setGlobalTruncOrder(S::EXACT);
                        require(result.getTruncOrder() == precision, "stored precision lost");
                        require(source.getTruncOrder() == precision, "source precision changed");
                        const double expected = range[0] <= 2 && range[1] >= 2 ? 1 : 0;
                        const auto sum        = result + S::makeVariable(variables, 0);
                        const double actual   = sum.getMaxOrder() >= 2 ? sum[quadratic] : 0;
                        require(actual == expected, "arithmetic lost retained quadratic");
                        require(result[0] == (range[0] == 0 ? 2 : 0), "constant selection");
                        require(source[quadratic] == 1 && source[0] == 2, "source changed");
                    }
        // Scalars and the truncate alias preserve precision as well.
        S::setGlobalTruncOrder(0);
        auto scalar    = S(2.).filter(0, 3);
        auto truncated = polynomial(1, 4).truncate(2);
        S::setGlobalTruncOrder(S::EXACT);
        require(scalar.getVariables() == 0 && scalar.getTruncOrder() == S::EXACT,
                "scalar stored precision");
        require(truncated.getTruncOrder() == 4 && truncated[2] == 1, "truncate precision");
        std::cout << "Stored filter precision: 120 range/cap/dimension cases passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
