#include <Kokkos_Core.hpp>
#include <iostream>
#include <stdexcept>
#include <string>
#include "tps/Algebra/Tps.hpp"

using namespace tps;
using TPS = Tps<double>;
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
int main(int argc, char** argv) {
    Kokkos::ScopeGuard kokkos;

    try {
        const int variables    = argc > 1 ? std::stoi(argv[1]) : 2;
        const bool seriesFirst = argc > 2;
        // Each CTest invocation is a fresh process; no prior positive-order cache.
        TpsMonomial zero(variables);
        if (seriesFirst) {
            TPS constant(0, variables);
            constant[0] = 7;
            require(constant[zero] == 7, "fresh zero-order TPS indexing");
        } else
            require(zero.getIndex() == 0, "fresh zero-order monomial index");
        TpsData* data = TpsData::getTpsData(0, variables);
        require(data->getSize(-1) == 0, "negative degree has no coefficients");
        require(data->getSize(0) == 1, "degree zero has one coefficient");
        require(data->getExponents(0).getOrder() == 0, "constant exponent table");
        require(data->getProductArray(0)[0] == 0, "constant product table");
        require(data->getSubTable().size() == 1, "constant substitution table");
        TPS a(0, variables);
        a[zero]       = 3;
        const TPS& ca = a;
        require(ca[zero] == 3, "const indexing");
        require((a * a)[zero] == 9, "zero-order multiplication");
        Vector<double> point(variables, 0.2);
        require(a.evaluate(point) == 3, "zero-order evaluation");
        Matrix<double> identity(variables, variables, 0.0);
        for (int j = 0; j < variables; ++j)
            identity(j, j) = 1;
        require(a.substitute(identity)[0] == 3, "zero-order substitution");
        // Grow the same cache and then reuse an existing zero-order object.
        auto x    = TPS::makeVariable(variables, 0);
        auto cube = x * x * x;
        TpsMonomial third(variables);
        third[0] = 3;
        require(cube[third] == 1, "cache growth");
        require(zero.getIndex() == 0 && a[zero] == 3, "constant index after growth");
        TpsMonomial empty;
        require(empty.getIndex() == 0, "zero-variable constant index");
        TPS scalar(4.0);
        require(scalar[empty] == 4, "scalar monomial indexing");
        std::cout << "Zero-order cache checks passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
