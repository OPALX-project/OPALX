#include <iostream>
#include <stdexcept>
#include <tps/tps.hpp>
int main() {
    try {
        bool rejected = false;
        try {
            tps::Array1D<int> invalid(1);
        } catch (const std::logic_error&) {
            rejected = true;
        }
        if (!rejected) return 1;
        {
            Kokkos::ScopeGuard runtime;
            auto x = tps::Tps<double>::makeVariable(6, 0);
            auto f = x * x;
            if (f.getMaxOrder() != 2) return 2;
            tps::Array1D<tps::Array1D<int>> nested(2);
            nested[0] = tps::Array1D<int>(3, 7);
            if (nested[0][2] != 7) return 3;
            // Static monomial cache owns views here; its finalize hook must release them.
        }
        rejected = false;
        try {
            tps::Array1D<int> invalid(1);
        } catch (const std::logic_error&) {
            rejected = true;
        }
        if (!rejected) return 4;
        std::cout << "Runtime and nested-view lifetime checks passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 5;
    }
}
