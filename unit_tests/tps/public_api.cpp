// Intentionally do not include the umbrella: TpsMath must supply its own dependencies.
#include <tps/Algebra/TpsMath.h>
#include <tps/Error.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>

static_assert(std::is_base_of_v<tps::Error, tps::SizeError>);
static_assert(std::is_base_of_v<std::runtime_error, tps::Error>);
int main() {
    Kokkos::ScopeGuard runtime;
    try {
        // Force elementary-function templates to instantiate outside tps.hpp.
        const auto x = tps::Tps<double>::makeVariable(1, 0);
        const auto f = tps::sqrt(x + 4.0);
        if (std::abs(f[0] - 2.0) > 2e-13 || std::abs(f[1] - 0.25) > 2e-13)
            throw std::runtime_error("standalone elementary-function header");
        try {
            throw tps::SizeError("test operation", "dimension mismatch");
        } catch (const tps::Error& error) {
            if (error.where() != "test operation" ||
                std::string(error.what()) != "dimension mismatch")
                throw std::runtime_error("algebra exception metadata");
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
