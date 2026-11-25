#include "Utilities/LogicalError.h"
#include "Utilities/OpalException.h"
// Intentionally do not include the umbrella: TpsMath must supply its own dependencies.
#include <tps/Algebra/TpsMath.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>

static_assert(std::is_base_of_v<::OpalException, ::LogicalError>);
static_assert(!std::is_base_of_v<std::exception, ::OpalException>);
int main() {
    Kokkos::ScopeGuard runtime;
    try {
        // Force elementary-function templates to instantiate outside tps.hpp.
        const auto x = tps::Tps<double>::makeVariable(1, 0);
        const auto f = tps::sqrt(x + 4.0);
        if (std::abs(f[0] - 2.0) > 2e-13 || std::abs(f[1] - 0.25) > 2e-13)
            throw std::runtime_error("standalone elementary-function header");
        try {
            throw ::LogicalError("test operation", "dimension mismatch");
        } catch (const ::OpalException& error) {
            if (error.where() != "test operation" ||
                std::string(error.what()) != "dimension mismatch")
                throw std::runtime_error("algebra exception metadata");
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
