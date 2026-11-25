#include <Kokkos_Core.hpp>
#include <type_traits>
#ifdef TPS_REVERSE_INCLUDE_ORDER
#include <tps/tps.hpp>
#include "Utilities/DivideError.h"
#else
#include "Utilities/DivideError.h"
#include <tps/tps.hpp>
#endif
static_assert(std::is_base_of_v<::OpalException, ::DivideError>);
int main() {
    Kokkos::ScopeGuard kokkos;

    bool opalCaught = false, tpsCaught = false;
    try {
        throw ::DivideError("coexistence");
    } catch (const ::ArithmeticError&) {
        opalCaught = true;
    }
    try {
        (void)(tps::Tps<double>(1.0) / 0.0);
    } catch (const ::OpalException&) {
        tpsCaught = true;
    }
    if (!opalCaught || !tpsCaught)
        return 2;
    auto x = tps::Tps<double>::makeVariable(2, 0);
    const auto f = (x + 2.0) * (x + 2.0);
    return f.getCoefficient(1) == 4.0 ? 0 : 1;
}
