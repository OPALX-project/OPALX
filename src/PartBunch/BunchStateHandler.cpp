#include "PartBunch/BunchStateHandler.h"

#include "Ippl.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <cmath>
#include <functional>
#include <utility>

namespace {
    // Converge a local bool across all MPI ranks using logical-OR: any rank reporting the
    // unitless state forces the whole bunch to adopt it.
    inline bool syncOr(bool v) {
        bool out = v;
        ippl::Comm->allreduce(v, out, 1, std::logical_or<bool>());
        return out;
    }
}  // namespace

// -- BunchStateHandler::ContainerState -------------------------------------

void BunchStateHandler::ContainerState::setUnitlessPositions(bool v) {
    unitlessPositions = Options::aggressiveStateSync ? syncOr(v) : v;
}

// -- BunchStateHandler -----------------------------------------------------

std::shared_ptr<BunchStateHandler::ContainerState> BunchStateHandler::registerContainer() {
    auto state = std::make_shared<ContainerState>();
    registered_m.emplace_back(state);
    return state;
}

void BunchStateHandler::setFixedCartesianDomain(
        std::array<double, 3> lower, std::array<double, 3> upper) {
    for (std::size_t dimension = 0; dimension < lower.size(); ++dimension) {
        if (!std::isfinite(lower[dimension]) || !std::isfinite(upper[dimension])
            || lower[dimension] >= upper[dimension]) {
            throw OpalException(
                    "BunchStateHandler::setFixedCartesianDomain",
                    "Fixed Cartesian domain bounds must be finite and strictly increasing.");
        }
    }

    FixedCartesianDomainState requested{std::move(lower), std::move(upper)};
    if (fixedCartesianDomain_m.has_value()) {
        if (*fixedCartesianDomain_m == requested) {
            return;
        }
        throw OpalException(
                "BunchStateHandler::setFixedCartesianDomain",
                "Clear the active fixed Cartesian domain before replacing its bounds.");
    }
    fixedCartesianDomain_m = std::move(requested);
}
