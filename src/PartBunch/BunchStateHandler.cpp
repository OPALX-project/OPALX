#include "PartBunch/BunchStateHandler.h"

#include "Ippl.h"
#include "Utilities/Options.h"

#include <algorithm>
#include <functional>

namespace {
    // Converge a local bool across all MPI ranks using logical-OR. Intentionally
    // chosen over logical-AND so that any rank reporting the "set" state forces
    // the whole bunch to adopt it: dirty wins over clean, unitless wins over
    // physical, first-repartition wins over done.
    inline bool syncOr(bool v) {
        bool out = v;
        ippl::Comm->allreduce(v, out, 1, std::logical_or<bool>());
        return out;
    }

    inline double syncMin(double v) {
        double out = v;
        ippl::Comm->allreduce(v, out, 1, std::less<double>());
        return out;
    }
}  // namespace

// -- BunchStateHandler::ContainerState -------------------------------------

void BunchStateHandler::ContainerState::setUnitlessPositions(bool v) {
    unitlessPositions = Options::aggressiveStateSync ? syncOr(v) : v;
}

void BunchStateHandler::ContainerState::markMomentsDirty() {
    momentsDirty = Options::aggressiveStateSync ? syncOr(true) : true;
}

void BunchStateHandler::ContainerState::markMomentsClean() {
    // With AGGRESSIVE_STATE_SYNC, "clear" only takes effect if every rank
    // agrees the cache is clean; any rank still dirty keeps the whole bunch
    // dirty.
    momentsDirty = Options::aggressiveStateSync ? syncOr(false) : false;
}

// -- BunchStateHandler -----------------------------------------------------

std::shared_ptr<BunchStateHandler::ContainerState> BunchStateHandler::registerContainer() {
    auto state = std::make_shared<ContainerState>();
    registered_m.emplace_back(state);
    return state;
}

void BunchStateHandler::setFirstRepartition(bool v) {
    firstRepartition_m = Options::aggressiveStateSync ? syncOr(v) : v;
}

void BunchStateHandler::setEmissionMeshProgress(bool active, double emittedFraction) {
    const bool syncedActive = Options::aggressiveStateSync ? syncOr(active) : active;
    double fraction         = active ? std::clamp(emittedFraction, 0.0, 1.0) : 1.0;

    if (Options::aggressiveStateSync) {
        fraction = syncMin(fraction);
    }

    emissionMeshStretchEnabled_m   = syncedActive;
    emissionMeshProgressFraction_m = syncedActive ? fraction : 1.0;
}
