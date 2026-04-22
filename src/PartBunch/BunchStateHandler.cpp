#include "PartBunch/BunchStateHandler.h"

#include "Ippl.h"
#include "Utilities/Options.h"

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
}

std::shared_ptr<BunchStateHandler::ContainerState> BunchStateHandler::registerContainer() {
    auto state = std::make_shared<ContainerState>();
    registered_m.emplace_back(state);
    return state;
}

void BunchStateHandler::setUnitlessPositions(ContainerState& s, bool v) {
    s.unitlessPositions = Options::aggressiveStateSync ? syncOr(v) : v;
}

void BunchStateHandler::markMomentsDirty(ContainerState& s) {
    s.momentsDirty = Options::aggressiveStateSync ? syncOr(true) : true;
}

void BunchStateHandler::clearMomentsDirty(ContainerState& s) {
    // With AGGRESSIVE_STATE_SYNC, "clear" only takes effect if every rank
    // agrees the cache is clean; any rank still dirty keeps the whole bunch
    // dirty.
    s.momentsDirty = Options::aggressiveStateSync ? syncOr(false) : false;
}

void BunchStateHandler::setFirstRepartition(bool v) {
    firstRepartition_m = Options::aggressiveStateSync ? syncOr(v) : v;
}
