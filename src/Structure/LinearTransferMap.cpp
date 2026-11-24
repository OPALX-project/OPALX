// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Structure/LinearTransferMap.h"
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <utility>
#include "Algorithms/TpsLinearTransferMap.h"
#include "Utilities/OpalException.h"

struct StoredTransferMap::State {
    std::optional<tps::Map<double>> map;
    explicit State(tps::Map<double> value) : map(std::move(value)) {}
};

std::vector<std::weak_ptr<StoredTransferMap::State>>& StoredTransferMap::registry() {
    static std::vector<std::weak_ptr<State>> states;
    return states;
}

StoredTransferMap::StoredTransferMap(tps::Map<double> map) {
    if (!Kokkos::is_initialized() || Kokkos::is_finalized())
        throw OpalException("StoredTransferMap", "An active Kokkos runtime is required.");
    if (map.context().input().coordinates.size() != 6
        || map.context().output().coordinates.size() != 6)
        throw OpalException("StoredTransferMap", "Expected a six-dimensional transfer map.");
    static bool hookRegistered = false;
    if (!hookRegistered) {
        Kokkos::push_finalize_hook([] {
            releaseTransferMapStorage();
            hookRegistered = false;
        });
        hookRegistered = true;
    }
    auto& states = registry();
    // Do not retain old control blocks through repeated map recalculation.
    std::erase_if(states, [](const auto& state) {
        return state.expired();
    });
    state_m = std::make_shared<State>(std::move(map));
    states.push_back(state_m);
}

void releaseTransferMapStorage() {
    for (const auto& weak : StoredTransferMap::registry())
        if (auto state = weak.lock()) state->map.reset();
    StoredTransferMap::registry().clear();
}

bool StoredTransferMap::available() const {
    return state_m && state_m->map.has_value() && Kokkos::is_initialized()
           && !Kokkos::is_finalized();
}

const tps::Map<double>& StoredTransferMap::polynomial() const {
    if (!available())
        throw OpalException(
                "StoredTransferMap::polynomial", "Map payload was released or runtime ended.");
    return *state_m->map;
}

matrix6x6_t StoredTransferMap::matrix() const {
    const auto derivative = polynomial().jacobian(std::array<double, 6>{});
    matrix6x6_t result(0.0);
    for (unsigned i = 0; i < 6; ++i)
        for (unsigned j = 0; j < 6; ++j)
            result(i, j) = derivative(i, j);
    return result;
}

LinearTransferMap::LinearTransferMap(
        const LinearTransferMapDraft& draft, const std::string& occurrenceId)
    : metadata_m(draft),
      id_m(occurrenceId),
      payload_m(opalx::maps::toTpsLinearTransferMap(draft, occurrenceId)) {}

TransferMapPass::TransferMapPass(
        std::string id, std::vector<LinearTransferMapHandle> segments, tps::Map<double> combined,
        double determinant, double symplectic)
    : id_m(std::move(id)),
      segments_m(std::move(segments)),
      combined_m(std::move(combined)),
      determinantResidual_m(determinant),
      symplecticResidual_m(symplectic) {}
