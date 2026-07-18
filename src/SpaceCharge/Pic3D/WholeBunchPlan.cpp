/**
 * @file WholeBunchPlan.cpp
 * @brief Implements direct whole-bunch PIC iteration.
 */

#include "SpaceCharge/Pic3D/WholeBunchPlan.h"

#include "Utilities/OpalException.h"

#include <array>
#include <limits>
#include <string>

namespace opalx::spacecharge {

    template <typename T, unsigned Dim>
    PreparedIteration WholeBunchPlan<T, Dim>::prepare(
            generation_type particleGeneration, BinConfigurationObserver*) {
        advanceEpoch();
        unitPending_m = true;
        return PreparedIteration{
                IterationKind::WholeBunch, particleGeneration, 0, preparationEpoch_m};
    }

    template <typename T, unsigned Dim>
    bool WholeBunchPlan<T, Dim>::hasNext(
            const PreparedIteration& prepared, generation_type currentParticleGeneration) {
        validatePrepared(prepared, currentParticleGeneration);
        return unitPending_m;
    }

    template <typename T, unsigned Dim>
    std::optional<typename WholeBunchPlan<T, Dim>::Unit> WholeBunchPlan<T, Dim>::next(
            const PreparedIteration& prepared, generation_type currentParticleGeneration) {
        validatePrepared(prepared, currentParticleGeneration);
        if (!unitPending_m) {
            return std::nullopt;
        }
        unitPending_m = false;

        const size_type localCount  = particles_m.getLocalNum();
        const size_type globalCount = particles_m.getTotalNum();
        return Unit(
                IterationKind::WholeBunch, SolveUnitFieldMode::Direct, 0, localCount, globalCount,
                Selection::direct(0, localCount), true, std::array<double, Dim>{}, 1.0);
    }

    template <typename T, unsigned Dim>
    void WholeBunchPlan<T, Dim>::validatePrepared(
            const PreparedIteration& prepared, generation_type currentParticleGeneration) const {
        if (prepared.kind != IterationKind::WholeBunch
            || prepared.preparationEpoch != preparationEpoch_m) {
            throw OpalException(
                    "WholeBunchPlan::next",
                    "The prepared iteration does not belong to the current whole-bunch epoch.");
        }
        if (prepared.particleGeneration != currentParticleGeneration) {
            throw OpalException(
                    "WholeBunchPlan::next",
                    "Particle storage changed after the whole-bunch iteration was prepared.");
        }
        if (prepared.mergedBinCount != 0) {
            throw OpalException(
                    "WholeBunchPlan::next", "The whole-bunch iteration token is malformed.");
        }
    }

    template <typename T, unsigned Dim>
    void WholeBunchPlan<T, Dim>::advanceEpoch() {
        if (preparationEpoch_m == std::numeric_limits<std::uint64_t>::max()) {
            throw OpalException(
                    "WholeBunchPlan::prepare", "The iteration preparation epoch overflowed.");
        }
        ++preparationEpoch_m;
    }

    template class WholeBunchPlan<double, 3>;

}  // namespace opalx::spacecharge
