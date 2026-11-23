//
// Class EmittingSamplers
//   The samplers that emit particles into the containers of a bunch during tracking.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPALX_EMITTING_SAMPLERS_H
#define OPALX_EMITTING_SAMPLERS_H

#include "OPALTypes.h"
#include "SpaceCharge/SpaceChargeSolveContext.h"

#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

class OpalBeamline;
class SamplingBase;

/**
 * @brief The samplers that emit particles into the containers of a bunch during tracking.
 *
 * Holds, per container, the samplers of the emission sources that emit during tracking
 * (buildTracking() makes them). A sampler appends the particles born in a step in the beamline
 * frame; emit() brings them into the frame of their container. The tracker owns this object and
 * passes it the bunch and the beamline it borrows.
 */
class EmittingSamplers {
public:
    /// Per container, its emitting samplers.
    using PerContainer = std::vector<std::vector<std::shared_ptr<SamplingBase>>>;

    EmittingSamplers() = default;
    explicit EmittingSamplers(PerContainer samplers);

    /// Emission offsets R0 [m] and P0 [beta*gamma] in the beamline frame, from the first sampler
    /// of a container that reports each.
    struct InitialReference {
        std::optional<Vector_t<double, 3>> position;
        std::optional<Vector_t<double, 3>> momentum;
    };

    /// The emission offsets the samplers of container @p ci report; empty if none does.
    [[nodiscard]] InitialReference initialReference(std::size_t ci) const;

    /**
     * @brief Let the samplers of each container not stopped at its s-stop append the particles
     *        born in the step from @p t [s] over @p dt [s], bring the new particles from the
     *        beamline frame into the container's frame and switch on the containers that now
     *        hold particles.
     * @throws OpalException if the local particle storage of a container was resized during
     *         emission, which loses the particles tracked so far.
     */
    void emit(PartBunch_t& bunch, const OpalBeamline& beamline, double t, double dt);

    /// Switch on each container with a sampler that still emits at time @p t [s].
    void activateEmittingContainers(PartBunch_t& bunch, double t) const;

    /// The smallest positive emission time step [s] of the samplers that still emit at time
    /// @p t [s]; empty if none asks for one.
    [[nodiscard]] std::optional<double> emissionTimeStep(double t) const;

    /**
     * @brief Whether a sampler still emits at time @p t [s], and the least emitted fraction of
     *        those, for the space-charge mesh stretch during emission.
     * @param aggressiveStateSync OPTION AGGRESSIVE_STATE_SYNC: reduce both over the ranks
     *        (collective).
     */
    [[nodiscard]] opalx::spacecharge::EmissionProgress spaceChargeEmissionProgress(
            double t, bool aggressiveStateSync) const;

private:
    PerContainer samplers_m;
};

#endif  // OPALX_EMITTING_SAMPLERS_H
