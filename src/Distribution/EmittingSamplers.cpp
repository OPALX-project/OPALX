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
#include "Distribution/EmittingSamplers.h"

#include "Distribution/SamplingBase.hpp"
#include "Elements/OpalBeamline.h"
#include "PartBunch/PartBunch.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <functional>
#include <limits>
#include <utility>

EmittingSamplers::EmittingSamplers(PerContainer samplers) : samplers_m(std::move(samplers)) {}

EmittingSamplers::InitialReference EmittingSamplers::initialReference(std::size_t ci) const {
    InitialReference offsets;
    if (ci < samplers_m.size()) {
        for (const auto& sampler : samplers_m[ci]) {
            if (!sampler) {
                continue;
            }
            if (!offsets.momentum && sampler->hasInitialReferenceMomentum()) {
                offsets.momentum = sampler->getInitialReferenceMomentum();
            }
            if (!offsets.position && sampler->hasInitialReferencePosition()) {
                offsets.position = sampler->getInitialReferencePosition();
            }
        }
    }
    return offsets;
}

void EmittingSamplers::emit(PartBunch_t& bunch, const OpalBeamline& beamline, double t, double dt) {
    const auto& containers = bunch.getParticleContainers();
    for (size_t ci = 0; ci < containers.size(); ++ci) {
        const auto& pc = containers[ci];
        if (!pc) {
            continue;
        }
        if (bunch.pcAtSStop(ci)) {
            continue;
        }

        // Record the extent of the position array before emission. If an internal resize
        // (Kokkos::realloc) happens during emission, the extent of R will change and we can
        // flag this as an error.
        const size_t extentBeforeEmission = pc->R.size();

        // Samplers append new particles in the beamline frame and do not touch the others, so
        // only the new particles change frame.
        const size_t nBefore = pc->getLocalNum();
        if (ci < samplers_m.size()) {
            for (const auto& sampler : samplers_m[ci]) {
                if (sampler) {
                    sampler->emitParticles(t, dt);
                }
            }
        }
        pc->adoptParticles(beamline.beamlineFrame(*pc), nBefore);

        // Sanity guard: the total number of macroparticles in the bunch must
        // never exceed the globally configured BEAM::NALLOC value. Overshooting
        // this limit would trigger internal reallocations in the particle
        // container and silently drop already-tracked particles/delete their data in the particle
        // attributes. This is only a check for the number of local particles.
        const size_t extentAfterEmission = pc->R.size();
        if (extentAfterEmission != extentBeforeEmission) {
            throw OpalException(
                    "EmittingSamplers::emit",
                    "Local particle storage was resized during emission (likely due to "
                    "over-emission causing a Kokkos::realloc). This leads to loss of "
                    "previously tracked particles. Please increase the total number of "
                    "macroparticles or adjust NPARTDIST / the emission profile. If you are "
                    "using emission sources, please check the emission profile and adjust "
                    "the number of particles emitted.");
        }
    }
    bunch.refreshPcActiveAfterEmit();
}

void EmittingSamplers::activateEmittingContainers(PartBunch_t& bunch, double t) const {
    for (size_t i = 0; i < samplers_m.size(); ++i) {
        for (const auto& sampler : samplers_m[i]) {
            if (sampler && !sampler->isEmissionDone(t)) {
                bunch.setPcActive(i);
                break;
            }
        }
    }
}

std::optional<double> EmittingSamplers::emissionTimeStep(double t) const {
    double emissionDt  = std::numeric_limits<double>::max();
    bool hasEmissionDt = false;

    for (const auto& samplers : samplers_m) {
        for (const auto& sampler : samplers) {
            if (!sampler || sampler->isEmissionDone(t)) {
                continue;
            }
            const double samplerDt = sampler->getEmissionTimeStep();
            if (samplerDt > 0.0) {
                emissionDt    = std::min(emissionDt, samplerDt);
                hasEmissionDt = true;
            }
        }
    }

    if (hasEmissionDt) {
        return emissionDt;
    }
    return std::nullopt;
}

opalx::spacecharge::EmissionProgress EmittingSamplers::spaceChargeEmissionProgress(
        double t, bool aggressiveStateSync) const {
    opalx::spacecharge::EmissionProgress progress;

    // While emission is active, CartesianPIC3D stretches its beam-frame mesh over the full source
    // pulse length, matching old OPAL. For example, 5 percent emission produces approximately a
    // factor-20 longitudinal stretch. Without it, early charge is compressed onto an artificially
    // short mesh and can receive excessive self-field kicks, including kicks back into the source.
    // The least-complete active source determines the required stretch.
    for (const auto& samplers : samplers_m) {
        for (const auto& sampler : samplers) {
            if (!sampler || sampler->isEmissionDone(t)) {
                continue;
            }
            const double samplerFraction = sampler->getEmittedFraction();
            if (samplerFraction < 1.0) {
                progress.active   = true;
                progress.fraction = std::min(progress.fraction, samplerFraction);
            }
        }
    }

    if (aggressiveStateSync) {
        // Preserve legacy state convergence across ranks: any active source keeps stretching
        // enabled, and the minimum emitted fraction supplies the most conservative domain length.
        bool globallyActive = progress.active;
        ippl::Comm->allreduce(progress.active, globallyActive, 1, std::logical_or<bool>());
        const double localFraction =
                progress.active ? std::clamp(progress.fraction, 0.0, 1.0) : 1.0;
        double globalFraction = localFraction;
        ippl::Comm->allreduce(localFraction, globalFraction, 1, std::less<double>());
        progress.active   = globallyActive;
        progress.fraction = globallyActive ? globalFraction : 1.0;
    } else {
        progress.fraction = progress.active ? std::clamp(progress.fraction, 0.0, 1.0) : 1.0;
    }
    return progress;
}
