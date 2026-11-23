//
// Class BunchReference
//   The reference particles of the containers of a bunch: their start, the reference pass of the
//   primary beam and their advance in each step.
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
#ifndef OPALX_BUNCH_REFERENCE_H
#define OPALX_BUNCH_REFERENCE_H

#include "Algorithms/CavityAutophaser.h"
#include "OPALTypes.h"
#include "Utility/IpplTimings.h"

#include <functional>
#include <memory>
#include <optional>
#include <vector>

class Inform;
class OpalBeamline;
class SamplingBase;
class StepSizeConfig;

/**
 * @brief The reference particles of the containers of a bunch.
 *
 * Every container that is not a decay daughter has a reference particle: a lab-frame position
 * and momentum, a path length and the reference-to-lab transform of its particles, all stored in
 * the container. A decay daughter container copies those of its root. This class sets them for a
 * fresh run, moves them to the start position, traces the primary reference once through the
 * lattice and advances them in every step.
 *
 * Borrows the bunch and the beamline of the tracker; owns the cavity phases that autophase sets
 * and a checkpoint stores, and whether a reference particle hit an element.
 */
class BunchReference {
public:
    /**
     * @param beamReferenceMomenta Per container, the reference momentum [beta*gamma] set on BEAM;
     *        an empty entry uses the bunch mean.
     * @param autoPhase Value of OPTION AUTOPHASE.
     * @param cavityPhases Phases a checkpoint restored; autophase adds to them.
     */
    BunchReference(
            PartBunch_t& bunch, OpalBeamline& beamline,
            std::vector<std::optional<Vector_t<double, 3>>> beamReferenceMomenta, int autoPhase,
            CavityPhases cavityPhases);

    /**
     * @brief Set the lab-frame pose of each beam's reference particle and its reference-to-lab
     *        transform for a fresh run; decay daughter containers are skipped. Collective.
     *
     * Position and momentum come from the bunch mean, else from the emission offsets (R0, P0)
     * of the sampler, else from the design pose; PC, ENERGY or GAMMA on BEAM sets the momentum.
     * @param emittingSamplers Per container, the samplers that emit during tracking.
     */
    void setInitialPoses(
            const std::vector<std::vector<std::shared_ptr<SamplingBase>>>& emittingSamplers);

    /**
     * @brief Integrate the references in time from t = 0 until the primary reaches path length
     *        @p sStart, autophasing the cavities on the way; the last step is fractional.
     * @param stepSizes Tracking schedule; the search walks a copy of it from its first segment.
     * @param selectTimeStep Sets the bunch time step for the DT of a segment; the search calls it
     *        when it enters a segment.
     * @return False if the primary is already at or past @p sStart and nothing moved.
     */
    bool findStartPositions(
            double sStart, const StepSizeConfig& stepSizes,
            const std::function<void(double)>& selectTimeStep);

    /**
     * @brief Trace the primary beam's reference particle once through the lattice
     *        (ReferenceTracer: autophase, element-length and reach checks), starting a bunch
     *        length behind it. Collective.
     * @param minTimeStep Time step of the pass [s], the smallest DT of @p stepSizes.
     */
    void tracePrimary(const StepSizeConfig& stepSizes, double minTimeStep);

    /**
     * @brief Advance the reference particles of the tracked containers over the bunch time step
     *        and move their reference-to-lab transforms onto them; decay daughter containers
     *        take their root's frame change. The references move only once the bunch time is
     *        positive.
     */
    void advance();

    /**
     * @brief Copy each decay daughter container's root reference particle, frame and path
     *   length into it, without moving particles.
     * @note Used after the reference particles are set up; during tracking the daughter
     *   containers follow their root in advance().
     */
    void copyFollowers();

    /// Whether a reference particle hit an element in advance() on this rank; not reduced over
    /// the ranks.
    [[nodiscard]] bool hitElement() const { return hitElement_m; }

    /// Phases autophase set (and a restart restored), for the checkpoint.
    [[nodiscard]] const CavityPhases& cavityPhases() const { return cavityPhases_m; }

private:
    /// Advance reference positions and momenta through the beamline for one step.
    void updateReferenceParticles();

    /// Move each container's reference-to-lab transform onto its advanced reference particle.
    void updateRefToLabCSTrafo();

    /// Autophase the RFCAVITY and TRAVELINGWAVE elements at the primary reference.
    void autophaseCavities();

    /// Log the reference state of each container at track start. Collective.
    void printInitialContainerRefs(Inform& m) const;

    PartBunch_t& bunch_m;
    OpalBeamline& beamline_m;
    std::vector<std::optional<Vector_t<double, 3>>> beamReferenceMomenta_m;
    int autoPhase_m;
    CavityPhases cavityPhases_m;
    bool hitElement_m = false;
    IpplTimings::TimerRef referencePassTimer_m;
};

#endif  // OPALX_BUNCH_REFERENCE_H
