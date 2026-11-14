//
// buildTracking(): make the particle bunch, output, samplers, space-charge solver and step-size
// schedule of one RUN.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#ifndef OPALX_BUILD_TRACKING_H
#define OPALX_BUILD_TRACKING_H

#include <memory>
#include <vector>

#include "Algorithms/TrackingSetup.h"
#include "PartBunch/PartBunch.h"
#include "SpaceCharge/SpaceChargeSolver.h"
#include "Structure/H5PartWrapper.h"
#include "Track/RunDescription.h"

class ElementBase;

namespace opalx::run {

    /**
     * @brief The objects buildTracking() makes for one RUN; the tracker borrows them.
     *
     * Members are destroyed in reverse order: the space-charge solver before the bunch it
     * borrows particles from, and both before the phase-space files.
     */
    struct RunObjects {
        /// Phase-space files, one per beam; the DataSink borrows them.
        std::vector<std::unique_ptr<H5PartWrapper>> phaseSpaceSinks;
        std::unique_ptr<PartBunch_t> bunch;
        std::unique_ptr<spacecharge::SpaceChargeSolver> spaceChargeSolver;
    };

    /**
     * @brief Build everything the tracker needs for this RUN: beams, particle containers,
     *        output, samplers or the checkpoint, the space-charge solver, the step-size
     *        schedule and the stop rules.
     * @param description What freeze() read for this RUN, checked by validate().
     * @param elements The elements buildElements() made for the tracked line.
     * @param objects Receives the bunch, the solver and the phase-space files; it must outlive
     *        the tracker.
     */
    TrackingSetup buildTracking(
            const RunDescription& description, std::vector<std::shared_ptr<ElementBase>> elements,
            RunObjects& objects);

}  // namespace opalx::run

#endif  // OPALX_BUILD_TRACKING_H
