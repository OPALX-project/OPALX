//
// validate(): every check of a RUN that needs only the input values, in one place.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#ifndef OPALX_VALIDATE_RUN_H
#define OPALX_VALIDATE_RUN_H

#include <string>
#include <vector>

#include "Track/RunDescription.h"

namespace opalx::run {

    /**
     * @brief Check a run description before anything is allocated.
     *
     * Holds the checks that TrackRun::execute() and ParallelTracker::execute() used to make one
     * by one: missing objects, photon and decay beams, restart limits and beam stops. Each
     * message is the one the old check threw.
     *
     * @note Checks that need built elements stay where they are: ElementBase::checkRegions()
     *       needs placement and field-map extents, and ReferenceTracer::checkElementLengths()
     *       needs the traced reference momentum.
     * @return All errors found, empty if the run is valid.
     */
    std::vector<std::string> validate(const RunDescription& run);

    /// Join @p errors into one message, one error per line.
    std::string joinErrors(const std::vector<std::string>& errors);

}  // namespace opalx::run

#endif  // OPALX_VALIDATE_RUN_H
