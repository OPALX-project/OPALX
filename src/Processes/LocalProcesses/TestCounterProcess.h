//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
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
#ifndef OPAL_TESTCOUNTERPROCESS_H
#define OPAL_TESTCOUNTERPROCESS_H

#include "Processes/LocalProcesses/LocalProcess.h"

/**
 * @brief Physics-free local process that counts particles inside the element body.
 *
 * Used to validate the local-process plumbing (attachment, tracker invocation, frame
 * guarantees) without touching particle state. Counts particles whose element-local
 * longitudinal position lies in [0, element length).
 *
 * @note Counters are per-rank. apply() never modifies particle state and returns 0.
 */
class TestCounterProcess : public LocalProcess {
public:
    explicit TestCounterProcess(const std::string& name);

    std::size_t apply(
            Component& element, ParticleContainer<double, 3>& pc,
            const LocalProcessContext& context) override;

    /// Number of times apply() has been called.
    long long getInvocationCount() const;

    /// Particles inside the element body during the last apply() call.
    std::size_t getLastInsideCount() const;

    /// Sum of inside counts over all apply() calls.
    std::size_t getTotalInsideCount() const;

    /// Global track step passed in the context of the last apply() call (-1 if never called).
    long long getLastGlobalTrackStep() const;

private:
    /// Number of apply() invocations.
    long long invocationCount_m = 0;

    /// Inside count of the most recent invocation.
    std::size_t lastInsideCount_m = 0;

    /// Accumulated inside count over all invocations.
    std::size_t totalInsideCount_m = 0;

    /// Global track step of the most recent invocation.
    long long lastGlobalTrackStep_m = -1;
};

#endif
