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
#include "Processes/LocalProcesses/TestCounterProcess.h"

#include "AbsBeamline/Component.h"
#include "PartBunch/ParticleContainer.hpp"

namespace {
    using pc_size_type = ippl::detail::size_type;
}

TestCounterProcess::TestCounterProcess(const std::string& name) : LocalProcess(name) {}

std::size_t TestCounterProcess::apply(
        Component& element, ParticleContainer<double, 3>& pc,
        const LocalProcessContext& context) {
    // Host-side parameter extraction; only PODs are captured by the kernel.
    const double length        = element.getElementLength();
    const pc_size_type nLocal  = pc.getLocalNum();
    auto Rview                 = pc.R.getView();

    pc_size_type insideCount = 0;
    Kokkos::parallel_reduce(
            "TestCounterProcess::countInside", nLocal,
            KOKKOS_LAMBDA(const pc_size_type i, pc_size_type& count) {
                const double z = Rview(i)[2];
                count += (z >= 0.0 && z < length) ? 1 : 0;
            },
            insideCount);
    Kokkos::fence();

    ++invocationCount_m;
    lastInsideCount_m = static_cast<std::size_t>(insideCount);
    totalInsideCount_m += static_cast<std::size_t>(insideCount);
    lastGlobalTrackStep_m = context.globalTrackStep;

    return 0;
}

long long TestCounterProcess::getInvocationCount() const {
    return invocationCount_m;
}

std::size_t TestCounterProcess::getLastInsideCount() const {
    return lastInsideCount_m;
}

std::size_t TestCounterProcess::getTotalInsideCount() const {
    return totalInsideCount_m;
}

long long TestCounterProcess::getLastGlobalTrackStep() const {
    return lastGlobalTrackStep_m;
}
