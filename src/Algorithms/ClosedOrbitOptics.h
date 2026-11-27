#ifndef OPAL_CLOSED_ORBIT_OPTICS_H
#define OPAL_CLOSED_ORBIT_OPTICS_H

#include <cstdint>
#include "Algorithms/ClosedOrbitInitialState.h"
#include "Algorithms/ClosedOrbitSolver.h"

/** Successful static-magnetic RING solve, owned by value without Kokkos storage.
 * The matrix differentiates the fixed-energy section return in (x,px,y,py),
 * metres and mechanical p/(mc). normalMagneticField is B dot section-normal [T]
 * at the solved launch, sampled with the same evaluator as the ray tracker.
 * The generation identifies a particular execution of the named COF. A snapshot
 * describes that execution, never a live view of subsequently edited elements.
 * Failed reruns and unexecuted clones have no accessible optics snapshot.
 */
struct ClosedOrbitOptics {
    ClosedOrbitInitialState launch;
    OneTurnMap::Coordinates coordinates{}, residual{};
    OneTurnMap::Matrix matrix{};
    OneTurnMap::Coordinates scales{1, 1, 1, 1};
    /// Section Bx, By [T], retained for section/common-time sampling conversion.
    std::array<double, 2> transverseMagneticField{};
    double normalMagneticField = 0;
    double relativeEnergyDrift = 0;
    std::string command;
    std::uint64_t generation = 0;
};
#endif
