// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPALX_PASSIVE_RING_PROBE_H
#define OPALX_PASSIVE_RING_PROBE_H

#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include "Algorithms/PassiveProbe.h"
#include "PartBunch/PartBunchFwd.h"

/** @brief Internal opt-in, passive selected-particle diagnostics for bare rings.
 *
 * This recorder reads accepted, synchronized endpoints. It never changes live
 * particle attributes, field evaluation, reference advancement or timestep
 * decisions. Selection uses stable nonnegative particle IDs, independent of
 * local array order and MPI ownership. Only selected endpoints leave the device.
 * MPI collectives assemble compact endpoints; rank zero interpolates and writes.
 *
 * OPALX_TEST_PASSIVE_RING_PROBES names a strict configuration file, with paths
 * relative to the run working directory. This is an internal experiment, not a
 * RUN input attribute. The caller must check eligibility (fresh aligned analytic
 * bare ring); cyclotron tracking must never instantiate this recorder.
 *
 * Configuration grammar (blank lines and whole-line # comments are allowed):
 * @code
 * OUTPUT passive-probes.csv
 * IDS 0 1 2
 * PLANE 0 0 0 0 0 0 1
 * @endcode
 * PLANE is followed by its integer ID, origin [m], normal, and an optional
 * positive arming tolerance [m] (default 1 nm). Normal length is arbitrary.
 * OUTPUT accepts one unquoted or double-quoted pathname. Existing output files
 * are rejected. Exactly one OUTPUT and IDS and at least one PLANE are required.
 *
 * Missing IDs are written explicitly with status=missing and NaN phase space;
 * their history resets, so reappearance cannot bridge an unobserved interval.
 * Reappearance preserves the observed crossing count, not an inferred turn count.
 * Unresolved intra-step crossings remain the caller's timestep-convergence duty.
 */
class PassiveRingProbe {
public:
    struct NamedPlane {
        std::int64_t id = 0;
        passive_probe::Plane plane;
    };
    struct Config {
        std::vector<std::int64_t> ids;
        std::vector<NamedPlane> planes;
        std::string output;
    };

    /// Strict, dependency-free parser. Throws std::invalid_argument on errors.
    static Config parseConfig(const std::string& text);
    /// Collective environment/configuration load. Unset environment returns null.
    static std::unique_ptr<PassiveRingProbe> fromEnvironment(bool eligible);
    /// Collective construction with identical validated config on every rank.
    explicit PassiveRingProbe(Config config);
    /// Collective observation; errors are propagated to every rank before throwing.
    void observe(const std::shared_ptr<ParticleContainer_t>& particles, double acceptedTime);

private:
    Config config_m;
    Kokkos::View<std::int64_t*> ids_m;
    Kokkos::View<double*> coordinates_m;
    Kokkos::View<int*> ownership_m;
    std::vector<passive_probe::State> histories_m;
    std::vector<bool> missing_m;
    std::ofstream output_m;
};

#endif
