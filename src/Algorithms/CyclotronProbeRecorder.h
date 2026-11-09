// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPALX_CYCLOTRON_PROBE_RECORDER_H
#define OPALX_CYCLOTRON_PROBE_RECORDER_H
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "Algorithms/CyclotronProbe.h"
#include "PartBunch/PartBunchFwd.h"

class Probe;
class LossDataSink;

namespace cyclotron_probe::detail {
    /** Internal compact device-to-host payload, without host-only ownership.
     * Namespace scope keeps captured View value types accessible to NVCC extended
     * device lambdas. Recorder state and output ownership remain private.
     */
    struct Hit {
        Sample sample;
        std::int64_t id = -1;
        double charge = 0, mass = 0;
    };
}  // namespace cyclotron_probe::detail

/** @brief Shared passive recorder for global cyclotron and FFA PROBE segments.
 *
 * The historical class name does not restrict the crossing geometry: both
 * trackers supply accepted laboratory X-Z endpoints to the same device kernel.
 * A detector is an observer outside the magnetic lattice, so it cannot change
 * field selection, element placement, integration trials or the step schedule.
 *
 * Histories live in optional particle attributes, so migration and compaction
 * preserve ID/history association. Only compact hit records leave device memory.
 * No particle coordinates, fields, reference state or timestep are writable here.
 * save() is collective and explicit, never called from a destructor. No-hit runs
 * write zero histogram and header-only peaks when radial summaries are enabled.
 * FFA uses phase-space-only output: lab-origin radius is not the machine radius
 * of a translated ring, and no ring centre can be inferred from a segment.
 * The standard loss writer produces
 * no phase-space file when there are no hits. Restart/follow-up is rejected by
 * the tracker until diagnostic histories and output rewind are checkpointed.
 */
class CyclotronProbeRecorder {
public:
    struct Config {
        cyclotron_probe::Segment segment;
        double step = 0;              ///< Radial histogram bin width [m].
        std::string output;           ///< Suffix-free output basename.
        bool radialSummaries = true;  ///< Cyclotron lab-origin .hist/.peaks; false for FFA.
    };
    /// Validate geometry, positive STEP, and histogram size when summaries are enabled.
    static void validate(const Config& config);
    static Config configuration(const Probe& probe);
    explicit CyclotronProbeRecorder(std::vector<Config> configs);
    ~CyclotronProbeRecorder();
    /// Initialize/reset diagnostic histories for this tracking run (collective callers).
    void initialize(const std::vector<std::shared_ptr<ParticleContainer_t>>& containers);
    /// Observe synchronized accepted endpoints in lab coordinates, time [s].
    void observe(
            const std::shared_ptr<ParticleContainer_t>& pc, double time, short bunch, int turn);
    /// Collective flush: append pending HDF5/ASCII hits, replace cumulative
    /// .hist/.peaks summaries when enabled. Safe to repeat; histories persist.
    void save();

private:
    Kokkos::View<cyclotron_probe::detail::Hit*> hits_m;
    Kokkos::View<std::size_t> count_m{"PROBE hit count"};
    struct Diagnostic {
        Config config;
        double minimum = 0, maximum = 0;
        std::vector<double> histogram;
        // Reference turn -> (sum of radii [mm], particle count).
        std::map<int, std::pair<double, double>> peaks;
        std::unique_ptr<LossDataSink> sink;
    };
    std::vector<Diagnostic> diagnostics_m;
};
#endif
