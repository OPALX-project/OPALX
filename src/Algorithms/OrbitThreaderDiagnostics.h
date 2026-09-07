// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_ORBIT_THREADER_DIAGNOSTICS_H
#define OPAL_ORBIT_THREADER_DIAGNOSTICS_H

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <ostream>
#include "Utility/IpplTimings.h"

/// Internal host-only instrumentation; does not change the ray-tracking public API.
namespace orbit_threader_diagnostics {
    enum Phase { reference, segmentation, maps, output, phaseCount };
    inline constexpr std::array<const char*, phaseCount> phaseNames{
        "OT reference", "OT segmentation", "OT maps", "OT output"};

    /** Work performed in a phase, including discarded localization trials.
     * nominalSteps counts outer tracking iterations; advances also includes event
     * localization calls. acceptedSteps counts recursive leaves (or direct steps),
     * including leaves of trajectories later discarded by an enclosing bisection.
     * fieldSamples counts summed-field evaluator calls; elementFields counts the
     * individual element applications. Lookup tests count all scanned elements.
     */
    struct Work {
        std::uint64_t nominalSteps{0}, advances{0}, trials{0}, acceptedSteps{0};
        std::uint64_t capSplits{0}, supportSplits{0}, fieldSamples{0}, elementFields{0};
        std::uint64_t supportLookups{0}, bodyLookups{0}, elementTests{0};
        std::uint64_t membershipReuses{0}, exitIterations{0};
        std::uint64_t samples{0}, segments{0}, rays{0}, logRows{0};
        unsigned maxDepth{0};
        double minAcceptedDt{std::numeric_limits<double>::infinity()};
        double stepCap{std::numeric_limits<double>::infinity()};
        double raySeconds{0.0}, logSeconds{0.0};
        double rayTransportSeconds{0.0}, rayExitSeconds{0.0};
    };
    struct Report { std::array<Work, phaseCount> work; };

    // A null context makes shared beamline/ray calls outside threading uninstrumented.
    // No atomics, MPI collectives, device memory, or OpenMP regions are introduced.
    // If ray tracking is parallelized later, worker-local reports must be merged explicitly.
    inline thread_local Report* activeReport = nullptr;
    inline thread_local Work* activeWork = nullptr;

    inline void count(std::uint64_t Work::*member, std::uint64_t n = 1) {
        if (activeWork) activeWork->*member += n;
    }
    inline void accepted(double dt) {
        if (activeWork) {
            ++activeWork->acceptedSteps;
            activeWork->minAcceptedDt = std::min(activeWork->minAcceptedDt, std::abs(dt));
        }
    }

    class Session {
        Report* previousReport = activeReport;
        Work* previousWork = activeWork;
    public:
        explicit Session(Report& report) { activeReport = &report; activeWork = nullptr; }
        ~Session() { activeReport = previousReport; activeWork = previousWork; }
        Session(const Session&) = delete;
        Session& operator=(const Session&) = delete;
    };

    /// Coarse IPPL timers accumulate across species/TRACK calls in timing.dat.
    /// They use the existing configured timer-fence policy; never used per substep.
    class TimedPhase {
        Work* previousWork = activeWork;
        IpplTimings::TimerRef timer{};
        bool running = false;
    public:
        explicit TimedPhase(Phase phase) {
            if (!activeReport) return;
            timer = IpplTimings::getTimer(phaseNames[phase]);
            activeWork = &activeReport->work[phase];
            IpplTimings::startTimer(timer);
            running = true;
        }
        void stop() {
            if (!running) return;
            IpplTimings::stopTimer(timer);
            activeWork = previousWork;
            running = false;
        }
        ~TimedPhase() { stop(); }
        TimedPhase(const TimedPhase&) = delete;
        TimedPhase& operator=(const TimedPhase&) = delete;
    };

    /// Local wall-clock subset of a phase, without IPPL fences or stored samples.
    class Duration {
        using Clock = std::chrono::steady_clock;
        double* total;
        Clock::time_point start;
    public:
        explicit Duration(double Work::*member)
            : total(activeWork ? &(activeWork->*member) : nullptr),
              start(total ? Clock::now() : Clock::time_point{}) {}
        void stop() {
            if (!total) return;
            *total += std::chrono::duration<double>(Clock::now() - start).count();
            total = nullptr;
        }
        ~Duration() { stop(); }
        Duration(const Duration&) = delete;
        Duration& operator=(const Duration&) = delete;
    };

    inline void print(std::ostream& out, const Report& report) {
        out << "* OrbitThreader work (this pass, rank 0 local; accepted includes localization trials):\n";
        for (unsigned phase = 0; phase < phaseCount; ++phase) {
            const auto& w = report.work[phase];
            out << "* " << phaseNames[phase]
                << ": nominal_steps=" << w.nominalSteps << " advances=" << w.advances
                << " trials=" << w.trials << " accepted_steps=" << w.acceptedSteps
                << " cap_splits=" << w.capSplits << " support_splits=" << w.supportSplits
                << " max_depth=" << w.maxDepth
                << " min_accepted_dt_s=" << (w.acceptedSteps ? w.minAcceptedDt : 0.0)
                << " step_cap_s=" << (w.advances ? w.stepCap : 0.0) << '\n'
                << "* " << phaseNames[phase]
                << ": support_lookups=" << w.supportLookups << " body_lookups=" << w.bodyLookups
                << " element_tests=" << w.elementTests << " field_samples=" << w.fieldSamples
                << " element_fields=" << w.elementFields
                << " membership_reuses=" << w.membershipReuses << '\n'
                << "* " << phaseNames[phase]
                << ": reference_samples=" << w.samples << " segments=" << w.segments
                << " rays=" << w.rays << " ray_tracking_s=" << w.raySeconds
                << " ray_transport_s=" << w.rayTransportSeconds
                << " ray_exit_s=" << w.rayExitSeconds << " exit_iterations=" << w.exitIterations
                << " log_rows=" << w.logRows << " log_write_s=" << w.logSeconds << '\n';
        }
    }
}
#endif
