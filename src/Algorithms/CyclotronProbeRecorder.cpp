// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/CyclotronProbeRecorder.h"
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <stdexcept>
#include "AbsBeamline/Probe.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"

namespace {
    std::pair<double, double> radialBounds(const cyclotron_probe::Segment& segment) {
        const cyclotron_probe::Vector d = segment.end - segment.start;
        const double u                  = std::clamp(-dot(segment.start, d) / dot(d, d), 0.0, 1.0);
        const cyclotron_probe::Vector closest = segment.start + u * d;
        return {std::hypot(closest(0), closest(2)),
                std::max(
                        std::hypot(segment.start(0), segment.start(2)),
                        std::hypot(segment.end(0), segment.end(2)))};
    }

}  // namespace

void CyclotronProbeRecorder::validate(const Config& config) {
    cyclotron_probe::State state;
    cyclotron_probe::Sample sample;
    if (config.output.empty() || !std::isfinite(config.step) || config.step <= 0
        || cyclotron_probe::update(config.segment, state, 0, {}, sample)
                   != cyclotron_probe::Status::Initialized)
        throw std::invalid_argument(
                "PROBE requires finite distinct endpoints, positive STEP and output name");
    const auto [minimum, maximum] = radialBounds(config.segment);
    // Bound allocation and MPI counts before conversion. 10 million bins already
    // require 80 MB per rank per probe; smaller STEP is not a useful default.
    const double bins = std::ceil((maximum - minimum) / config.step) + 1;
    if (!std::isfinite(bins) || bins > 10000000)
        throw std::invalid_argument("PROBE histogram exceeds 10000000 bins; increase STEP");
}

CyclotronProbeRecorder::Config CyclotronProbeRecorder::configuration(const Probe& probe) {
    Config config;
    config.segment.start = cyclotron_probe::Vector(probe.getXStart(), 0.0, probe.getYStart());
    config.segment.end   = cyclotron_probe::Vector(probe.getXEnd(), 0.0, probe.getYEnd());
    config.step          = probe.getStep();
    config.output        = probe.getOutputFN();
    validate(config);
    return config;
}

CyclotronProbeRecorder::CyclotronProbeRecorder(std::vector<Config> configs) {
    std::set<std::string> names;
    // Validate all entries before constructing any output sinks.
    for (const auto& config : configs) {
        validate(config);
        if (!names.insert(std::filesystem::absolute(config.output).lexically_normal().string())
                     .second)
            throw std::invalid_argument("PROBE output basenames must be unique: " + config.output);
    }
    // Fresh-run output must not leave an earlier hit file behind when this
    // run has no hits. Match the passive recorder's refusal to reuse outputs.
    int existing = 0;
    if (ippl::Comm->rank() == 0)
        for (const auto& config : configs)
            for (const auto* suffix : {".h5", ".loss", ".hist", ".peaks"}) {
                std::error_code error;
                existing |= std::filesystem::exists(config.output + suffix, error) || bool(error);
            }
    MPI_Bcast(&existing, 1, MPI_INT, 0, ippl::Comm->getCommunicator());
    if (existing)
        throw std::runtime_error(
                "PROBE output already exists or cannot be inspected; use a fresh basename");
    for (auto& config : configs) {
        Diagnostic diagnostic;
        diagnostic.config             = std::move(config);
        const auto [minimum, maximum] = radialBounds(diagnostic.config.segment);
        diagnostic.minimum            = minimum;
        diagnostic.maximum            = maximum;
        diagnostic.histogram.resize(
                static_cast<std::size_t>(std::ceil((maximum - minimum) / diagnostic.config.step))
                + 1);
        diagnostic.sink =
                std::make_unique<LossDataSink>(diagnostic.config.output, !Options::asciidump);
        OpalData::getInstance()->checkAndAddOutputFileName(diagnostic.config.output + ".hist");
        OpalData::getInstance()->checkAndAddOutputFileName(diagnostic.config.output + ".peaks");
        diagnostics_m.push_back(std::move(diagnostic));
    }
}

CyclotronProbeRecorder::~CyclotronProbeRecorder() = default;

void CyclotronProbeRecorder::initialize(
        const std::vector<std::shared_ptr<ParticleContainer_t>>& containers) {
    if (containers.size() > static_cast<std::size_t>(std::numeric_limits<short>::max()))
        throw std::invalid_argument("Too many PROBE particle containers for bunchNumber");
    for (const auto& pc : containers)
        for (std::size_t j = 0; j < diagnostics_m.size(); ++j)
            Kokkos::deep_copy(pc->probeHistory(j).getView(), cyclotron_probe::State{});
}

void CyclotronProbeRecorder::observe(
        const std::shared_ptr<ParticleContainer_t>& pc, double time, short bunch, int turn) {
    if (turn < 0) throw std::invalid_argument("Negative PROBE reference turn");
    const auto positions = pc->R.getView();
    const auto momenta   = pc->P.getView();
    const auto ids       = pc->ID.getView();
    const auto charge    = pc->getQView();
    const auto mass      = pc->getMView();
    const bool sharedQM = pc->getQMStorageMode() == ParticleContainer_t::QMStorageMode::SingleValue;
    const auto origin   = pc->getToLabTrafo().getOrigin();
    const auto rotation = pc->getToLabTrafo().getRotationMatrix();
    if (hits_m.extent(0) < pc->getLocalNum())
        hits_m = Kokkos::View<Hit*>("PROBE compact hits", pc->getLocalNum());
    for (std::size_t j = 0; j < diagnostics_m.size(); ++j) {
        auto& diagnostic   = diagnostics_m[j];
        const auto segment = diagnostic.config.segment;
        const auto states  = pc->probeHistory(j).getView();
        const auto hits    = hits_m;
        const auto count   = count_m;
        Kokkos::deep_copy(count, std::size_t(0));
        int invalid = 0;
        Kokkos::parallel_reduce(
                "PROBE accepted endpoints", pc->getLocalNum(),
                KOKKOS_LAMBDA(std::size_t i, int& errors) {
                    cyclotron_probe::Endpoint endpoint;
                    const cyclotron_probe::Vector offset = positions(i) - origin;
                    endpoint.position                    = prod_vector(rotation, offset);
                    endpoint.momentum                    = prod_vector(rotation, momenta(i));
                    endpoint.time                        = time;
                    cyclotron_probe::Sample sample;
                    const auto status =
                            cyclotron_probe::update(segment, states(i), ids(i), endpoint, sample);
                    if (status == cyclotron_probe::Status::Crossing) {
                        const auto index = Kokkos::atomic_fetch_add(&count(), std::size_t(1));
                        hits(index)      = Hit{
                                sample, ids(i), charge(sharedQM ? 0 : i), mass(sharedQM ? 0 : i)};
                    } else if (
                            status != cyclotron_probe::Status::Initialized
                            && status != cyclotron_probe::Status::Advanced)
                        ++errors;
                },
                invalid);
        MPI_Allreduce(MPI_IN_PLACE, &invalid, 1, MPI_INT, MPI_MAX, ippl::Comm->getCommunicator());
        if (invalid)
            throw std::runtime_error("PROBE encountered invalid or decreasing-time endpoint");
        const auto n    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), count)();
        const auto host = Kokkos::create_mirror_view_and_copy(
                Kokkos::HostSpace(), Kokkos::subview(hits, std::make_pair(std::size_t(0), n)));
        for (std::size_t i = 0; i < n; ++i) {
            const auto& hit = host(i);
            const auto& e   = hit.sample.crossing;
            diagnostic.sink->addParticle(
                    OpalParticle(hit.id, e.position, e.momentum, e.time, hit.charge, hit.mass),
                    std::make_pair(turn, bunch));
            const double radius = std::hypot(e.position(0), e.position(2));
            const auto bin      = static_cast<std::size_t>(std::max(
                    0.0, std::floor((radius - diagnostic.minimum) / diagnostic.config.step)));
            ++diagnostic.histogram[std::min(bin, diagnostic.histogram.size() - 1)];
            auto& peak = diagnostic.peaks[turn];
            peak.first += 1000 * radius;
            ++peak.second;
        }
    }
}

void CyclotronProbeRecorder::save() {
    const auto comm = ippl::Comm->getCommunicator();
    for (auto& diagnostic : diagnostics_m) {
        // Basenames are verified fresh at construction. APPEND creates a new
        // file for the first nonempty batch and preserves earlier batches.
        diagnostic.sink->save(1, OpalData::OpenMode::APPEND);
        std::vector<double> histogram(diagnostic.histogram.size());
        MPI_Reduce(
                diagnostic.histogram.data(), histogram.data(), static_cast<int>(histogram.size()),
                MPI_DOUBLE, MPI_SUM, 0, comm);
        int turns = diagnostic.peaks.empty() ? -1 : diagnostic.peaks.rbegin()->first;
        MPI_Allreduce(MPI_IN_PLACE, &turns, 1, MPI_INT, MPI_MAX, comm);
        int failed = 0;
        std::ofstream peaks, hist;
        if (ippl::Comm->rank() == 0) {
            peaks.open(diagnostic.config.output + ".peaks");
            hist.open(diagnostic.config.output + ".hist");
            peaks << std::setprecision(17) << "# Peak Radii (mm)\n";
            hist << std::setprecision(17) << "# Histogram bin counts (min, max, nbins, binsize) "
                 << 1000 * diagnostic.minimum << " mm " << 1000 * diagnostic.maximum << " mm "
                 << histogram.size() << ' ' << 1000 * diagnostic.config.step << " mm\n";
            for (const auto value : histogram)
                hist << value << '\n';
        }
        for (long long turn = 0; turn <= turns; ++turn) {
            const auto it   = diagnostic.peaks.find(static_cast<int>(turn));
            double local[2] = {0, 0}, global[2] = {0, 0};
            if (it != diagnostic.peaks.end()) {
                local[0] = it->second.first;
                local[1] = it->second.second;
            }
            MPI_Reduce(local, global, 2, MPI_DOUBLE, MPI_SUM, 0, comm);
            if (ippl::Comm->rank() == 0 && global[1] > 0) peaks << global[0] / global[1] << '\n';
        }
        if (ippl::Comm->rank() == 0) {
            hist.close();
            peaks.close();
            failed = !hist || !peaks;
        }
        MPI_Bcast(&failed, 1, MPI_INT, 0, comm);
        if (failed) throw std::runtime_error("Cannot write PROBE histogram/peaks files");
    }
}
