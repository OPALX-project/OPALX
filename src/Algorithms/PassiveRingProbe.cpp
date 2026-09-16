// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/PassiveRingProbe.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iterator>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <mpi.h>
#include "Ippl.h"
#include "PartBunch/ParticleContainer.hpp"

namespace {
void broadcastString(std::string& value) {
    int length = static_cast<int>(value.size());
    const auto communicator = ippl::Comm->getCommunicator();
    MPI_Bcast(&length, 1, MPI_INT, 0, communicator);
    value.resize(length);
    MPI_Bcast(value.data(), length, MPI_CHAR, 0, communicator);
}

void collectiveError(std::string error) {
    broadcastString(error);
    if (!error.empty()) throw std::runtime_error("PassiveRingProbe: " + error);
}

std::int64_t parseId(const std::string& token) {
    std::int64_t value = -1;
    const auto result = std::from_chars(token.data(), token.data() + token.size(), value);
    if (token.empty() || token.front() < '0' || token.front() > '9'
        || result.ec != std::errc{} || result.ptr != token.data() + token.size() || value < 0)
        throw std::invalid_argument("Passive probe IDs must be nonnegative signed 64-bit integers");
    return value;
}

void validateConfig(PassiveRingProbe::Config& config) {
    if (config.output.empty() || config.ids.empty() || config.planes.empty())
        throw std::invalid_argument("Passive probes require OUTPUT, IDS and at least one PLANE");
    if (config.ids.size() > static_cast<std::size_t>(std::numeric_limits<int>::max() / 6)
        || config.planes.size() > std::numeric_limits<std::size_t>::max() / config.ids.size())
        throw std::invalid_argument("Passive probe configuration exceeds collective array limits");
    std::sort(config.ids.begin(), config.ids.end());
    if (config.ids.front() < 0
        || std::adjacent_find(config.ids.begin(), config.ids.end()) != config.ids.end())
        throw std::invalid_argument("Passive probe particle IDs must be unique and nonnegative");
    std::set<std::int64_t> planeIds;
    for (const auto& named : config.planes) {
        passive_probe::State state;
        passive_probe::Sample sample;
        if (named.id < 0 || !planeIds.insert(named.id).second
            || passive_probe::update(named.plane, state, passive_probe::Endpoint{}, sample)
                       != passive_probe::Status::Initialized)
            throw std::invalid_argument("Invalid or duplicate passive probe plane");
    }
}
} // namespace

PassiveRingProbe::Config PassiveRingProbe::parseConfig(const std::string& text) {
    Config config;
    bool hasOutput = false, hasIds = false;
    std::istringstream input(text);
    std::string line;
    while (std::getline(input, line)) {
        std::istringstream fields(line);
        std::string command;
        if (!(fields >> command) || command.front() == '#') continue;
        if (command == "OUTPUT") {
            if (hasOutput || !(fields >> std::quoted(config.output)))
                throw std::invalid_argument("Invalid or repeated passive OUTPUT");
            hasOutput = true;
        } else if (command == "IDS") {
            if (hasIds) throw std::invalid_argument("Repeated passive IDS");
            hasIds = true;
            std::string token;
            while (fields >> token) config.ids.push_back(parseId(token));
        } else if (command == "PLANE") {
            NamedPlane named;
            std::string id;
            if (!(fields >> id)) throw std::invalid_argument("Missing passive plane ID");
            named.id = parseId(id);
            for (unsigned d = 0; d < 3; ++d)
                if (!(fields >> named.plane.origin(d)))
                    throw std::invalid_argument("Invalid passive plane origin");
            for (unsigned d = 0; d < 3; ++d)
                if (!(fields >> named.plane.normal(d)))
                    throw std::invalid_argument("Invalid passive plane normal");
            fields >> std::ws;
            if (!fields.eof() && !(fields >> named.plane.tolerance))
                throw std::invalid_argument("Invalid passive plane tolerance");
            config.planes.push_back(named);
        } else throw std::invalid_argument("Unknown passive probe directive: " + command);
        std::string extra;
        if (fields >> extra) throw std::invalid_argument("Trailing passive probe configuration text");
    }
    validateConfig(config);
    return config;
}

std::unique_ptr<PassiveRingProbe> PassiveRingProbe::fromEnvironment(bool eligible) {
    std::string source, error;
    int enabled = 0;
    if (ippl::Comm->rank() == 0) {
        const char* path = std::getenv("OPALX_TEST_PASSIVE_RING_PROBES");
        enabled = path != nullptr;
        if (enabled) {
            try {
                if (!*path) throw std::invalid_argument("Empty OPALX_TEST_PASSIVE_RING_PROBES path");
                std::ifstream file(path);
                if (!file) throw std::runtime_error("Cannot open passive probe configuration");
                source.assign(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
                if (file.bad() || source.size() > 1024 * 1024)
                    throw std::runtime_error("Unreadable or excessive passive probe configuration");
                parseConfig(source);
            } catch (const std::exception& exception) { error = exception.what(); }
        }
    }
    const auto communicator = ippl::Comm->getCommunicator();
    MPI_Bcast(&enabled, 1, MPI_INT, 0, communicator);
    if (!enabled) return nullptr;
    collectiveError(error);
    int allEligible = eligible ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &allEligible, 1, MPI_INT, MPI_MIN, communicator);
    if (!allEligible)
        throw std::runtime_error("PassiveRingProbe requires a fresh aligned analytic bare ring");
    broadcastString(source);
    return std::make_unique<PassiveRingProbe>(parseConfig(source));
}

PassiveRingProbe::PassiveRingProbe(Config config) : config_m(std::move(config)) {
    validateConfig(config_m);
    const auto count = config_m.ids.size();
    ids_m = Kokkos::View<std::int64_t*>("passive stable IDs", count);
    coordinates_m = Kokkos::View<double*>("passive selected endpoints", 6 * count);
    ownership_m = Kokkos::View<int*>("passive selected ownership", count);
    auto idsHost = Kokkos::create_mirror_view(ids_m);
    for (std::size_t i = 0; i < count; ++i) idsHost(i) = config_m.ids[i];
    Kokkos::deep_copy(ids_m, idsHost);
    std::string error;
    if (ippl::Comm->rank() == 0) {
        try {
            if (std::filesystem::exists(config_m.output))
                throw std::runtime_error("Passive probe output already exists: " + config_m.output);
            output_m.open(config_m.output);
            if (!output_m) throw std::runtime_error("Cannot create passive probe output: " + config_m.output);
            output_m << "particle_id,plane_id,turn,time_s,x_m,y_m,z_m,px_mc,py_mc,pz_mc,status,"
                        "bracket_start_time_s,bracket_end_time_s,fraction\n" << std::setprecision(17);
            output_m.flush();
            if (!output_m) throw std::runtime_error("Cannot write passive probe header");
            histories_m.resize(count * config_m.planes.size());
            missing_m.resize(count, false);
        } catch (const std::exception& exception) { error = exception.what(); }
    }
    collectiveError(error);
}

void PassiveRingProbe::observe(
        const std::shared_ptr<ParticleContainer_t>& particles, double acceptedTime) {
    const auto communicator = ippl::Comm->getCommunicator();
    int invalid = !particles || !std::isfinite(acceptedTime);
    MPI_Allreduce(MPI_IN_PLACE, &invalid, 1, MPI_INT, MPI_MAX, communicator);
    if (invalid) throw std::runtime_error("PassiveRingProbe received invalid container or time");
    double minimumTime = acceptedTime, maximumTime = acceptedTime;
    MPI_Allreduce(MPI_IN_PLACE, &minimumTime, 1, MPI_DOUBLE, MPI_MIN, communicator);
    MPI_Allreduce(MPI_IN_PLACE, &maximumTime, 1, MPI_DOUBLE, MPI_MAX, communicator);
    if (minimumTime != maximumTime)
        throw std::runtime_error("PassiveRingProbe requires synchronized accepted times");

    Kokkos::deep_copy(coordinates_m, 0.0);
    Kokkos::deep_copy(ownership_m, 0);
    const auto selectedIds = ids_m;
    const auto coordinates = coordinates_m;
    const auto ownership = ownership_m;
    const auto ids = particles->ID.getView();
    const auto positions = particles->R.getView();
    const auto momenta = particles->P.getView();
    const auto origin = particles->getToLabTrafo().getOrigin();
    const auto rotation = particles->getToLabTrafo().getRotationMatrix();
    const std::size_t count = config_m.ids.size();
    Kokkos::parallel_for("passive selected endpoint gather", particles->getLocalNum(),
        KOKKOS_LAMBDA(std::size_t i) {
            std::size_t first = 0, last = count;
            while (first < last) {
                const std::size_t middle = first + (last - first) / 2;
                if (selectedIds(middle) < ids(i)) first = middle + 1;
                else last = middle;
            }
            if (first == count || selectedIds(first) != ids(i)) return;
            // Only the first local owner writes; duplicates are errors after the
            // reduction, without a concurrent-write data race in the meantime.
            if (Kokkos::atomic_fetch_add(&ownership(first), 1) != 0) return;
            const passive_probe::Vector offset = positions(i) - origin;
            const auto position = prod_vector(rotation, offset);
            const auto momentum = prod_vector(rotation, momenta(i));
            for (unsigned d = 0; d < 3; ++d) {
                coordinates(6 * first + d) = position(d);
                coordinates(6 * first + 3 + d) = momentum(d);
            }
        });
    const auto localCoordinates = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), coordinates_m);
    const auto localOwnership = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ownership_m);
    std::vector<double> globalCoordinates(6 * count);
    std::vector<int> globalOwnership(count);
    MPI_Allreduce(localCoordinates.data(), globalCoordinates.data(), static_cast<int>(6 * count),
                  MPI_DOUBLE, MPI_SUM, communicator);
    MPI_Allreduce(localOwnership.data(), globalOwnership.data(), static_cast<int>(count),
                  MPI_INT, MPI_SUM, communicator);
    std::string error;
    if (ippl::Comm->rank() == 0) {
        try {
            for (std::size_t i = 0; i < count; ++i) {
                if (globalOwnership[i] > 1)
                    throw std::runtime_error("Duplicate selected particle ID " + std::to_string(config_m.ids[i]));
                passive_probe::Endpoint endpoint;
                endpoint.time = acceptedTime;
                for (unsigned d = 0; d < 3; ++d) {
                    endpoint.position(d) = globalCoordinates[6 * i + d];
                    endpoint.momentum(d) = globalCoordinates[6 * i + 3 + d];
                }
                for (std::size_t j = 0; j < config_m.planes.size(); ++j) {
                    auto& history = histories_m[i * config_m.planes.size() + j];
                    const auto& named = config_m.planes[j];
                    if (globalOwnership[i] == 0) {
                        if (!missing_m[i])
                            output_m << config_m.ids[i] << ',' << named.id << ',' << history.turns
                                     << ',' << acceptedTime << ",nan,nan,nan,nan,nan,nan,missing,nan,nan,nan\n";
                        history.initialized = false;
                        history.armed = false;
                        continue;
                    }
                    passive_probe::Sample sample;
                    const double beforeTime = history.previous.time;
                    const auto status = passive_probe::update(named.plane, history, endpoint, sample);
                    if (status == passive_probe::Status::Crossing) {
                        output_m << config_m.ids[i] << ',' << named.id << ',' << sample.turn << ','
                                 << sample.crossing.time;
                        for (unsigned d = 0; d < 3; ++d) output_m << ',' << sample.crossing.position(d);
                        for (unsigned d = 0; d < 3; ++d) output_m << ',' << sample.crossing.momentum(d);
                        output_m << ",crossing," << beforeTime << ',' << acceptedTime << ','
                                 << sample.fraction << '\n';
                    } else if (status != passive_probe::Status::Initialized
                               && status != passive_probe::Status::Advanced) {
                        const char* statusName = status == passive_probe::Status::InvalidInput
                                ? "InvalidInput" : status == passive_probe::Status::DecreasingTime
                                ? "DecreasingTime" : status == passive_probe::Status::TurnOverflow
                                ? "TurnOverflow" : "UnexpectedStatus";
                        std::ostringstream message;
                        message << std::setprecision(std::numeric_limits<double>::max_digits10)
                                << "Invalid passive observation: particle ID " << config_m.ids[i]
                                << ", plane ID " << named.id << ", status=" << statusName
                                << " (" << static_cast<int>(status) << ')'
                                << ", previous_time_s=" << beforeTime
                                << ", current_time_s=" << acceptedTime;
                        throw std::runtime_error(message.str());
                    }
                }
                missing_m[i] = globalOwnership[i] == 0;
            }
            output_m.flush();
            if (!output_m) throw std::runtime_error("Cannot write passive probe samples");
        } catch (const std::exception& exception) { error = exception.what(); }
    }
    collectiveError(error);
}
