//
// validate(): every check of a RUN that needs only the input values, in one place.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include "Track/ValidateRun.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "Physics/Physics.h"
#include "Track/DecayChain.h"
#include "Utilities/OpalException.h"

namespace opalx::run {

    namespace {

        const AttributeEntry* attribute(const ObjectDescription& object, const std::string& name) {
            for (const AttributeEntry& entry : object.attributes) {
                if (entry.name == name) return &entry;
            }
            return nullptr;
        }

        /// True if the input set the attribute (Attribute::defaultUsed() is false).
        bool isSet(const ObjectDescription& object, const std::string& name) {
            const AttributeEntry* entry = attribute(object, name);
            return entry != nullptr && !entry->defaultUsed;
        }

        double real(const ObjectDescription& object, const std::string& name) {
            const AttributeEntry* entry = attribute(object, name);
            if (entry == nullptr || !entry->hasValue || !entry->error.empty()) return 0.0;
            const double* value = std::get_if<double>(&entry->value);
            return value != nullptr ? *value : 0.0;
        }

        bool logical(const ObjectDescription& object, const std::string& name) {
            const AttributeEntry* entry = attribute(object, name);
            if (entry == nullptr || !entry->hasValue || !entry->error.empty()) return false;
            const bool* value = std::get_if<bool>(&entry->value);
            return value != nullptr && *value;
        }

        std::string text(const ObjectDescription& object, const std::string& name) {
            const AttributeEntry* entry = attribute(object, name);
            if (entry == nullptr || !entry->hasValue || !entry->error.empty()) return "";
            const std::string* value = std::get_if<std::string>(&entry->value);
            return value != nullptr ? *value : "";
        }

        /// The element types ParallelTracker tracks with localized TURNS
        /// (ParallelTracker::execute(), before this check moved here).
        bool supportedForTurns(const ElementDescription& element) {
            const std::string& type = element.elementType;
            return type == "Drift" || type == "Marker" || type == "Monitor" || type == "Multipole"
                   || type == "MULTIPOLET" || type == "SBEND" || type == "RBEND"
                   || type == "SOLENOID" || type == "CyclotronSector" || type == "ScalingFFAMagnet"
                   || element.object.keyword == "VARIABLE_RF_CAVITY";
        }

        /// True for elements that reach OpalBeamline: MARKER is dropped, PROBE only configures
        /// a recorder, LASER throws, and the elements of type ANY other than
        /// VARIABLE_RF_CAVITY (ENGE, TANH, ASYMMETRIC_ENGE, VERTICALFFAMAGNET) have no tracked
        /// copy.
        bool isTracked(const ElementDescription& element) {
            const std::string& type = element.elementType;
            if (type == "Marker" || type == "Probe" || type == "Laser") return false;
            if (type == "Any") return element.object.keyword == "VARIABLE_RF_CAVITY";
            return true;
        }

        /// An RF cavity tracked as a single cyclotron gap (OpalCavity::build() sets the profile
        /// for TYPE = SINGLEGAP with a FMAPFN).
        bool isCyclotronGap(const ElementDescription& element) {
            return element.elementType == "RFCavity" && text(element.object, "TYPE") == "SINGLEGAP"
                   && !text(element.object, "FMAPFN").empty();
        }

    }  // namespace

    std::vector<std::string> validate(const RunDescription& run) {
        std::vector<std::string> errors;
        const TrackDescription& track = run.track;
        const ObjectDescription& cmd  = track.run;
        const bool restart            = run.context.restart;
        const bool isRing             = run.line.topology == "RING";
        const bool spectralTunes      = logical(cmd, "SPECTRALTUNES");
        const double kineticStop      = track.kineticEnergyStopGeV;

        // Field solver.
        if (attribute(cmd, "FIELDSOLVER") == nullptr) {
            errors.push_back("\"FIELDSOLVER\" must be set in \"RUN\" command.");
        } else if (!run.fieldSolver) {
            errors.push_back("FieldSolverCmd \"" + text(cmd, "FIELDSOLVER") + "\" not found.");
        }

        // Beams, in TRACK BEAMS order.
        const std::vector<std::string>& names = track.beamNames;
        std::vector<const BeamDescription*> beams;
        if (names.empty()) {
            errors.push_back("No beam specified: set TRACK::BEAM or TRACK::BEAMS.");
        }
        for (const std::string& name : names) {
            const BeamDescription* found = nullptr;
            for (const BeamDescription& beam : run.beams) {
                if (beam.object.name == name) found = &beam;
            }
            if (name.empty()) {
                errors.push_back("Empty beam name in resolved beam list.");
            } else if (found == nullptr) {
                errors.push_back("Beam \"" + name + "\" not found.");
            }
            beams.push_back(found);
        }
        const bool beamsFound =
                !beams.empty() && std::find(beams.begin(), beams.end(), nullptr) == beams.end();
        if (!beamsFound) return errors;

        for (const BeamDescription* beam : beams) {
            if (beam->isPhoton) {
                errors.push_back(
                        "TRACK does not support BEAM, PARTICLE=PHOTON yet. Photon beams may be "
                        "defined for future OPALX features, but they are currently rejected "
                        "during tracking.");
            }
        }

        // Decay chains: a beam named as another beam's DAUGHTERBEAM follows the first beam of
        // its chain and has no momentum or emission sources of its own.
        std::optional<std::vector<std::size_t>> roots;
        {
            std::vector<std::string> daughters;
            std::vector<bool> hasDecay;
            for (const BeamDescription* beam : beams) {
                daughters.push_back(beam->daughterBeam);
                hasDecay.push_back(
                        std::find(
                                beam->globalProcesses.begin(), beam->globalProcesses.end(), "DECAY")
                        != beam->globalProcesses.end());
            }
            try {
                roots = resolveDecayRoots(names, daughters, hasDecay);
            } catch (const OpalException& error) {
                errors.push_back(error.what());
            }
        }
        if (roots && (*roots)[0] != 0) {
            errors.push_back(
                    "The first beam in BEAMS (\"" + names[0]
                    + "\") receives decay daughters (DAUGHTERBEAM of another beam). List beam \""
                    + names[(*roots)[0]] + "\" first.");
        }
        for (std::size_t i = 0; roots && i < beams.size(); ++i) {
            if ((*roots)[i] == i) continue;
            if (beams[i]->hasExplicitEnergy) {
                errors.push_back(
                        "\"GAMMA\", \"ENERGY\" and \"PC\" are not allowed on beam \"" + names[i]
                        + "\": it receives decay daughters (DAUGHTERBEAM of another beam) and "
                          "uses the reference particle of beam \""
                        + names[(*roots)[i]] + "\".");
            }
            if (beams[i]->hasEmissionSources) {
                errors.push_back(
                        "\"SOURCES\" is not allowed on beam \"" + names[i]
                        + "\": it receives particles only from decays (DAUGHTERBEAM of another "
                          "beam).");
            }
        }

        // Emission sources of every beam that starts its own decay chain.
        for (std::size_t i = 0; roots && i < beams.size(); ++i) {
            if ((*roots)[i] != i) continue;
            const BeamDescription& beam = *beams[i];
            if (!beam.emissionSourceListError.empty()) {
                errors.push_back(beam.emissionSourceListError);
            } else if (!beam.emissionSourceListFound) {
                errors.push_back(
                        "EmissionSourceList \"" + beam.emissionSourceList + "\" not found.");
            } else if (beam.sources.empty()) {
                errors.push_back(
                        "Emission sources list for beam '" + names[i]
                        + "' must contain at least one EMISSIONSOURCE.");
            }
        }

        // Restart: only one-shot sources at t = 0 and no stochastic processes.
        if (restart) {
            for (std::size_t i = 0; roots && i < beams.size(); ++i) {
                if ((*roots)[i] == i) {
                    for (const SourceDescription& source : beams[i]->sources) {
                        if (!source.distributionFound) {
                            errors.push_back(
                                    "Distribution \"" + source.distribution + "\" not found.");
                        } else if (
                                source.emitted || source.t0 > 0.0
                                || source.distributionType == "EMITTEDFROMFILE") {
                            errors.push_back(
                                    "Checkpoint restart does not yet support time-dependent or "
                                    "delayed emission sources (beam '"
                                    + names[i] + "', distribution '" + source.distribution + "').");
                        }
                    }
                }
                if (!beams[i]->globalProcesses.empty()) {
                    errors.push_back(
                            "Checkpoint restart does not yet support stochastic global processes "
                            "such as DECAY (beam '"
                            + names[i] + "').");
                }
            }
        }

        // RING length, EKINSTOP and TURNS.
        const double ringPeriod = isRing ? run.line.length : 0.0;
        if (isRing && (!(ringPeriod > 0.0) || !std::isfinite(ringPeriod))) {
            errors.push_back("A RING requires a positive circumference.");
        }
        const bool oneSegment =
                track.dT.size() == 1 && track.maxSteps.size() == 1 && track.zStop.size() == 1;
        if (kineticStop > 0) {
            if (!isRing || restart || isSet(cmd, "TURNS") || !oneSegment
                || !(track.dT.front() > 0)) {
                errors.push_back(
                        "EKINSTOP requires a non-restarted RING, one positive DT segment and no "
                        "explicit TURNS.");
            }
        }
        const bool localizedTurns = isSet(cmd, "TURNS") && isRing && !spectralTunes;
        if (localizedTurns) {
            const auto turnsError = [&]() -> std::optional<std::string> {
                const double requested = real(cmd, "TURNS");
                const double rounded   = std::round(requested);
                if (!std::isfinite(requested)
                    || requested >= static_cast<double>(
                               std::numeric_limits<unsigned long long>::max())
                    || requested < 1.0 || std::abs(requested - rounded) > 1.0e-12) {
                    return "TURNS must be a positive integer for a RING.";
                }
                if (!oneSegment) {
                    return "RING tracking with explicit TURNS currently requires one "
                           "DT/MAXSTEPS/ZSTOP segment.";
                }
                if (beams.size() != 1) return "Localized TURNS requires one beam.";
                if (restart) {
                    return "Explicit RING TURNS restart requires persisted return-plane counters "
                           "and is not supported yet.";
                }
                const double distancePerStep =
                        Physics::c * std::abs(track.dT.front()) * beams.front()->reference.beta;
                if (!(distancePerStep > 0.0)) {
                    return "A RING tracked with TURNS requires positive circumference, DT, and "
                           "reference velocity.";
                }
                const double safetySteps =
                        2.0 * std::ceil(rounded * ringPeriod / distancePerStep) + 100.0;
                if (!std::isfinite(safetySteps)
                    || safetySteps >= static_cast<double>(
                               std::numeric_limits<unsigned long long>::max())) {
                    return "RING TURNS step budget overflows.";
                }
                return std::nullopt;
            }();
            if (turnsError) errors.push_back(*turnsError);
        }

        // INITIALORBIT (the orbit itself is checked by TrackRun against the COF result).
        if (track.hasInitialOrbit && (beams.size() != 1 || restart || spectralTunes)) {
            errors.push_back(
                    "INITIALORBIT: Requires one fresh beam and ordinary particle tracking.");
        }

        // SPECTRALTUNES.
        if (spectralTunes) {
            const double turns                 = real(cmd, "TURNS");
            const double sample                = real(cmd, "TUNESAMPLE");
            const BeamDescription& beam        = *beams.front();
            const std::vector<double>& initial = beam.tuneInitial;
            if (!isRing || restart || kineticStop > 0 || beams.size() != 1
                || track.fieldSolverType != "NONE" || beam.particleName != "PROTON"
                || track.dT.size() != 1 || initial.empty() || initial.size() % 3
                || !(turns >= 1 && turns <= 100000 && turns == std::floor(turns))
                || !(sample >= 1 && sample <= 100000 && sample == std::floor(sample))
                || track.stepsPerTurn <= 0 || !(track.dT.front() > 0)
                || !beam.globalProcesses.empty()) {
                errors.push_back(
                        "SPECTRALTUNES requires non-restarted proton RING, one beam/DT, "
                        "FIELDSOLVER=NONE, TUNEINITIAL triples and positive integer "
                        "TURNS/TUNESAMPLE.");
            } else if (
                    std::size_t(static_cast<unsigned>(turns)) * track.stepsPerTurn
                    > track.maxSteps.front()) {
                errors.push_back("Spectral tune steps exceed TRACK MAXSTEPS.");
            }
        }

        // Beam stops (TRACK ZSTOP for the first beam, BEAM ZSTOP / UNTILMAXSTEPS for others).
        const auto quoted = [&](std::size_t i) {
            return "\"" + names[i] + "\"";
        };
        for (std::size_t i = 0; i < beams.size(); ++i) {
            const BeamDescription& beam = *beams[i];
            if ((beam.zStop || beam.untilMaxSteps) && isRing) {
                errors.push_back(
                        "BEAM ZSTOP and UNTILMAXSTEPS apply to LINE tracking only (beam "
                        + quoted(i) + ").");
            }
            if (beam.zStop && i == 0) {
                errors.push_back(
                        "Beam " + quoted(i)
                        + " is the first beam of the TRACK; set its stop with TRACK ZSTOP, not "
                          "BEAM ZSTOP.");
            }
            if (beam.zStop && roots && (*roots)[i] != i) {
                errors.push_back(
                        "Beam " + quoted(i)
                        + " receives decay daughters and does not start at the start of the "
                          "line, so it cannot have a ZSTOP.");
            }
            const bool zStop = beam.zStop.has_value() || (i == 0 && !isRing && track.zStopGiven);
            if (zStop && beam.untilMaxSteps) {
                errors.push_back(
                        "Beam " + quoted(i) + " has both a ZSTOP" + (i == 0 ? " (from TRACK)" : "")
                        + " and UNTILMAXSTEPS; give only one of them.");
            }
        }

        // Checks ParallelTracker made on the element list.
        bool hasProbe = false, hasCyclotronSector = false, hasCyclotronGaps = false;
        for (const ElementDescription& element : run.line.elements) {
            hasProbe |= element.elementType == "Probe";
            hasCyclotronSector |= element.elementType == "CyclotronSector";
            hasCyclotronGaps |= isCyclotronGap(element);
        }
        if (hasProbe) {
            if (restart) {
                errors.push_back(
                        "PROBE: PROBE does not yet support restart or tracking continuation.");
            }
            if ((!hasCyclotronSector && ringPeriod <= 0) || spectralTunes) {
                errors.push_back(
                        "PROBE: PROBE requires cyclotron-sector or RING particle tracking; "
                        "spectral-tune trial tracking is unsupported.");
            }
        }
        if (hasCyclotronGaps) {
            const BeamDescription& beam = *beams.front();
            if (beams.size() != 1) {
                errors.push_back("SINGLEGAP currently requires one container.");
            } else if (
                    beam.hasPolarization || beam.reference.charge != 1 || restart
                    || std::abs(beam.reference.mass / (Physics::m_p * 1e9) - 1) > 1e-12) {
                errors.push_back(
                        "SINGLEGAP requires an unpolarized positive proton and does not support "
                        "restart.");
            }
        }
        if (localizedTurns && beams.size() == 1) {
            for (const ElementDescription& element : run.line.elements) {
                if (isTracked(element) && !supportedForTurns(element)) {
                    errors.push_back(
                            "Localized TURNS requires supported magnetic elements or "
                            "VARIABLE_RF_CAVITY; unsupported element "
                            + element.object.name);
                }
            }
        }
        if (kineticStop > 0 && (!hasCyclotronGaps || beams.size() != 1)) {
            errors.push_back("EKINSTOP requires the single-container cyclotron gap path.");
        }
        return errors;
    }

    std::string joinErrors(const std::vector<std::string>& errors) {
        std::string joined;
        for (const std::string& error : errors) {
            if (!joined.empty()) joined += "\n";
            joined += error;
        }
        return joined;
    }

}  // namespace opalx::run
