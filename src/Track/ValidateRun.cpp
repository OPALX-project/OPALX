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

        /// Value of a string attribute; empty if it has none.
        std::string text(const ObjectDescription& object, const std::string& name) {
            const AttributeEntry* entry = attribute(object, name);
            if (entry == nullptr || !entry->hasValue || !entry->error.empty()) return "";
            const std::string* value = std::get_if<std::string>(&entry->value);
            return value != nullptr ? *value : "";
        }

    }  // namespace

    std::vector<std::string> validate(const RunDescription& run) {
        std::vector<std::string> errors;
        const TrackDescription& track = run.track;
        const ObjectDescription& cmd  = track.run;
        const bool restart            = run.context.restart;

        if (text(cmd, "METHOD").empty()) {
            errors.push_back("The attribute \"METHOD\" isn't set for the \"RUN\" command");
        }

        // Field solver.
        if (text(cmd, "FIELDSOLVER").empty()) {
            errors.push_back("\"FIELDSOLVER\" must be set in \"RUN\" command.");
        } else if (!run.fieldSolver) {
            errors.push_back("FieldSolverCmd \"" + text(cmd, "FIELDSOLVER") + "\" not found.");
        }

        // There is no command that defines a boundary geometry.
        if (const std::string geometry = text(cmd, "BOUNDARYGEOMETRY");
            !geometry.empty() && geometry != "NONE") {
            errors.push_back(
                    "BOUNDARYGEOMETRY \"" + geometry
                    + "\" is not supported; OPALX has no boundary geometries. Use NONE.");
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

        // The distribution of every source; a restart allows only one-shot sources at t = 0 and no
        // stochastic processes.
        for (std::size_t i = 0; roots && i < beams.size(); ++i) {
            if ((*roots)[i] == i) {
                for (const SourceDescription& source : beams[i]->sources) {
                    if (!source.distributionFound) {
                        errors.push_back(
                                "Distribution \"" + source.distribution
                                + "\" not found; emission source \"" + source.name
                                + "\" names it as its DISTRIBUTION.");
                    } else if (
                            restart
                            && (source.emitted || source.t0 > 0.0
                                || source.distributionType == "EMITTEDFROMFILE")) {
                        errors.push_back(
                                "Checkpoint restart does not yet support time-dependent or "
                                "delayed emission sources (beam '"
                                + names[i] + "', distribution '" + source.distribution + "').");
                    }
                }
            }
            if (restart && !beams[i]->globalProcesses.empty()) {
                errors.push_back(
                        "Checkpoint restart does not yet support stochastic global processes "
                        "such as DECAY (beam '"
                        + names[i] + "').");
            }
        }

        // Beam stops (TRACK ZSTOP for the first beam, BEAM ZSTOP / UNTILMAXSTEPS for others).
        const auto quoted = [&](std::size_t i) {
            return "\"" + names[i] + "\"";
        };
        for (std::size_t i = 0; i < beams.size(); ++i) {
            const BeamDescription& beam = *beams[i];
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
            const bool zStop = beam.zStop.has_value() || (i == 0 && track.zStopGiven);
            if (zStop && beam.untilMaxSteps) {
                errors.push_back(
                        "Beam " + quoted(i) + " has both a ZSTOP" + (i == 0 ? " (from TRACK)" : "")
                        + " and UNTILMAXSTEPS; give only one of them.");
            }
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
