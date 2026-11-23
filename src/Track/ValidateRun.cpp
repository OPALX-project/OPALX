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
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "AbsBeamline/ElementBase.h"
#include "BuildInfo.h"
#include "Physics/Physics.h"
#include "Track/DecayChain.h"
#include "Utilities/OpalException.h"
#include "changes.h"

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

        // An input stating an older OPTION, VERSION is refused while a later release changed the
        // meaning of some input; without OPTION, VERSION it states this release.
        const int currentVersion =
                (buildinfo::version_major * 100 + buildinfo::version_minor) * 100;
        if (run.options.version < currentVersion) {
            const unsigned int fileVersion = run.options.version / 100;
            std::string changes;
            for (const auto& [version, change] : Versions::changes) {
                if (version > fileVersion) changes += change + "\n";
            }
            if (!changes.empty()) {
                errors.push_back(
                        "Version mismatch: the input states OPTION, VERSION = "
                        + std::to_string(run.options.version)
                        + ", and later releases changed the meaning of the input:\n" + changes
                        + "Make sure you understand these changes and adjust the input file "
                          "accordingly. Then add\nOPTION, VERSION = "
                        + std::to_string(currentVersion) + ";\nto the input file.");
            }
        }

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

        // Every MONITOR of the line writes its own loss file, as HDF5 or as ASCII.
        std::map<std::string, int> monitorFiles;
        for (const ElementDescription& element : run.line.elements) {
            if (element.elementType != ElementBase::getTypeString(ElementType::MONITOR)) continue;
            if (monitorFiles.empty() && !run.options.asciidump && !run.options.enableHDF5) {
                errors.push_back(
                        "You must select an OPTION to save Loss data files\n"
                        "Please, choose 'ENABLEHDF5=TRUE' or 'ASCIIDUMP=TRUE'");
            }
            const std::string file =
                    ElementBase::outputFileName(element.object.name, text(element.object, "OUTFN"));
            if (++monitorFiles[file] == 2) {
                errors.push_back("Duplicate file name for output, '" + file + "', detected");
            }
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

        // A fresh run samples the distributions: generated ones take the reference momentum from
        // GAMMA, ENERGY or PC, and REFDIRECTION applies to file distributions only.
        for (std::size_t i = 0; !restart && roots && i < beams.size(); ++i) {
            bool fileSource = false, generatedSource = false;
            if ((*roots)[i] == i) {
                for (const SourceDescription& source : beams[i]->sources) {
                    if (source.distributionType == "FROMFILE"
                        || source.distributionType == "EMITTEDFROMFILE") {
                        fileSource = true;
                    } else if (!source.distributionType.empty()) {
                        generatedSource = true;
                    }
                }
            }
            if (generatedSource && !beams[i]->hasExplicitEnergy) {
                errors.push_back(
                        "The energy of beam \"" + names[i]
                        + "\" hasn't been set. Set either \"GAMMA\", \"ENERGY\" or \"PC\" on "
                          "the BEAM command.");
            }
            if (beams[i]->hasReferenceDirection && !fileSource) {
                errors.push_back(
                        "\"REFDIRECTION\" on beam \"" + names[i]
                        + "\" needs a FROMFILE or EMITTEDFROMFILE distribution.");
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
