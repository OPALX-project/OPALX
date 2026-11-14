//
// freeze(): copy the values one RUN command tracks with out of the parsed objects.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include "Track/FreezeRun.h"

#include <algorithm>
#include <exception>
#include <string>

#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/Object.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartData.h"
#include "Attributes/Attributes.h"
#include "Distribution/Distribution.h"
#include "Lines/EmissionSourceList.h"
#include "Lines/Line.h"
#include "Lines/Ring.h"
#include "Structure/Beam.h"
#include "Structure/BinningCmd.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/EmissionSource.h"
#include "Structure/FieldSolverCmd.h"
#include "Track/BuildElements.h"
#include "Track/Track.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

namespace opalx::run {

    namespace {

        AttributeValue evaluate(const Attribute& attribute) {
            const std::string& type = attribute.getType();
            if (type == "real") return Attributes::getReal(attribute);
            if (type == "logical") return Attributes::getBool(attribute);
            if (type == "string") return Attributes::getString(attribute);
            if (type == "real array") return Attributes::getRealArray(attribute);
            if (type == "logical array") return Attributes::getBoolArray(attribute);
            if (type == "string array") return Attributes::getStringArray(attribute);
            return attribute.getImage();
        }

        AttributeEntry describeAttribute(const Attribute& attribute) {
            AttributeEntry entry;
            entry.name        = attribute.getName();
            entry.type        = attribute.getType();
            entry.defaultUsed = attribute.defaultUsed();
            entry.hasValue    = attribute.isBaseAllocated();
            if (!entry.hasValue) return entry;
            try {
                entry.value = evaluate(attribute);
            } catch (const OpalException& error) {
                entry.error = error.what();
            } catch (const std::exception& error) {
                entry.error = error.what();
            }
            return entry;
        }

        /// Look up a parsed object by name; nullptr if there is none or it has another type.
        template <class T>
        const T* findAs(const std::string& name) {
            if (name.empty()) return nullptr;
            return dynamic_cast<const T*>(OpalData::getInstance()->find(name));
        }

        /**
         * Append the members of @p line to @p out, in order. A nested LINE is entered; its
         * members are reflected when the nested LINE is, unless they are reflected themselves.
         * Reflection does not change the order. @p path holds the lines being entered, to
         * refuse a LINE that contains itself.
         */
        void addMembers(
                const Line& line, bool reflected, std::vector<ElementDescription>& out,
                std::vector<std::string>& path) {
            for (const Line::Member& member : line.getMembers()) {
                const bool memberReflected = reflected != member.reflected;
                if (const auto* nested = findAs<Line>(member.name)) {
                    if (std::find(path.begin(), path.end(), member.name) != path.end()) {
                        throw OpalException(
                                "opalx::run::describeLine",
                                "LINE \"" + member.name + "\" contains itself.");
                    }
                    path.push_back(member.name);
                    addMembers(*nested, memberReflected, out, path);
                    path.pop_back();
                    continue;
                }
                ElementDescription description;
                description.occurrence = out.size();
                description.reflected  = memberReflected;
                if (const auto* input = findAs<Element>(member.name)) {
                    description.object      = describeObject(*input);
                    description.elementType = elementTypeOf(description.object.keyword);
                } else {
                    description.object.name = member.name;
                }
                out.push_back(std::move(description));
            }
        }

        OptionValues currentOptions() {
            OptionValues o;
            o.echo                    = Options::echo;
            o.info                    = Options::info;
            o.mtrace                  = Options::mtrace;
            o.warn                    = Options::warn;
            o.fft2d5BatchSize         = Options::fft2d5BatchSize;
            o.seed                    = Options::seed;
            o.psDumpFreq              = Options::psDumpFreq;
            o.statDumpFreq            = Options::statDumpFreq;
            o.checkpointFreq          = Options::checkpointFreq;
            o.stepInfoFreq            = Options::stepInfoFreq;
            o.printRankDistrFreq      = Options::printRankDistrFreq;
            o.useQMAttributes         = Options::useQMAttributes;
            o.repartFreq              = Options::repartFreq;
            o.minBinEmitted           = Options::minBinEmitted;
            o.rhoDump                 = Options::rhoDump;
            o.ebDump                  = Options::ebDump;
            o.rankDump                = Options::rankDump;
            o.autoPhase               = Options::autoPhase;
            o.enableHDF5              = Options::enableHDF5;
            o.enableVTK               = Options::enableVTK;
            o.elementMappingGrid      = Options::elementMappingGrid;
            o.elementMappingReach     = Options::elementMappingReach;
            o.asciidump               = Options::asciidump;
            o.boundpDestroy           = Options::boundpDestroy;
            o.beamHaloBoundary        = Options::beamHaloBoundary;
            o.version                 = Options::version;
            o.haloShift               = Options::haloShift;
            o.computePercentiles      = Options::computePercentiles;
            o.aggressiveStateSync     = Options::aggressiveStateSync;
            o.loadBalancingThreshold  = Options::loadBalancingThreshold;
            o.loadBalancingPairWeight = Options::loadBalancingPairWeight;
            return o;
        }

        std::string openModeName(OpalData::OpenMode mode) {
            switch (mode) {
                case OpalData::OpenMode::WRITE:
                    return "WRITE";
                case OpalData::OpenMode::APPEND:
                    return "APPEND";
                default:
                    return "UNDEFINED";
            }
        }

        /// Run @p read; a parse object that cannot answer (OpalException) leaves the field at its
        /// default and stores the message in @p error, if given. TrackRun reports such errors.
        template <class F>
        void tryRead(F&& read, std::string* error = nullptr) {
            try {
                read();
            } catch (const OpalException& e) {
                if (error != nullptr) *error = e.what();
            }
        }

        /// Attribute of @p object by name, or nullptr.
        const Attribute* attributeNamed(const Object& object, const std::string& name) {
            for (const Attribute& attribute : object.itsAttr) {
                if (attribute && attribute.getName() == name) return &attribute;
            }
            return nullptr;
        }

        SourceDescription describeSource(const EmissionSource& source) {
            SourceDescription description;
            description.name = source.getOpalName();
            tryRead([&] {
                description.distribution = source.getDistributionName();
            });
            tryRead([&] {
                description.t0 = source.getT0();
            });
            if (const auto* dist = findAs<Distribution>(description.distribution)) {
                description.distributionFound = true;
                if (const Attribute* type = attributeNamed(*dist, "TYPE")) {
                    tryRead([&] {
                        description.distributionType = Attributes::getString(*type);
                    });
                }
                if (const Attribute* emitted = attributeNamed(*dist, "EMITTED")) {
                    tryRead([&] {
                        description.emitted = Attributes::getBool(*emitted);
                    });
                }
            }
            return description;
        }

        BeamDescription describeBeam(const Beam& beam) {
            BeamDescription d;
            d.object                  = describeObject(beam);
            const PartData& reference = beam.getReference();
            d.reference.mass          = reference.getM();
            d.reference.charge        = reference.getQ();
            d.reference.momentum      = reference.getP();
            d.reference.gamma         = reference.getGamma();
            d.reference.beta          = reference.getBeta();
            d.chargePerParticle       = beam.getChargePerParticle();
            d.massPerParticle         = beam.getMassPerParticle();
            tryRead(
                    [&] {
                        d.emissionSourceList = beam.getEmissionSourceListName();
                    },
                    &d.emissionSourceListError);

            tryRead([&] {
                d.particleName = beam.getParticleName();
            });
            tryRead([&] {
                d.isPhoton = beam.isPhoton();
            });
            tryRead([&] {
                d.daughterBeam = beam.getDaughterBeamName();
            });
            tryRead([&] {
                d.globalProcesses = beam.getGlobalProcessNames();
            });
            tryRead([&] {
                d.hasExplicitEnergy = beam.hasExplicitEnergy();
            });
            tryRead([&] {
                d.hasEmissionSources = beam.hasEmissionSources();
            });
            tryRead([&] {
                d.hasPolarization = beam.hasPolarization();
            });
            tryRead([&] {
                d.zStop = beam.getZStop();
            });
            tryRead([&] {
                d.untilMaxSteps = beam.getUntilMaxSteps();
            });
            tryRead([&] {
                d.tuneInitial = beam.getTuneInitial();
            });

            if (const auto* list = findAs<EmissionSourceList>(d.emissionSourceList)) {
                d.emissionSourceListFound = true;
                for (const EmissionSource* source : list->fetchSources()) {
                    if (source != nullptr) d.sources.push_back(describeSource(*source));
                }
            }
            return d;
        }

    }  // namespace

    ObjectDescription describeObject(const Object& object) {
        ObjectDescription description;
        description.keyword = object.getKeyword();
        description.name    = object.getOpalName();
        description.attributes.reserve(object.itsAttr.size());
        for (const Attribute& attribute : object.itsAttr) {
            description.attributes.push_back(describeAttribute(attribute));
        }
        return description;
    }

    RunDescription freeze(const Object& runCommand, int runIndex) {
        RunDescription run;
        OpalData* opal = OpalData::getInstance();

        RunContext& context      = run.context;
        context.runIndex         = runIndex;
        context.inputBasename    = opal->getInputBasename();
        context.restart          = opal->inRestartRun();
        context.restartFileName  = context.restart ? opal->getRestartFileName() : "";
        context.restartStep      = opal->getRestartStep();
        context.openMode         = openModeName(opal->getOpenMode());
        context.p0               = opal->getP0();
        context.globalPhaseShift = opal->getGlobalPhaseShift();

        run.options = currentOptions();

        const Track* track = Track::block;
        if (track == nullptr) {
            throw OpalException("opalx::run::freeze", "RUN is only valid inside TRACK.");
        }
        TrackDescription& trackDescription    = run.track;
        trackDescription.beamNames            = track->beamNames_m;
        trackDescription.dT                   = track->dT;
        trackDescription.maxSteps             = track->localTimeSteps;
        trackDescription.zStop                = track->zstop;
        trackDescription.zStopGiven           = track->zStopGiven;
        trackDescription.zStart               = track->zstart;
        trackDescription.t0                   = track->t0_m;
        trackDescription.stepsPerTurn         = track->stepsPerTurn;
        trackDescription.kineticEnergyStopGeV = track->kineticEnergyStopGeV;
        trackDescription.hasInitialOrbit      = track->initialOrbit.has_value();
        trackDescription.run                  = describeObject(runCommand);

        // Beams, then the emission sources and distributions of each beam, in beam order.
        for (const std::string& name : track->beamNames_m) {
            const Beam* beam = findAs<Beam>(name);
            if (beam == nullptr) continue;
            run.beams.push_back(describeBeam(*beam));
            const auto* list = findAs<EmissionSourceList>(run.beams.back().emissionSourceList);
            if (list == nullptr) continue;
            for (const EmissionSource* source : list->fetchSources()) {
                if (source == nullptr) continue;
                run.emissionSources.push_back(describeObject(*source));
                if (const auto* dist = findAs<Distribution>(source->getDistributionName())) {
                    run.distributions.push_back(describeObject(*dist));
                }
            }
        }

        if (const Attribute* solver = attributeNamed(runCommand, "FIELDSOLVER")) {
            if (const auto* fs = findAs<FieldSolverCmd>(Attributes::getString(*solver))) {
                run.fieldSolver = describeObject(*fs);
                tryRead([&] {
                    run.track.fieldSolverType = fs->getType();
                });
                if (const BinningCmd* binning = fs->getBinningCmd()) {
                    run.binning = describeObject(*binning);
                }
            }
        }
        if (const Attribute* geometry = attributeNamed(runCommand, "BOUNDARYGEOMETRY")) {
            if (const auto* bg = findAs<BoundaryGeometry>(Attributes::getString(*geometry))) {
                run.boundaryGeometry = describeObject(*bg);
            }
        }

        if (const BeamSequence* use = track->use) run.line = describeLine(*use);
        return run;
    }

    LineDescription describeLine(const BeamSequence& sequence) {
        LineDescription description;
        description.name = sequence.getOpalName();
        if (const auto* line = dynamic_cast<const Line*>(&sequence)) {
            description.topology  = dynamic_cast<const Ring*>(line) ? "RING" : "LINE";
            const auto& origin    = line->getOrigin();
            const auto& direction = line->getDirection();
            description.origin    = {origin[0], origin[1], origin[2]};
            description.direction = {direction[0], direction[1], direction[2], direction[3]};
            std::vector<std::string> path{line->getOpalName()};
            addMembers(*line, false, description.elements, path);
        }
        return description;
    }

}  // namespace opalx::run
