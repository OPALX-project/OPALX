//
// freeze(): copy the values one RUN command tracks with out of the parsed objects.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include "Track/FreezeRun.h"

#include <exception>
#include <string>

#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/Object.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/PartData.h"
#include "Attributes/Attributes.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "Distribution/Distribution.h"
#include "Lines/EmissionSourceList.h"
#include "Structure/Beam.h"
#include "Structure/BinningCmd.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/EmissionSource.h"
#include "Structure/FieldSolverCmd.h"
#include "Track/Track.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

namespace opalx::run {

    namespace {

        /// Keyword an object was made from: the prototype at the root of its parent chain.
        /// For `Q2: Q1` the parent of Q2 is Q1, whose parent is the QUADRUPOLE prototype.
        std::string keywordOf(const Object& object) {
            const Object* parent = object.getParent();
            if (parent == nullptr) return object.getOpalName();
            while (parent->getParent() != nullptr)
                parent = parent->getParent();
            return parent->getOpalName();
        }

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
         * Walks the tracked line in the order ParallelTracker visits it: nested lines are
         * entered, and the reflection flag of a member (-NAME) does not change the order,
         * because ParallelTracker::visitBeamline iterates every line left to right.
         */
        class LineWalker : public DefaultVisitor {
        public:
            LineWalker(const Beamline& line, std::vector<ElementDescription>& out)
                : DefaultVisitor(line, false, false), out_m(out) {}

            void visitBeamline(const Beamline& line) override { line.iterate(*this, false); }

            void visitFlaggedElmPtr(const FlaggedElmPtr& entry) override {
                // A member of a reflected line is reflected unless it is reflected itself.
                const bool outer = reflected_m;
                reflected_m      = outer != entry.getReflectionFlag();
                DefaultVisitor::visitFlaggedElmPtr(entry);
                reflected_m = outer;
            }

        private:
            void applyDefault(const ElementBase& element) override {
                ElementDescription description;
                description.elementType = element.getTypeString();
                description.occurrence  = out_m.size();
                description.reflected   = reflected_m;
                if (const auto* input = findAs<Element>(element.getName())) {
                    description.object = describeObject(*input);
                } else {
                    description.object.name = element.getName();
                }
                out_m.push_back(std::move(description));
            }

            std::vector<ElementDescription>& out_m;
            bool reflected_m = false;
        };

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
        description.keyword = keywordOf(object);
        description.name    = object.getOpalName();
        description.attributes.reserve(object.itsAttr.size());
        for (const Attribute& attribute : object.itsAttr) {
            if (attribute) description.attributes.push_back(describeAttribute(attribute));
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

        if (const BeamSequence* use = track->use) {
            run.line.name = use->getOpalName();
            if (const Beamline* line = use->fetchLine()) {
                // The same test TrackRun uses to tell a RING from a LINE.
                run.line.topology =
                        line->getBeamlineTopology() == BeamlineTopology::RING ? "RING" : "LINE";
                // Only a RING length is read: TrackRun never evaluates a LINE's length, and
                // evaluating it can fail.
                if (run.line.topology == "RING") {
                    tryRead([&] {
                        run.line.length = use->getLength();
                    });
                }
                LineWalker walker(*line, run.line.elements);
                line->accept(walker);
            }
        }
        return run;
    }

}  // namespace opalx::run
