//
// buildTracking(): make the particle bunch, output, samplers, space-charge solver and step-size
// schedule of one RUN.
//
// Copyright (c) 200x - 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Track/BuildTracking.h"

#include "Moments/BeamMoments.h"

#include "AbsBeamline/RFCavity.h"
#include "Algorithms/DefaultVisitor.h"

#include "Beamlines/TBeamline.h"

#include "Distribution/DistributionSettings.h"

#include "Distribution/EmittedFromFile.h"

#include "Distribution/Gaussian.h"

#include "Distribution/MultiVariateGaussian.h"

#include "Distribution/FlatTop.h"

#include "Distribution/FromFile.h"

#include "Distribution/OpalFlatTop.h"

#include "Distribution/Uniform.h"

#include "Physics/Physics.h"
#include "Physics/Units.h"

#include "Physics/ParticleProperties.h"
#include "Processes/GlobalProcesses/GlobalProcess.h"
#include "Processes/GlobalProcesses/MuonDecay.h"
#include "Processes/GlobalProcesses/PionDecay.h"

#include "Track/DecayChain.h"

#include "Utilities/LogicalError.h"
#include "Utilities/OpalException.h"

#include "Structure/CheckpointFile.h"
#include "Structure/DataSink.h"
#include "Structure/EmissionSource.h"
#include "Structure/FieldSolverCmd.h"
#include "Structure/H5PartWrapper.h"
#include "Structure/H5PartWrapperForPT.h"

#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeConfigBuilder.h"
#include "SpaceCharge/SpaceChargeFactory.h"
#include "SpaceCharge/SpaceChargeSolver.h"

#include "Track/ObjectValues.h"

#include "Utility/Inform.h"

#include "Utilities/BiMap.h"
#include "Utilities/Util.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

extern Inform* gmsg;

namespace opalx::run {

    namespace {
        /**
         * @brief Enforces unit macro weight
         * @note For now, for the moment calculation to give unbiased estimators of the
         * true moments, the macro weight needs to be 1.
         */
        void requireUnitMacroWeight(const BeamDescription& beam, const std::string& role) {
            const double partsPerMacro = beam.chargePerParticle / (beam.charge * Physics::q_e);
            if (std::abs(partsPerMacro - 1.0) > 1e-2) {
                throw OpalException(
                        "TrackRun::execute",
                        "DECAY requires one physical particle per macroparticle, but " + role
                                + " beam \"" + beam.object.name
                                + "\" has particles-per-macro = " + std::to_string(partsPerMacro)
                                + ". Set BCHARGE = NALLOC * |CHARGE| * q_e.");
            }
        }

        /**
         * @brief Builds a vector of processes for the given beam object
         */
        std::vector<std::unique_ptr<GlobalProcess>> makeGlobalProcessesForBeam(
                const BeamDescription& beam, std::size_t containerIndex) {
            std::vector<std::unique_ptr<GlobalProcess>> processes;
            const std::vector<std::string>& processNames = beam.globalProcesses;
            processes.reserve(processNames.size());

            for (const std::string& processName : processNames) {
                if (processName == "DECAY") {
                    const std::string particleName = beam.particleName;
                    const ParticleType pType = ParticleProperties::getParticleType(particleName);
                    const double tau         = ParticleProperties::getParticleLifetime(pType);
                    const double mass        = ParticleProperties::getParticleMass(pType);
                    const double parentQ     = beam.charge;
                    const int parentSign     = (parentQ > 0.0) - (parentQ < 0.0);

                    requireUnitMacroWeight(beam, "parent");

                    switch (pType) {
                        case ParticleType::MUON:
                            if (!beam.hasPolarization) {
                                throw OpalException(
                                        "TrackRun::execute",
                                        "Muon decay requires spin tracking: the differential decay "
                                        "rate is polarization-dependent. Set POLARIZATION = "
                                        "{Px, Py, Pz} on the muon BEAM (this enables spin "
                                        "tracking).");
                            }
                            processes.push_back(
                                    std::make_unique<MuonDecay>(
                                            tau, containerIndex, mass, parentSign));
                            break;
                        case ParticleType::PION:
                            processes.push_back(
                                    std::make_unique<PionDecay>(
                                            tau, containerIndex, mass, parentSign));
                            break;
                        default:
                            throw OpalException(
                                    "TrackRun::execute",
                                    "No decay implementation for PARTICLE=" + particleName
                                            + ". Supported: MUON, PION.");
                    }
                    continue;
                }

                throw OpalException(
                        "TrackRun::execute",
                        "Unknown global process \"" + processName + "\". Supported values: DECAY.");
            }

            return processes;
        }

        DataSink* initDataSink(
                size_t numParticleContainers, const RunContext& context,
                const OptionValues& options, RunObjects& objects) {
            objects.phaseSpaceSinks.clear();
            objects.phaseSpaceSinks.reserve(numParticleContainers);

            const bool isRestart        = context.restart;
            const std::string inputBase = context.inputBasename;

            DataSinkOptions sinkOptions;
            sinkOptions.enableHDF5         = options.enableHDF5;
            sinkOptions.psDumpFreq         = options.psDumpFreq;
            sinkOptions.statDumpFreq       = options.statDumpFreq;
            sinkOptions.stepInfoFreq       = options.stepInfoFreq;
            sinkOptions.printRankDistrFreq = options.printRankDistrFreq;
            sinkOptions.rankDump           = options.rankDump;
            sinkOptions.ebDump             = options.ebDump;
            sinkOptions.rhoDump            = options.rhoDump;
            sinkOptions.beamHaloBoundary   = options.beamHaloBoundary;

            for (size_t i = 0; sinkOptions.enableHDF5 && i < numParticleContainers; ++i) {
                const std::string stem =
                        DataSink::diagnosticStemForContainer(inputBase, numParticleContainers, i);
                const std::string dest = stem + std::string(".h5");

                if (isRestart && std::filesystem::exists(dest)) {
                    objects.phaseSpaceSinks.push_back(
                            std::make_unique<H5PartWrapperForPT>(
                                    dest, sinkOptions, -1, dest, H5_O_RDWR));
                } else {
                    objects.phaseSpaceSinks.push_back(
                            std::make_unique<H5PartWrapperForPT>(dest, sinkOptions, H5_O_WRONLY));
                }
            }

            std::vector<H5PartWrapper*> sinks;
            for (const auto& sink : objects.phaseSpaceSinks)
                sinks.push_back(sink.get());
            // A restart continues the diagnostics under their original basename: existing
            // phase-space H5 files are opened for append, absent ones are started. A later RUN
            // of the same input appends to the statistics files of the earlier ones.
            objects.dataSink = std::make_unique<DataSink>(
                    sinkOptions, sinks, numParticleContainers, inputBase, context.inputFileName,
                    isRestart || context.runIndex > 1);
            return objects.dataSink.get();
        }

        /// Attach the global processes of each beam to its particle container.
        void setupGlobalProcesses(
                PartBunch_t& bunch,
                std::vector<std::vector<std::unique_ptr<GlobalProcess>>> globalProcessesLists) {
            const auto& particleContainers = bunch.getParticleContainers();
            if (particleContainers.size() != globalProcessesLists.size()) {
                throw OpalException(
                        "TrackRun::setupGlobalProcesses",
                        "Mismatch between number of particle containers and global process lists.");
            }

            for (size_t i = 0; i < particleContainers.size(); ++i) {
                if (!particleContainers[i]) {
                    continue;
                }
                particleContainers[i]->setGlobalProcesses(std::move(globalProcessesLists[i]));
            }
        }

        /**
         * @brief Stop condition of every beam: TRACK ZSTOP for the first, BEAM ZSTOP or
         *        UNTILMAXSTEPS for the others.
         * @param restored Planes and done flags from a checkpoint; empty for a fresh run.
         */
        std::vector<BeamStop> makeBeamStops(
                const std::vector<BeamDescription>& beams, const TrackDescription& track,
                const std::vector<BeamStop>& restored) {
            // validate() checked the combinations of TRACK ZSTOP, BEAM ZSTOP and UNTILMAXSTEPS.
            std::vector<BeamStop> stops(beams.size());
            for (std::size_t i = 0; i < beams.size(); ++i) {
                stops[i].zStop         = beams[i].zStop;
                stops[i].untilMaxSteps = beams[i].untilMaxSteps;
            }
            if (track.zStopGiven) stops[0].zStop = track.zStop.back();
            for (std::size_t i = 0; i < stops.size() && i < restored.size(); ++i) {
                // A restart continues with the planes and done flags of the checkpoint.
                stops[i].planeSet    = restored[i].planeSet;
                stops[i].planeOrigin = restored[i].planeOrigin;
                stops[i].planeNormal = restored[i].planeNormal;
                stops[i].done        = restored[i].done;
            }
            return stops;
        }

        /**
         * @brief Wire daughter containers to cross-container processes (e.g. muon decay ->
         * electron).
         * @note Must be called after setupGlobalProcesses - reads processes from the containers.
         */
        void wireDaughterContainers(
                PartBunch_t& bunch, const std::vector<BeamDescription>& beams,
                const std::vector<std::string>& beamNames,
                const std::vector<std::size_t>& decayRoots) {
            const auto& containers = bunch.getParticleContainers();

            for (std::size_t i = 0; i < beams.size(); ++i) {
                const std::string& daughterName = beams[i].daughterBeam;
                if (daughterName.empty()) {
                    continue;
                }

                // Find the container index whose beam name matches DAUGHTERBEAM (resolveDecayRoots
                // has already checked that it exists).
                auto it = std::find(beamNames.begin(), beamNames.end(), daughterName);
                const std::size_t daughterIdx =
                        static_cast<std::size_t>(std::distance(beamNames.begin(), it));

                // Use the physical rest mass from the Beam definition (in GeV), not the
                // macro-particle mass from the container.
                const double daughterMass = beams[daughterIdx].mass;
                for (const auto& proc : containers[i]->getGlobalProcesses()) {
                    auto* decayProc = dynamic_cast<Decay*>(proc.get());
                    if (decayProc) {
                        requireUnitMacroWeight(beams[daughterIdx], "daughter");
                        const double parentCharge   = beams[i].charge;
                        const double daughterCharge = beams[daughterIdx].charge;
                        if (daughterCharge != parentCharge) {
                            std::ostringstream msg;
                            msg << "Daughter beam \"" << daughterName
                                << "\" has CHARGE = " << daughterCharge
                                << ", but its parent beam \"" << beamNames[i]
                                << "\" has CHARGE = " << parentCharge
                                << ". A decay keeps the charge: set the same CHARGE on both beams "
                                   "(for "
                                   "mu+ -> e+, CHARGE = +1 on the ELECTRON beam).";
                            throw OpalException("TrackRun::wireDaughterContainers", msg.str());
                        }
                        // A muon daughter (e.g. from pion decay) receives a per-particle
                        // polarization from the decay, so its container must have spin storage -
                        // which is enabled by setting POLARIZATION on the daughter muon BEAM.
                        const ParticleType daughterType = ParticleProperties::getParticleType(
                                beams[daughterIdx].particleName);
                        if (daughterType == ParticleType::MUON
                            && !beams[daughterIdx].hasPolarization) {
                            throw OpalException(
                            "TrackRun::wireDaughterContainers",
                            "Decay produces muons in daughter beam \"" + daughterName
                                    + "\", whose polarization must be tracked. Set POLARIZATION "
                                      "= {Px, Py, Pz} on that BEAM to enable spin tracking.");
                        }
                        decayProc->setDaughterContainer(containers[daughterIdx], daughterMass);
                        *gmsg << level2 << "* Wired decay on beam \"" << beamNames[i]
                              << "\" to daughter beam \"" << daughterName << "\" (container "
                              << daughterIdx << ")." << endl;
                    }
                }
            }

            // Decay daughter containers use the reference particle, frame, path length and orbit of
            // the first container of their decay chain.
            for (std::size_t i = 0; i < beams.size(); ++i) {
                bunch.setPcRoot(i, decayRoots[i]);
                if (decayRoots[i] != i) {
                    *gmsg << level2 << "* Beam \"" << beamNames[i]
                          << "\" uses the reference particle of beam \"" << beamNames[decayRoots[i]]
                          << "\"." << endl;
                }
            }
        }

        /// Set the phases of the RF cavities that a checkpoint stored; autophase leaves these
        /// cavities alone.
        void restoreCavityPhases(OpalBeamline& beamline, const CavityPhases& phases) {
            ElementList cavities       = beamline.getElementByType(ElementType::RFCAVITY);
            ElementList travelingWaves = beamline.getElementByType(ElementType::TRAVELINGWAVE);
            cavities.insert(cavities.end(), travelingWaves.begin(), travelingWaves.end());
            for (auto it = phases.begin(); it != phases.end(); ++it) {
                for (const auto& cavity : cavities) {
                    if (cavity->getName() != it->first) continue;
                    auto* element = static_cast<RFCavity*>(cavity.get());
                    element->setPhasem(it->second);
                    element->setAutophaseVeto();
                    *gmsg << level3
                          << "Restored cavity phase from the h5 file. Name: " << element->getName()
                          << ", phase: " << it->second << " rad" << endl;
                    break;
                }
            }
        }

        /**
         * @brief Hand the elements of the line to an OpalBeamline in line order (addElement()
         *        calls initialise()), place them, write the 3D lattice files and build the
         *        element mapping.
         *
         * Every element gets the OPTION values and the run context first (setRunValues()); an
         * element without REACH gets OPTION ELEMENTMAPPINGREACH before it is added.
         * @param cavityPhases The cavity phases a checkpoint stored; empty for a fresh run.
         */
        OpalBeamline buildBeamline(
                const LineDescription& line,
                const std::vector<std::shared_ptr<ElementBase>>& elements, PartBunch_t& bunch,
                const CavityPhases& cavityPhases, const OptionValues& options,
                const RunContext& context) {
            OpalBeamline beamline(
                    Vector_t<double, 3>(line.origin[0], line.origin[1], line.origin[2]),
                    Quaternion(
                            line.direction[0], line.direction[1], line.direction[2],
                            line.direction[3]));
            for (const auto& built : elements) {
                built->setRunValues(options, context);
                switch (built->trackedAs()) {
                    case TrackedAs::FieldElement:
                        if (!built->getReach()) built->setReach(options.elementMappingReach);
                        beamline.addElement(built, bunch);
                        break;
                    case TrackedAs::Unsupported:
                        throw LogicalError(
                                "opalx::run::buildBeamline",
                                "Tracking of the \"" + Util::toUpper(built->getTypeString())
                                        + "\" element is not implemented yet.");
                    case TrackedAs::Ignored:
                        break;
                }
            }

            // Sorts the elements by starting position and places them.
            beamline.prepareSections();
            beamline.compute3DLattice();

            beamline.save3DLattice(context.inputBasename, context.openMode == "APPEND");
            beamline.save3DInput(context.inputFileName, context.inputBasename);

            restoreCavityPhases(beamline, cavityPhases);

            // Placement, misalignment and field-map extents are final here: each element's
            // aperture must lie inside its reach before the grid is built from the reaches.
            for (const auto& element : beamline.getElementByType(ElementType::ANY))
                element->checkRegions();
            beamline.buildElementMapping(options.elementMappingGrid, options.elementMappingReach);
            return beamline;
        }

        /// The description named @p name in @p objects; freeze() recorded it.
        const ObjectDescription& named(
                const std::vector<ObjectDescription>& objects, const std::string& name) {
            for (const ObjectDescription& object : objects) {
                if (object.name == name) return object;
            }
            throw LogicalError("opalx::run::buildTracking", "No description of \"" + name + "\".");
        }

        /// Compute total number of macroparticles for the bunch from BEAM::NALLOC and
        /// optional per-distribution NPARTDIST values on the emission sources.
        size_t computeTotalAllocationForBunch(
                const BeamDescription& beam, const std::vector<ObjectValues>& sources,
                const std::vector<ObjectDescription>& distributions) {
            size_t beamAllocSize = beam.numAlloc;

            size_t totalFromDists = 0;
            for (const ObjectValues& src : sources) {
                totalFromDists += DistributionSettings::numParticles(ObjectValues(
                        named(distributions, EmissionSource::getDistributionName(src))));
            }

            if (totalFromDists > 0) {
                *gmsg << level3
                      << "* Sum of per-distribution NPARTDIST over all emission sources = "
                      << totalFromDists << ", BEAM::NALLOC = " << beamAllocSize << endl;
                if (totalFromDists > beamAllocSize) {
                    *gmsg << level1 << "* WARNING: Sum of NPARTDIST over all distributions ("
                          << totalFromDists << ") exceeds BEAM::NALLOC (" << beamAllocSize
                          << "). Allocation baseline may be insufficient; "
                          << "macro-charge per particle is still derived from BEAM::NALLOC."
                          << endl;
                }
                return totalFromDists;
            }

            return beamAllocSize;
        }

        /// Build samplers for all emission sources, perform initial sampling for t0 == 0
        /// sources, and populate emittingSamplers for time-dependent or delayed sources.
        /// Applied to particle container [index]
        /// Returns the reference momentum [beta*gamma] set on BEAM for a file distribution
        /// (PC, ENERGY or GAMMA along REFDIRECTION), or nothing to use the bunch mean.
        /// The samplers get OPTION SEED @p seed and the input file directory @p inputDirectory.
        std::optional<Vector_t<double, 3>> setupDistributionsAndSamplers(
                PartBunch_t& bunch, const std::vector<ObjectValues>& sources,
                const std::vector<ObjectDescription>& distributions, const BeamDescription& beam,
                std::vector<std::shared_ptr<SamplingBase>>& emittingSamplers, size_t index,
                double& globalPhaseShift, int seed, const std::string& inputDirectory) {
            static IpplTimings::TimerRef samplingTime = IpplTimings::getTimer("samplingTime");

            IpplTimings::startTimer(samplingTime);

            // Common containers / parameters used by all samplers.
            auto pc               = bunch.getParticleContainer(index);
            Vector_t<int, Dim> nr = bunch.nr_m;
            const double avrgpz   = beam.reference.getP() * Units::eV2GeV / beam.mass;

            emittingSamplers.clear();
            bool hasFileSource = false;

            for (const ObjectValues& src : sources) {
                const DistributionSettings opalDist(
                        ObjectValues(
                                named(distributions, EmissionSource::getDistributionName(src))),
                        avrgpz, seed, inputDirectory);

                // File-based distributions carry absolute momenta. PC/ENERGY/GAMMA on the BEAM is
                // optional for them and only sets the reference momentum (see the return value).
                hasFileSource = hasFileSource || opalDist.getType() == DistributionType::FROMFILE
                                || opalDist.getType() == DistributionType::EMITTEDFROMFILE;

                // Build a sampler instance for this emission source.
                std::shared_ptr<SamplingBase> sampler;
                switch (opalDist.getType()) {
                    case DistributionType::UNIFORM:
                        sampler = std::make_shared<Uniform>(pc, opalDist);
                        break;
                    case DistributionType::GAUSS:
                        sampler = std::make_shared<Gaussian>(pc, opalDist);
                        break;
                    case DistributionType::MULTIVARIATEGAUSS:
                        sampler = std::make_shared<MultiVariateGaussian>(pc, opalDist);
                        break;
                    case DistributionType::FLATTOP:
                        sampler = std::make_shared<FlatTop>(pc, opalDist);
                        break;
                    case DistributionType::OPALFLATTOP:
                        sampler = std::make_shared<OpalFlatTop>(pc, opalDist);
                        break;
                    case DistributionType::FROMFILE:
                        sampler = std::make_shared<FromFile>(pc, opalDist);
                        break;
                    case DistributionType::EMITTEDFROMFILE:
                        sampler = std::make_shared<EmittedFromFile>(pc, opalDist);
                        break;
                    default:
                        throw OpalException(
                                "Distribution::create", "Unknown \"TYPE\" of \"DISTRIBUTION\"");
                }

                // Per-source emission state. The distribution remains a shape sampler; source
                // physics is folded into the sampler setup here before particles are generated.
                const auto R0               = EmissionSource::getR0(src);
                auto P0                     = EmissionSource::getP0(src);
                const double t0             = EmissionSource::getT0(src);
                const std::string emitModel = EmissionSource::getEmissionModel(src);
                const bool generatedEmission =
                        opalDist.emitting_m && opalDist.getType() != DistributionType::FROMFILE
                        && opalDist.getType() != DistributionType::EMITTEDFROMFILE;
                const double eKin                = EmissionSource::getKineticEnergy(src);
                double emissionMomentumMagnitude = 0.0;
                if (generatedEmission && eKin > 0.0) {
                    // EKIN is source-side kinetic energy in eV. Convert it to normalized momentum.
                    // For EMISSIONMODEL=NONE it is added to P0Z, reproducing old OPAL's
                    // cathode-energy increment. For ASTRA, the normalized half-sphere samples are
                    // stretched to this magnitude in the sampler, then the P0 offset is added to
                    // every sampled momentum.
                    const double thermalP = Util::getBetaGamma(eKin, beam.mass * Units::GeV2eV);
                    if (emitModel == "NONE") {
                        P0[2] += thermalP;
                    } else if (emitModel == "ASTRA") {
                        emissionMomentumMagnitude = thermalP;
                    }
                }
                sampler->setEmissionOffsets(R0, P0, t0, emitModel, emissionMomentumMagnitude);

                // Initial polarization from BEAM (ignored if container has no spin attribute).
                const std::vector<double>& pol = beam.polarization;
                sampler->setInitialPolarization({pol[0], pol[1], pol[2]});

                const size_t Ndist = opalDist.getNumParticles();
                size_t Nmutable    = Ndist;

                // Always call generateParticles once per source; time-independent samplers
                // will internally early-return when t0 > 0, while FlatTop will set up its
                // emission structures irrespective of t0.
                sampler->generateParticles(Nmutable, nr);

                globalPhaseShift = std::max(globalPhaseShift, sampler->getGlobalTimeShift());

                // Time-dependent (emitted) distributions (e.g. FlatTop) and delayed
                // one-shot injectors (t0 > 0) participate in emitParticles(t, dt)
                // during tracking.
                if (opalDist.emitting_m || t0 > 0.0
                    || opalDist.getType() == DistributionType::EMITTEDFROMFILE) {
                    emittingSamplers.push_back(sampler);
                    *gmsg << level2 << "* Configured emitting source of type "
                          << opalDist.getTypeofDistribution() << " with NPARTDIST = " << Ndist
                          << ", t0 = " << t0 << endl;
                }
            }

            *gmsg << level2 << "* Particle sampling / sampler setup for all emission sources done."
                  << endl;
            IpplTimings::stopTimer(samplingTime);

            if (!hasFileSource || !beam.hasExplicitEnergy) {
                return std::nullopt;
            }
            const std::vector<double>& dir = beam.referenceDirection;
            return Vector_t<double, 3>(avrgpz * dir[0], avrgpz * dir[1], avrgpz * dir[2]);
        }

    }  // namespace

    TrackingSetup buildTracking(
            const RunDescription& description, std::vector<std::shared_ptr<ElementBase>> elements,
            RunObjects& objects) {
        const ObjectValues run(description.track.run);

        // validate() checked that the RUN names a FIELDSOLVER.
        const ObjectValues fieldSolver(*description.fieldSolver);
        std::optional<ObjectValues> binning;
        if (description.binning) binning.emplace(*description.binning);
        FieldSolverCmd::printInfo(*gmsg << level1, fieldSolver) << endl;
        if (binning) {
            BinningCmd::printInfo(*gmsg << level1, *binning) << endl;
        }

        const TrackDescription& track             = description.track;
        const RunContext& context                 = description.context;
        const OptionValues& options               = description.options;
        const std::vector<std::string>& beamNames = track.beamNames;

        // Create vector of BEAMs (validate() checked that they exist).
        const std::vector<BeamDescription>& beams = description.beams;
        // Decay chains: a beam named as another beam's DAUGHTERBEAM uses the reference particle of
        // the first beam of its chain. It has no momentum and no emission sources of its own.
        std::vector<std::size_t> decayRoots;
        {
            std::vector<std::string> daughterNames;
            std::vector<bool> hasDecay;
            for (const BeamDescription& b : beams) {
                daughterNames.push_back(b.daughterBeam);
                const std::vector<std::string>& processes = b.globalProcesses;
                hasDecay.push_back(
                        std::find(processes.begin(), processes.end(), "DECAY") != processes.end());
            }
            decayRoots = resolveDecayRoots(beamNames, daughterNames, hasDecay);
        }
        *gmsg << level1 << "* RUN resolved beams: ";
        for (size_t i = 0; i < beamNames.size(); ++i) {
            *gmsg << beamNames[i] << (i + 1 < beamNames.size() ? ", " : "");
        }
        *gmsg << endl;
        // Print the BEAM banner for each resolved beam.
        for (const BeamDescription& b : beams) {
            *gmsg << level1 << b.banner << endl;
        }

        // Vectors for each species
        std::vector<double> macrocharges;
        std::vector<double> macromasses;
        std::vector<std::vector<ObjectValues>> emissionSourcesLists;
        std::vector<std::vector<std::unique_ptr<GlobalProcess>>> globalProcessesLists;
        macrocharges.reserve(beams.size());
        macromasses.reserve(beams.size());
        emissionSourcesLists.reserve(beams.size());
        globalProcessesLists.resize(beams.size());

        // Fill macro quantities and emissionSourceList per container (beam)
        for (size_t i = 0; i < beams.size(); ++i) {
            const BeamDescription& b = beams[i];

            const double macrocharge = b.chargePerParticle;
            const double macromass   = b.massPerParticle;
            macrocharges.push_back(macrocharge);
            macromasses.push_back(macromass);

            const double part_per_macro_ratio = macrocharge / (b.charge * Physics::q_e);
            *gmsg << level2 << "* Beam[" << i << "] " << beamNames[i]
                  << " macro charge per particle [C]: " << macrocharge << endl;
            *gmsg << level2 << "* Beam[" << i << "] " << beamNames[i]
                  << " macro mass per particle [GeV/c^2]: " << macromass << endl;
            *gmsg << level2 << "* Beam[" << i << "] " << beamNames[i]
                  << " particles per macro particle: " << part_per_macro_ratio << endl
                  << endl;

            if (decayRoots[i] != i) {
                // A decay daughter beam gets particles only from the decay.
                emissionSourcesLists.emplace_back();
            } else {
                emissionSourcesLists.emplace_back();
                for (const SourceDescription& source : b.sources) {
                    emissionSourcesLists.back().emplace_back(
                            named(description.emissionSources, source.name));
                }
            }

            globalProcessesLists[i] = makeGlobalProcessesForBeam(b, i);
        }

        // validate() refused emitting or delayed sources and DECAY for a restart.
        const bool isRestart = context.restart;

        std::vector<ObjectValues> allSources;
        for (const auto& sources : emissionSourcesLists) {
            allSources.insert(allSources.end(), sources.begin(), sources.end());
        }
        auto spaceChargeConfig = opalx::spacecharge::buildSpaceChargeConfig(
                fieldSolver, binning, allSources,
                std::vector<ObjectValues>(
                        description.distributions.begin(), description.distributions.end()),
                options.repartFreq, options.loadBalancingThreshold,
                options.loadBalancingPairWeight);
        const opalx::spacecharge::DirichletPlaneConfig dirichletPlane =
                std::get<opalx::spacecharge::CartesianPIC3DConfig>(spaceChargeConfig)
                        .dirichletPlane;
        const auto cartesianDomainConfig =
                opalx::spacecharge::makeCartesianDomainConfig(spaceChargeConfig);

        /*
        Need the following units for mass and charge:
        - Charge per macro particle in [C], this should be macrocharge_m or q_m in the bunch.
          This will be used for the field calculations.
        - The pusher needs consistent units: eV for mass and elementary charges for charge.
          This will (hopefully) be handled inside the pusher routines!
        */

        DataSink* dataSink = initDataSink(beams.size(), context, options, objects);

        // Set total particles per container (beam)
        std::vector<size_t> totalParticlesPerBeam(beams.size());
        for (size_t i = 0; i < beams.size(); ++i) {
            totalParticlesPerBeam[i] = computeTotalAllocationForBunch(
                    beams[i], emissionSourcesLists[i], description.distributions);
        }

        // Create PartBunch with multiple particle containers.
        objects.bunch = std::make_unique<PartBunch_t>(
                macrocharges,                    // Macro charge [C]
                macromasses,                     // Macro Mass [GeV]
                beams,                           // Beam descriptions per container
                totalParticlesPerBeam,           // Per-beam particle counts for allocation
                options.loadBalancingThreshold,  // Load balancing threshold
                "LF2",                           // Integrator
                cartesianDomainConfig);          // Cartesian domain and particle-layout setup

        // Validate container setup produced by constructor
        PartBunch_t& bunch             = *objects.bunch;
        const auto& particleContainers = bunch.getParticleContainers();
        if (particleContainers.size() != beams.size()) {
            throw OpalException(
                    "TrackRun::execute",
                    "Mismatch between number of beams and particle containers.");
        }

        // Global processes
        setupGlobalProcesses(bunch, std::move(globalProcessesLists));
        wireDaughterContainers(bunch, beams, beamNames, decayRoots);

        // BC handler
        *gmsg << level2 << FieldSolverCmd::constructBCHandler(fieldSolver) << endl;

        // A fresh run samples its distributions. A restart restores the exact saved particle state
        // instead, so sampling it again would duplicate particles and consume random numbers.
        double globalPhaseShift = 0.0;
        std::vector<std::vector<std::shared_ptr<SamplingBase>>> emittingSamplersList(
                particleContainers.size());
        std::vector<std::optional<Vector_t<double, 3>>> beamReferenceMomenta(
                particleContainers.size());
        CheckpointFile::Metadata restartMetadata;
        if (isRestart) {
            restartMetadata = CheckpointFile::read(context.restartFileName, bunch);
            dataSink->rewindToCheckpoint(bunch);
            *gmsg << level1 << "* Restored checkpoint '" << context.restartFileName
                  << "' at global step " << restartMetadata.globalTrackStep
                  << ", t = " << Util::getTimeString(restartMetadata.time) << "." << endl;
        } else {
            const std::string inputDirectory =
                    std::filesystem::path(context.inputFileName).parent_path().string();
            for (size_t i = 0; i < particleContainers.size(); ++i) {
                beamReferenceMomenta[i] = setupDistributionsAndSamplers(
                        bunch, emissionSourcesLists[i], description.distributions, beams[i],
                        emittingSamplersList[i], i, globalPhaseShift, options.seed, inputDirectory);
            }
        }
        objects.spaceChargeSolver = opalx::spacecharge::makeSpaceChargeSolver(
                std::move(spaceChargeConfig), bunch, dataSink);

        if (!isRestart) {
            // Refresh the initial particle statistics after distribution setup.
            bunch.setCharge();
            bunch.setMass();
        }

        std::vector<opalx::moments::BeamMoments> initialMoments;
        for (const auto& pc : bunch.getParticleContainers()) {
            initialMoments.push_back(opalx::moments::compute(*pc, PartBunch_t::printedMoments));
        }
        bunch.print(*gmsg, initialMoments);

        // Set ZStart, ZStop, and dT
        double sStart = track.zStart;
        if (!isRestart && bunch.getParticleContainer()->getTotalNum() > 0) {
            double spos       = track.zStart;
            const auto& zstop = track.zStop;
            auto it           = track.dT.begin();

            unsigned int i = 0;
            while (i + 1 < zstop.size() && zstop[i + 1] < spos) {
                ++i;
                ++it;
            }

            bunch.setdT(*it);
        } else if (!isRestart) {
            sStart = 0.0;
        }

        const auto& maxSteps = track.maxSteps;
        const auto& sStop    = track.zStop;
        // Everything the tracker reads, in one struct.
        TrackingSetup setup;
        setup.bunch                = &bunch;
        setup.spaceChargeSolver    = objects.spaceChargeSolver.get();
        setup.dirichletPlane       = dirichletPlane;
        setup.dataSink             = dataSink;
        setup.dt                   = track.dT;
        setup.sStop                = sStop;
        setup.maxSteps             = maxSteps;
        setup.sStart               = sStart;
        setup.globalPhaseShift     = globalPhaseShift;
        setup.emittingSamplers     = EmittingSamplers(std::move(emittingSamplersList));
        setup.beamReferenceMomenta = std::move(beamReferenceMomenta);
        setup.restarting           = isRestart;
        setup.restartGlobalStep = static_cast<unsigned long long>(restartMetadata.globalTrackStep);
        setup.restartDt         = restartMetadata.dt;
        setup.restartPosition   = StepSizeConfig::ResumePosition{
                restartMetadata.stepSizeSegment, restartMetadata.stepsCompletedInSegment};
        setup.cavityPhases   = restartMetadata.cavityPhases;
        setup.checkpointPath = CheckpointFile::defaultPath(context.inputBasename);
        setup.options        = description.options;
        setup.emFieldDumps   = description.emFieldDumps;

        // The tracker steps with a positive time step; elements read it in initialise().
        bunch.setdT(std::abs(bunch.getdT()));
        setup.beamline = buildBeamline(
                description.line, elements, bunch, restartMetadata.cavityPhases, options, context);
        setup.spaceChargeFieldUpdate = run.getString("SCFIELDUPDATE") == "PRESTEP"
                                               ? SpaceChargeFieldUpdate::PRESTEP
                                               : SpaceChargeFieldUpdate::MIDPOINT;
        setup.beamStops              = makeBeamStops(beams, track, restartMetadata.beamStops);
        return setup;
    }

}  // namespace opalx::run
