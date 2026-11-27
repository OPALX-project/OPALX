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
         * @param seed OPTION SEED, for the random streams of the decays.
         */
        std::vector<std::unique_ptr<GlobalProcess>> makeGlobalProcessesForBeam(
                const BeamDescription& beam, std::size_t containerIndex, int seed) {
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
                                            tau, containerIndex, mass, parentSign, seed));
                            break;
                        case ParticleType::PION:
                            processes.push_back(
                                    std::make_unique<PionDecay>(
                                            tau, containerIndex, mass, parentSign, seed));
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

        /// What buildTracking() prepares for one beam before the bunch exists.
        struct BeamSetup {
            std::size_t decayRoot = 0;  ///< The beam whose reference particle this one follows.
            std::vector<ObjectValues>
                    sources;  ///< Its emission sources; none for a decay daughter.
            std::vector<std::unique_ptr<GlobalProcess>> processes;  ///< Its global processes.
            BunchBeam bunchBeam;  ///< What PartBunch needs for the beam's container.
        };

        /// Give each container the global processes of its beam.
        void setupGlobalProcesses(PartBunch_t& bunch, std::vector<BeamSetup>& beamSetups) {
            for (std::size_t i = 0; i < beamSetups.size(); ++i) {
                bunch.getParticleContainer(i)->setGlobalProcesses(
                        std::move(beamSetups[i].processes));
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
         * @brief For every beam, the index of the beam whose reference particle it follows
         *        (resolveDecayRoots()): the first beam of its decay chain, or itself.
         * @param beams The beams of the TRACK, in the order of @p beamNames.
         */
        std::vector<std::size_t> decayRootsOf(
                const std::vector<BeamDescription>& beams,
                const std::vector<std::string>& beamNames) {
            std::vector<std::string> daughterNames;
            std::vector<bool> hasDecay;
            for (const BeamDescription& b : beams) {
                daughterNames.push_back(b.daughterBeam);
                const std::vector<std::string>& processes = b.globalProcesses;
                hasDecay.push_back(
                        std::find(processes.begin(), processes.end(), "DECAY") != processes.end());
            }
            return resolveDecayRoots(beamNames, daughterNames, hasDecay);
        }

        /**
         * @brief Wire daughter containers to cross-container processes (e.g. muon decay ->
         * electron).
         * @note Must be called after setupGlobalProcesses - reads processes from the containers.
         */
        void wireDaughterContainers(
                PartBunch_t& bunch, const std::vector<BeamDescription>& beams,
                const std::vector<std::string>& beamNames,
                const std::vector<BeamSetup>& beamSetups) {
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
                const std::size_t root = beamSetups[i].decayRoot;
                bunch.setPcRoot(i, root);
                if (root != i) {
                    *gmsg << level2 << "* Beam \"" << beamNames[i]
                          << "\" uses the reference particle of beam \"" << beamNames[root] << "\"."
                          << endl;
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
            // Common containers / parameters used by all samplers.
            auto pc             = bunch.getParticleContainer(index);
            const auto meshSize = bunch.cartesianDomain().layoutExtents();
            Vector_t<int, Dim> nr;
            for (unsigned d = 0; d < Dim; ++d) {
                nr[d] = static_cast<int>(meshSize[d]);
            }
            const double avrgpz = beam.reference.getP() * Units::eV2GeV / beam.mass;

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

            if (!hasFileSource || !beam.hasExplicitEnergy) {
                return std::nullopt;
            }
            const std::vector<double>& dir = beam.referenceDirection;
            return Vector_t<double, 3>(avrgpz * dir[0], avrgpz * dir[1], avrgpz * dir[2]);
        }

        /// Print the field solver with its boundary conditions, the binning and the beams of the
        /// RUN.
        void printRunSetup(
                const RunDescription& description, const ObjectValues& fieldSolver,
                const std::optional<ObjectValues>& binning) {
            FieldSolverCmd::printInfo(*gmsg << level1, fieldSolver) << endl;
            if (binning) {
                BinningCmd::printInfo(*gmsg << level1, *binning) << endl;
            }
            *gmsg << level2 << FieldSolverCmd::constructBCHandler(fieldSolver) << endl;

            const std::vector<std::string>& beamNames = description.track.beamNames;
            *gmsg << level1 << "* RUN resolved beams: ";
            for (size_t i = 0; i < beamNames.size(); ++i) {
                *gmsg << beamNames[i] << (i + 1 < beamNames.size() ? ", " : "");
            }
            *gmsg << endl;
            for (const BeamDescription& beam : description.beams) {
                *gmsg << level1 << beam.banner << endl;
            }
        }

        /// For each beam: its decay root, emission sources, global processes and what the bunch
        /// needs for its container, with the number of particles to reserve room for.
        std::vector<BeamSetup> setupBeams(const RunDescription& description) {
            const std::vector<BeamDescription>& beams = description.beams;
            const std::vector<std::string>& beamNames = description.track.beamNames;
            const std::vector<std::size_t> decayRoots = decayRootsOf(beams, beamNames);

            std::vector<BeamSetup> beamSetups(beams.size());
            for (size_t i = 0; i < beams.size(); ++i) {
                const BeamDescription& beam = beams[i];
                BeamSetup& setup            = beamSetups[i];

                const double partsPerMacro = beam.chargePerParticle / (beam.charge * Physics::q_e);
                *gmsg << level2 << "* Beam[" << i << "] " << beamNames[i]
                      << " macro charge per particle [C]: " << beam.chargePerParticle << endl;
                *gmsg << level2 << "* Beam[" << i << "] " << beamNames[i]
                      << " macro mass per particle [GeV/c^2]: " << beam.massPerParticle << endl;
                *gmsg << level2 << "* Beam[" << i << "] " << beamNames[i]
                      << " particles per macro particle: " << partsPerMacro << endl
                      << endl;

                setup.decayRoot = decayRoots[i];
                // A decay daughter beam gets particles only from the decay.
                if (setup.decayRoot == i) {
                    for (const SourceDescription& source : beam.sources) {
                        setup.sources.emplace_back(named(description.emissionSources, source.name));
                    }
                }
                setup.processes = makeGlobalProcessesForBeam(beam, i, description.options.seed);

                setup.bunchBeam.chargePerParticle = beam.chargePerParticle;
                setup.bunchBeam.massPerParticle   = beam.massPerParticle;
                setup.bunchBeam.reference         = beam.reference;
                setup.bunchBeam.particleName      = beam.particleName;
                setup.bunchBeam.hasSpin           = beam.hasPolarization;
                setup.bunchBeam.numParticles      = computeTotalAllocationForBunch(
                        beam, setup.sources, description.distributions);
            }
            return beamSetups;
        }

        /// Where the run starts and what each beam emits later.
        struct StartState {
            double sStart           = 0.0;  ///< Path length the run starts at [m].
            double globalPhaseShift = 0.0;  ///< The run starts at minus this time [s].
            std::vector<std::vector<std::shared_ptr<SamplingBase>>> emittingSamplers;
            /// Reference momentum along REFDIRECTION of a beam with a file distribution.
            std::vector<std::optional<Vector_t<double, 3>>> beamReferenceMomenta;
            CheckpointFile::Metadata checkpoint;  ///< Empty for a fresh run.
        };

        /**
         * @brief A restart: particles, reference particles, time and time step come from the
         *        checkpoint, and the output files are rewound to it. Nothing is sampled.
         */
        StartState restoreCheckpoint(
                const RunDescription& description, PartBunch_t& bunch, DataSink& dataSink) {
            const RunContext& context = description.context;
            StartState start;
            start.sStart = description.track.zStart;
            start.emittingSamplers.resize(bunch.getNumParticleContainers());
            start.beamReferenceMomenta.resize(bunch.getNumParticleContainers());
            start.checkpoint = CheckpointFile::read(context.restartFileName, bunch);
            dataSink.rewindToCheckpoint(bunch);
            *gmsg << level1 << "* Restored checkpoint '" << context.restartFileName
                  << "' at global step " << start.checkpoint.globalTrackStep
                  << ", t = " << Util::getTimeString(start.checkpoint.time) << "." << endl;
            return start;
        }

        /**
         * @brief A fresh run: sample the distributions of every beam. The run starts at ZSTART,
         *        or at 0 if the first beam has no particles yet, and the bunch gets the DT of the
         *        segment that holds the start.
         */
        StartState sampleBeams(
                const RunDescription& description, PartBunch_t& bunch,
                const std::vector<BeamSetup>& beamSetups) {
            const TrackDescription& track = description.track;
            StartState start;
            start.emittingSamplers.resize(beamSetups.size());
            start.beamReferenceMomenta.resize(beamSetups.size());
            const std::string inputDirectory =
                    std::filesystem::path(description.context.inputFileName).parent_path().string();
            for (size_t i = 0; i < beamSetups.size(); ++i) {
                start.beamReferenceMomenta[i] = setupDistributionsAndSamplers(
                        bunch, beamSetups[i].sources, description.distributions,
                        description.beams[i], start.emittingSamplers[i], i, start.globalPhaseShift,
                        description.options.seed, inputDirectory);
            }

            start.sStart = track.zStart;
            if (bunch.getParticleContainer()->getTotalNum() == 0) {
                start.sStart = 0.0;
            }
            // The schedule as the tracker builds it; advanceToPos() moves to the first segment
            // whose ZSTOP is at or past the start.
            StepSizeConfig schedule;
            for (std::size_t i = 0; i < track.zStop.size(); ++i) {
                schedule.push_back(track.dT[i], track.zStop[i], track.maxSteps[i]);
            }
            schedule.sortAscendingSStop();
            schedule.resetIterator();
            bunch.setdT(schedule.advanceToPos(start.sStart).getdT());
            return start;
        }

    }  // namespace

    TrackingSetup buildTracking(
            const RunDescription& description, std::vector<std::shared_ptr<ElementBase>> elements,
            RunObjects& objects) {
        // The steps are numbered as in the description of buildTracking() in BuildTracking.h.
        const TrackDescription& track = description.track;
        const RunContext& context     = description.context;
        const OptionValues& options   = description.options;

        // 1. Print the field solver, its boundary conditions, the binning and the beams.
        const ObjectValues fieldSolver(*description.fieldSolver);
        std::optional<ObjectValues> binning;
        if (description.binning) binning.emplace(*description.binning);
        printRunSetup(description, fieldSolver, binning);

        // 2. Per beam: decay root, emission sources, global processes and BunchBeam.
        std::vector<BeamSetup> beamSetups = setupBeams(description);

        // 3. Space-charge configuration from the FIELDSOLVER, the BINNING, all sources and all
        //    distributions; the field-debug dumps are named after the input. The Dirichlet
        //    plane goes to the tracker, the domain configuration to the bunch.
        std::vector<ObjectValues> allSources;
        for (const BeamSetup& beam : beamSetups) {
            allSources.insert(allSources.end(), beam.sources.begin(), beam.sources.end());
        }
        auto spaceChargeConfig = opalx::spacecharge::buildSpaceChargeConfig(
                fieldSolver, binning, allSources,
                std::vector<ObjectValues>(
                        description.distributions.begin(), description.distributions.end()),
                options.repartFreq, options.loadBalancingThreshold,
                options.loadBalancingPairWeight);
        auto& cartesianConfig =
                std::get<opalx::spacecharge::CartesianPIC3DConfig>(spaceChargeConfig);
        cartesianConfig.fieldDumpBasename = context.inputBasename;
        const opalx::spacecharge::DirichletPlaneConfig dirichletPlane =
                cartesianConfig.dirichletPlane;
        const auto cartesianDomainConfig =
                opalx::spacecharge::makeCartesianDomainConfig(spaceChargeConfig);

        // 4. Output: the .stat writers and the .h5 files, one per beam.
        DataSink* dataSink = initDataSink(beamSetups.size(), context, options, objects);

        // 5. The bunch: one particle container per beam, with the global processes of its beam;
        //    decay daughters are connected to their parents.
        std::vector<BunchBeam> bunchBeams;
        for (const BeamSetup& beam : beamSetups) {
            bunchBeams.push_back(beam.bunchBeam);
        }
        objects.bunch = std::make_unique<PartBunch_t>(
                bunchBeams, options.loadBalancingThreshold, cartesianDomainConfig,
                options.useQMAttributes, options.aggressiveStateSync);
        PartBunch_t& bunch = *objects.bunch;
        setupGlobalProcesses(bunch, beamSetups);
        wireDaughterContainers(bunch, description.beams, track.beamNames, beamSetups);

        // 6. The start. A fresh run samples its distributions and sets the start position and
        //    the bunch time step. A restart restores the exact saved particle state instead, so
        //    sampling it again would duplicate particles and consume random numbers.
        StartState start = context.restart ? restoreCheckpoint(description, bunch, *dataSink)
                                           : sampleBeams(description, bunch, beamSetups);

        // 7. The space-charge solver: fields and a warmup solve, also for TYPE = NONE.
        objects.spaceChargeSolver = opalx::spacecharge::makeSpaceChargeSolver(
                std::move(spaceChargeConfig), bunch, dataSink);

        // 8. Start summary of each container.
        std::vector<opalx::moments::BeamMoments> initialMoments;
        for (const auto& pc : bunch.getParticleContainers()) {
            initialMoments.push_back(opalx::moments::compute(*pc, PartBunch_t::printedMoments));
        }
        bunch.print(*gmsg, initialMoments);

        // 9. Everything the tracker reads, in one struct: schedule, start, emission, restart
        //    values, OPTION values, field dumps, the time integrator and the stop rule of each
        //    beam.
        const CheckpointFile::Metadata& checkpoint = start.checkpoint;
        TrackingSetup setup;
        setup.bunch                = &bunch;
        setup.spaceChargeSolver    = objects.spaceChargeSolver.get();
        setup.dirichletPlane       = dirichletPlane;
        setup.dataSink             = dataSink;
        setup.dt                   = track.dT;
        setup.sStop                = track.zStop;
        setup.maxSteps             = track.maxSteps;
        setup.sStart               = start.sStart;
        setup.globalPhaseShift     = start.globalPhaseShift;
        setup.emittingSamplers     = EmittingSamplers(std::move(start.emittingSamplers));
        setup.beamReferenceMomenta = std::move(start.beamReferenceMomenta);
        setup.restarting           = context.restart;
        setup.restartGlobalStep    = static_cast<unsigned long long>(checkpoint.globalTrackStep);
        setup.restartDt            = checkpoint.dt;
        setup.restartPosition      = StepSizeConfig::ResumePosition{
                checkpoint.stepSizeSegment, checkpoint.stepsCompletedInSegment};
        setup.cavityPhases   = checkpoint.cavityPhases;
        setup.checkpointPath = CheckpointFile::defaultPath(context.inputBasename);
        setup.options        = options;
        setup.emFieldDumps   = description.emFieldDumps;
        setup.timeIntegrator = track.timeIntegrator;
        setup.beamStops      = makeBeamStops(description.beams, track, checkpoint.beamStops);

        // 10. The beamline: initialise, place and check the elements and build the element
        //     grid. The tracker steps with a positive time step, and elements read it in
        //     initialise().
        bunch.setdT(std::abs(bunch.getdT()));
        setup.beamline = buildBeamline(
                description.line, elements, bunch, checkpoint.cavityPhases, options, context);
        return setup;
    }

}  // namespace opalx::run
