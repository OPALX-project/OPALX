//
// Class Option
//   The OPTION command.
//   The user interface allowing setting of OPAL options.
//   The actual option flags are contained in namespace Options.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "BasicActions/Option.h"
#include "Ippl.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/GeneralElementMapping.h"
#include "Attributes/Attributes.h"
#include "OpalParser/FileStream.h"
#include "Utilities/ExpressionRandom.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include "Utility/Inform.h"
#include "Utility/IpplInfo.h"

#include "Utilities/BiMap.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <limits>
#include <set>

extern Inform* gmsg;

using namespace Options;

namespace {
    // The attributes of class Option.
    enum {
        ECHO,
        INFO,
        TRACE,
        WARN,
        SEED,
        TELL,
        PSDUMPFREQ,
        STATDUMPFREQ,
        CHECKPOINTFREQ,
        STEPINFOFQ,
        PRINTRANKDISTRFQ,
        PSDUMPEACHTURN,
        PSDUMPFRAME,
        SPTDUMPFREQ,
        REPARTFREQ,
        REBINFREQ,
        SCSOLVEFREQ,
        MTSSUBSTEPS,
        REMOTEPARTDEL,
        RHODUMP,
        EBDUMP,
        RANKDUMP,
        CSRDUMP,
        AUTOPHASE,
        NUMBLOCKS,
        RECYCLEBLOCKS,
        NLHS,
        CZERO,
        RNGTYPE,
        ENABLEHDF5,
        ENABLEVTK,
        ENABLELINEARTRANSFERMAPS,
        ELEMENTMAPPINGGRID,
        ELEMENTMAPPINGREACH,
        ASCIIDUMP,
        BOUNDPDESTROY,
        BEAMHALOBOUNDARY,
        CLOTUNEONLY,
        IDEALIZED,
        LOGBENDTRAJECTORY,
        VERSION,
        MEMORYDUMP,
        HALOSHIFT,
        DELPARTFREQ,
        MINBINEMITTED,
        MINSTEPFORREBIN,
        COMPUTEPERCENTILES,
        QM_MODE,
        AGGRESSIVE_STATE_SYNC,
        LOADBALANCINGTHRESHOLD,
        LOADBALANCINGPAIRWEIGHT,
        SIZE
    };
}  // namespace

Option::Option()
    : Action(SIZE, "OPTION", "The \"OPTION\" statement defines OPAL execution options.") {
    itsAttr[ECHO] = Attributes::makeBool("ECHO", "If true, give echo of input", echo);

    itsAttr[INFO] = Attributes::makeBool("INFO", "If true, print information messages", info);

    itsAttr[TRACE] = Attributes::makeBool(
            "TRACE",
            "If true, print execution trace"
            "Must be the first option in the inputfile in "
            "order to render correct results",
            mtrace);

    itsAttr[WARN] = Attributes::makeBool("WARN", "If true, print warning messages", warn);

    itsAttr[SEED] = Attributes::makeReal(
            "SEED", "The seed for the random generator, -1 will use time(0) as seed ");

    itsAttr[TELL] = Attributes::makeBool(
            "TELL",
            "If true, print the current settings. "
            "Must be the last option in the inputfile in "
            "order to render correct results",
            false);

    itsAttr[PSDUMPFREQ] = Attributes::makeReal(
            "PSDUMPFREQ",
            "The frequency to dump the phase space, "
            "i.e. dump data when step%psDumpFreq==0; 0 disables phase space dumps. "
            "Its default value is 10.",
            psDumpFreq);

    itsAttr[STATDUMPFREQ] = Attributes::makeReal(
            "STATDUMPFREQ",
            "The frequency to dump statistical data "
            "(e.g. RMS beam quantities), i.e. dump data when step%statDumpFreq == 0, "
            "0 disables statistical dumps. Its default value is 10.",
            statDumpFreq);

    itsAttr[CHECKPOINTFREQ] = Attributes::makeReal(
            "CHECKPOINTFREQ",
            "The frequency to write a restart checkpoint after completed integration steps, "
            "i.e. checkpoint when GlobalTrackStep % CHECKPOINTFREQ == 0. "
            "A value of 0 disables checkpointing. Its default value is 0.",
            checkpointFreq);

    itsAttr[STEPINFOFQ] = Attributes::makeReal(
            "STEPINFOFQ",
            "The frequency to print per-step tracking status lines. "
            "A value of 0 disables these status lines; its default value is 1.",
            stepInfoFreq);

    itsAttr[PRINTRANKDISTRFQ] = Attributes::makeReal(
            "PRINTRANKDISTRFQ",
            "The frequency to print per-rank particle distribution tables. "
            "A value of 0 disables these tables; its default value is 0.",
            printRankDistrFreq);

    itsAttr[PSDUMPEACHTURN] = Attributes::makeBool(
            "PSDUMPEACHTURN",
            "If true, dump phase space after each "
            "turn ,only aviable for OPAL-cycl, its default value is false",
            false);

    itsAttr[SCSOLVEFREQ] = Attributes::makeReal(
            "SCSOLVEFREQ", "The frequency to solve space charge fields. its default value is 1");

    itsAttr[MTSSUBSTEPS] = Attributes::makeReal(
            "MTSSUBSTEPS",
            "How many small timesteps "
            "are inside the large timestep used in multiple "
            "time stepping (MTS) integrator");

    itsAttr[REMOTEPARTDEL] = Attributes::makeReal(
            "REMOTEPARTDEL",
            "Artifically delete the remote particle "
            "if its distance to the beam mass is larger than "
            "REMOTEPARTDEL times of the beam rms size, "
            "its default values is 0 (no delete) ",
            0.0);

    itsAttr[PSDUMPFRAME] = Attributes::makePredefinedString(
            "PSDUMPFRAME",
            "Controls the frame of phase space dump in "
            "stat file and h5 file. If 'GLOBAL' OPAL will dump in the "
            "lab (global) Cartesian frame; if 'BUNCH_MEAN' OPAL will "
            "dump in the local Cartesian frame of the beam mean; "
            "if 'REFERENCE'  OPAL will dump in the local Cartesian "
            "frame of the reference particle 0. Only available for "
            "OPAL-cycl.",
            {"BUNCH_MEAN", "REFERENCE", "GLOBAL"}, "GLOBAL");

    itsAttr[SPTDUMPFREQ] = Attributes::makeReal(
            "SPTDUMPFREQ",
            "The frequency to dump single "
            "particle trajectory of particles with ID = 0 & 1, "
            "its default value is 1.",
            1);

    itsAttr[REPARTFREQ] = Attributes::makeReal(
            "REPARTFREQ",
            "The frequency to do particles repartition "
            "for better load balance between nodes, its "
            "default value is "
                    + std::to_string(repartFreq) + ".",
            repartFreq);

    itsAttr[MINBINEMITTED] = Attributes::makeReal(
            "MINBINEMITTED",
            "The number of bins that have to be emitted before the bins are squashed into "
            "a single bin; the default value is "
                    + std::to_string(minBinEmitted) + ".",
            minBinEmitted);

    itsAttr[MINSTEPFORREBIN] = Attributes::makeReal(
            "MINSTEPFORREBIN",
            "The number of steps into the simulation before the bins are squashed into "
            "a single bin; the default value is 200.",
            200);

    itsAttr[REBINFREQ] = Attributes::makeReal(
            "REBINFREQ",
            "The frequency to reset energy bin ID for "
            "all particles, its default value is 100.",
            100);

    itsAttr[RHODUMP] = Attributes::makeBool(
            "RHODUMP",
            "If true, in addition to the phase "
            "space the scalar rho field is also dumped (H5Block)",
            rhoDump);

    itsAttr[EBDUMP] = Attributes::makeBool(
            "EBDUMP",
            "If true, in addition to the phase space the "
            "E and B field at each particle is also dumped into the H5 file)",
            ebDump);

    itsAttr[RANKDUMP] = Attributes::makeBool(
            "RANKDUMP",
            "If true, the current MPI rank for each particle is dumped into the H5 file. "
            "Default: false.",
            rankDump);

    itsAttr[CSRDUMP] = Attributes::makeBool(
            "CSRDUMP",
            "If true, the csr E field, line density "
            "and the line density derivative is dumped into the "
            "data directory)",
            false);

    itsAttr[AUTOPHASE] = Attributes::makeReal(
            "AUTOPHASE",
            "If greater than zero OPAL is scanning "
            "the phases of each rf structure in order to get maximum "
            "acceleration. Defines the number of refinements of the "
            "search range",
            autoPhase);

    itsAttr[CZERO] = Attributes::makeBool(
            "CZERO",
            "If set to true a symmetric distribution is "
            "created -> centroid == 0.0",
            false);

    itsAttr[RNGTYPE] = Attributes::makePredefinedString(
            "RNGTYPE",
            "Type of pseudo- or quasi-random number generator, "
            "see also Quasi-Random Sequences, GSL reference manual.",
            {"RANDOM", "HALTON", "SOBOL", "NIEDERREITER"}, "RANDOM");

    itsAttr[CLOTUNEONLY] = Attributes::makeBool(
            "CLOTUNEONLY",
            "If set to true stop after "
            "CLO and tune calculation ",
            false);

    itsAttr[NUMBLOCKS] = Attributes::makeReal(
            "NUMBLOCKS",
            "Maximum number of vectors in the Krylov "
            "space (for RCGSolMgr). Default value is 0 and BlockCGSolMgr will be used.");

    itsAttr[RECYCLEBLOCKS] = Attributes::makeReal(
            "RECYCLEBLOCKS",
            "Number of vectors in the recycle "
            "space (for RCGSolMgr). Default value is 0 and BlockCGSolMgr will be used.");

    itsAttr[NLHS] = Attributes::makeReal(
            "NLHS",
            "Number of stored old solutions for extrapolating "
            "the new starting vector. Default value is 1 and just the last solution is used.");

    itsAttr[ENABLEHDF5] =
            Attributes::makeBool("ENABLEHDF5", "If true, HDF5 actions are enabled", enableHDF5);

    itsAttr[ENABLEVTK] = Attributes::makeBool(
            "ENABLEVTK", "If true, writing of VTK files are enabled", enableVTK);

    itsAttr[ENABLELINEARTRANSFERMAPS] = Attributes::makeBool(
            "ENABLELINEARTRANSFERMAPS",
            "Linear transfer maps were removed; only FALSE is accepted.", false);

    itsAttr[ELEMENTMAPPINGGRID] = Attributes::makeRealArray(
            "ELEMENTMAPPINGGRID",
            "Cells of the lab-frame grid that maps a bunch to its active elements, along x, y "
            "and z: three positive integers, at most 2^24 cells in total. Default {8, 8, 8}.");
    Attributes::setRealArray(
            itsAttr[ELEMENTMAPPINGGRID],
            std::vector<double>(elementMappingGrid.begin(), elementMappingGrid.end()));
    itsAttr[ELEMENTMAPPINGREACH] = Attributes::makeReal(
            "ELEMENTMAPPINGREACH",
            "Minimum transverse half-width [m] of an element's box for the element mapping; "
            "an element's REACH overrides it. Parallel branches must be separated by more than "
            "the reach plus the bunch's transverse half-width, halo included. Default 0.5.",
            elementMappingReach);

    itsAttr[ASCIIDUMP] = Attributes::makeBool(
            "ASCIIDUMP", "If true, some of the elements dump in ASCII instead of HDF5", asciidump);

    itsAttr[BOUNDPDESTROY] = Attributes::makeReal(
            "BOUNDPDESTROY",
            "The number of sigmas away from the mean "
            "at which particles are deleted. Default 10",
            boundpDestroy);

    itsAttr[BEAMHALOBOUNDARY] = Attributes::makeReal(
            "BEAMHALOBOUNDARY",
            "Defines in terms of sigma where "
            "the halo starts. Default 0.0",
            beamHaloBoundary);

    itsAttr[IDEALIZED] = Attributes::makeBool(
            "IDEALIZED",
            "Using the hard edge model for the calculation "
            "of path length. Default: false",
            false);

    itsAttr[LOGBENDTRAJECTORY] = Attributes::makeBool(
            "LOGBENDTRAJECTORY",
            "Writing the trajectory of "
            "every bend to disk. Default: false",
            false);

    itsAttr[VERSION] = Attributes::makeReal(
            "VERSION", "Version of OPAL for which input file was written", version);

    itsAttr[MEMORYDUMP] =
            Attributes::makeBool("MEMORYDUMP", "If true, write memory to SDDS file", false);

    itsAttr[HALOSHIFT] = Attributes::makeReal(
            "HALOSHIFT", "Constant parameter to shift halo value (default: 0.0)", haloShift);

    itsAttr[DELPARTFREQ] = Attributes::makeReal(
            "DELPARTFREQ",
            "The frequency to delete particles, "
            "i.e. delete when step%1 == 0. Default: 1",
            1);

    itsAttr[COMPUTEPERCENTILES] = Attributes::makeBool(
            "COMPUTEPERCENTILES",
            "Flag to control whether the 68.27 "
            "(1 sigma for normal distribution), the 95.45 (2 sigmas), "
            "the 99.73 (3 sigmas) and the 99.994 (4 sigmas) percentiles "
            "for the beam size and the normalized emittance should "
            "be computed. Default: false",
            computePercentiles);

    itsAttr[QM_MODE] = Attributes::makeString(
            "QM_MODE",
            "Storage mode for particle charge/mass. "
            "SINGLE uses one shared value per container; "
            "ATTRIBUTES stores per-particle values (accessible via per-particle "
            "`Q`/`M` views).",
            useQMAttributes ? std::string("ATTRIBUTES") : std::string("SINGLE"));

    itsAttr[AGGRESSIVE_STATE_SYNC] = Attributes::makeBool(
            "AGGRESSIVE_STATE_SYNC",
            "If true, every mutation of the shared BunchStateHandler flags "
            "(moments-dirty, unitless-positions, emitting-now) "
            "performs an MPI allreduce so that all ranks converge to the same "
            "value. Guards against rank-local divergence at the cost of an extra "
            "collective on every state change. Default: false.",
            aggressiveStateSync);

    itsAttr[LOADBALANCINGTHRESHOLD] = Attributes::makeReal(
            "LOADBALANCINGTHRESHOLD",
            "The threshold for triggering load balancing. If the ratio difference of particles in "
            "a rank exceeds this threshold, load balancing will be triggered. Default is 0.05 "
            "(5%). This threshold is only tested every `repartFreq` steps.",
            loadBalancingThreshold);

    itsAttr[LOADBALANCINGPAIRWEIGHT] = Attributes::makeReal(
            "LOADBALANCINGPAIRWEIGHT",
            "Share of the ORB load-balancing weight that follows the P3M particle-particle pair "
            "count; the rest follows the particle count. 0 balances particle counts only, 1 pair "
            "counts only. Values above 0 require the P3M field solver. Default is 0.",
            loadBalancingPairWeight);

    registerOwnership(AttributeHandler::STATEMENT);

    FileStream::setEcho(echo);
}

Option::Option(const std::string& name, Option* parent) : Action(name, parent) {
    Attributes::setBool(itsAttr[ECHO], echo);
    Attributes::setBool(itsAttr[INFO], info);
    Attributes::setBool(itsAttr[TRACE], mtrace);
    Attributes::setBool(itsAttr[WARN], warn);
    Attributes::setReal(itsAttr[SEED], seed);
    Attributes::setReal(itsAttr[PSDUMPFREQ], psDumpFreq);
    Attributes::setReal(itsAttr[STATDUMPFREQ], statDumpFreq);
    Attributes::setReal(itsAttr[CHECKPOINTFREQ], checkpointFreq);
    Attributes::setReal(itsAttr[STEPINFOFQ], stepInfoFreq);
    Attributes::setReal(itsAttr[PRINTRANKDISTRFQ], printRankDistrFreq);
    Attributes::setReal(itsAttr[REPARTFREQ], repartFreq);
    Attributes::setReal(itsAttr[MINBINEMITTED], minBinEmitted);
    Attributes::setBool(itsAttr[RHODUMP], rhoDump);
    Attributes::setBool(itsAttr[EBDUMP], ebDump);
    Attributes::setBool(itsAttr[RANKDUMP], rankDump);
    Attributes::setReal(itsAttr[AUTOPHASE], autoPhase);
    Attributes::setBool(itsAttr[ENABLEHDF5], enableHDF5);
    Attributes::setBool(itsAttr[ENABLEVTK], enableVTK);
    Attributes::setBool(itsAttr[ENABLELINEARTRANSFERMAPS], false);
    Attributes::setRealArray(
            itsAttr[ELEMENTMAPPINGGRID],
            std::vector<double>(elementMappingGrid.begin(), elementMappingGrid.end()));
    Attributes::setReal(itsAttr[ELEMENTMAPPINGREACH], elementMappingReach);
    Attributes::setBool(itsAttr[ASCIIDUMP], asciidump);
    Attributes::setReal(itsAttr[BOUNDPDESTROY], boundpDestroy);
    Attributes::setReal(itsAttr[BEAMHALOBOUNDARY], beamHaloBoundary);
    Attributes::setReal(itsAttr[VERSION], version);
    Attributes::setReal(itsAttr[HALOSHIFT], haloShift);
    Attributes::setBool(itsAttr[COMPUTEPERCENTILES], computePercentiles);
    Attributes::setString(
            itsAttr[QM_MODE], useQMAttributes ? std::string("ATTRIBUTES") : std::string("SINGLE"));
    Attributes::setBool(itsAttr[AGGRESSIVE_STATE_SYNC], aggressiveStateSync);
    Attributes::setReal(itsAttr[LOADBALANCINGTHRESHOLD], loadBalancingThreshold);
    Attributes::setReal(itsAttr[LOADBALANCINGPAIRWEIGHT], loadBalancingPairWeight);
}

Option::~Option() {}

Option* Option::clone(const std::string& name) { return new Option(name, this); }

void Option::execute() {
    if (Attributes::getBool(itsAttr[ENABLELINEARTRANSFERMAPS]))
        throw OpalException(
                "Option::execute",
                "Linear transfer maps were removed from OPALX; ENABLELINEARTRANSFERMAPS must be "
                "FALSE.");

    const auto grid = Attributes::getRealArray(itsAttr[ELEMENTMAPPINGGRID]);
    std::array<unsigned, 3> mappingGrid{};
    double mappingCells = 1.0;
    bool validGrid      = grid.size() == 3;
    for (std::size_t d = 0; validGrid && d < 3; ++d) {
        validGrid = std::isfinite(grid[d]) && grid[d] >= 1.0 && std::floor(grid[d]) == grid[d];
        mappingCells *= grid[d];
        if (validGrid) mappingGrid[d] = static_cast<unsigned>(grid[d]);
    }
    if (!validGrid || mappingCells > GeneralElementMapping::maximumCells)
        throw OpalException(
                "Option::execute",
                "ELEMENTMAPPINGGRID requires three positive integers with at most 2^24 cells "
                "in total.");
    const double mappingReach = Attributes::getReal(itsAttr[ELEMENTMAPPINGREACH]);
    if (!std::isfinite(mappingReach) || mappingReach <= 0.0)
        throw OpalException("Option::execute", "ELEMENTMAPPINGREACH must be a positive length.");

    // Store the option flags.
    echo       = Attributes::getBool(itsAttr[ECHO]);
    info       = Attributes::getBool(itsAttr[INFO]);
    mtrace     = Attributes::getBool(itsAttr[TRACE]);
    warn       = Attributes::getBool(itsAttr[WARN]);
    rhoDump    = Attributes::getBool(itsAttr[RHODUMP]);
    ebDump     = Attributes::getBool(itsAttr[EBDUMP]);
    rankDump   = Attributes::getBool(itsAttr[RANKDUMP]);
    enableHDF5 = Attributes::getBool(itsAttr[ENABLEHDF5]);
    enableVTK  = Attributes::getBool(itsAttr[ENABLEVTK]);

    elementMappingGrid  = mappingGrid;
    elementMappingReach = mappingReach;

    asciidump = Attributes::getBool(itsAttr[ASCIIDUMP]);
    version   = Attributes::getReal(itsAttr[VERSION]);
    seed      = Attributes::getReal(itsAttr[SEED]);

    haloShift                = Attributes::getReal(itsAttr[HALOSHIFT]);
    computePercentiles       = Attributes::getBool(itsAttr[COMPUTEPERCENTILES]);
    const std::string qmMode = Attributes::getString(itsAttr[QM_MODE]);
    if (qmMode == "ATTRIBUTES") {
        useQMAttributes = true;
    } else if (qmMode == "SINGLE") {
        useQMAttributes = false;
    } else {
        throw OpalException(
                "Option::execute",
                "Unsupported QM_MODE '" + qmMode + "'. Use \"SINGLE\" or \"ATTRIBUTES\".");
    }

    aggressiveStateSync = Attributes::getBool(itsAttr[AGGRESSIVE_STATE_SYNC]);

    // Set threshold and check if in the valid range [0, 1]
    loadBalancingThreshold = Attributes::getReal(itsAttr[LOADBALANCINGTHRESHOLD]);
    if (loadBalancingThreshold < 0.0 || loadBalancingThreshold > 1.0) {
        throw OpalException(
                "Option::execute",
                "LOADBALANCINGTHRESHOLD must be in the range [0, 1]. Current value: "
                        + std::to_string(loadBalancingThreshold));
    }

    loadBalancingPairWeight = Attributes::getReal(itsAttr[LOADBALANCINGPAIRWEIGHT]);
    if (!(loadBalancingPairWeight >= 0.0 && loadBalancingPairWeight <= 1.0)) {
        throw OpalException(
                "Option::execute",
                "LOADBALANCINGPAIRWEIGHT must be in the range [0, 1]. Current value: "
                        + std::to_string(loadBalancingPairWeight));
    }

    /// note: rangen is used only for the random number generator in the OPAL language
    ///       not for the distributions

    if (Options::seed == -1)
        rangen.init55(time(0));
    else
        rangen.init55(seed);

    /// \todo check if we need to do this anymore
    /*
    IpplInfo::Info->on(info);
    IpplInfo::Warn->on(warn);
    */

    /// note: rangen is used only for the random number generator in the OPAL language
    ///       not for the distributions
    if (itsAttr[SEED]) {
        seed = int(Attributes::getReal(itsAttr[SEED]));
        if (seed == -1)
            rangen.init55(time(0));
        else
            rangen.init55(seed);
    }

    if (itsAttr[PSDUMPFREQ]) {
        psDumpFreq = int(Attributes::getReal(itsAttr[PSDUMPFREQ]));
        if (psDumpFreq < 0) {
            throw OpalException("Option::execute", "PSDUMPFREQ must be non-negative.");
        }
    }

    if (itsAttr[STATDUMPFREQ]) {
        statDumpFreq = int(Attributes::getReal(itsAttr[STATDUMPFREQ]));
        if (statDumpFreq < 0) {
            throw OpalException("Option::execute", "STATDUMPFREQ must be non-negative.");
        }
    }

    if (itsAttr[CHECKPOINTFREQ]) {
        checkpointFreq = int(Attributes::getReal(itsAttr[CHECKPOINTFREQ]));
        if (checkpointFreq < 0) {
            throw OpalException("Option::execute", "CHECKPOINTFREQ must be non-negative.");
        }
    }

    if (itsAttr[STEPINFOFQ]) {
        stepInfoFreq = int(Attributes::getReal(itsAttr[STEPINFOFQ]));
        if (stepInfoFreq < 0) stepInfoFreq = 0;
    }

    if (itsAttr[PRINTRANKDISTRFQ]) {
        printRankDistrFreq = int(Attributes::getReal(itsAttr[PRINTRANKDISTRFQ]));
        if (printRankDistrFreq < 0) printRankDistrFreq = 0;
    }

    if (itsAttr[REPARTFREQ]) {
        repartFreq = int(Attributes::getReal(itsAttr[REPARTFREQ]));
    }

    if (itsAttr[MINBINEMITTED]) {
        minBinEmitted = int(Attributes::getReal(itsAttr[MINBINEMITTED]));
    }

    if (itsAttr[AUTOPHASE]) {
        autoPhase = int(Attributes::getReal(itsAttr[AUTOPHASE]));
    }

    if (itsAttr[BOUNDPDESTROY]) {
        /*
         * Historically, BOUNDPDESTROYFQ was used as a positive frequency and
         * clamped to values >= 1 here. In OPALX, the same parameter now also
         * encodes the N-sigma boundary used in ParticleContainer::
         * markParticlesOutside(N). A value <= 0 is treated as "disabled"
         * by markParticlesOutside, so we must *not* clamp it to 1 anymore.
         */
        boundpDestroy = Attributes::getReal(itsAttr[BOUNDPDESTROY]);
    }

    if (itsAttr[BEAMHALOBOUNDARY]) {
        beamHaloBoundary = Attributes::getReal(itsAttr[BEAMHALOBOUNDARY]);
    } else {
        beamHaloBoundary = 0;
    }

    // These options are still accepted so that old input files read unchanged, but nothing
    // uses them. Say so once per run and option.
    static std::set<std::string> reported;
    for (const int attr :
         {PSDUMPEACHTURN, PSDUMPFRAME, REMOTEPARTDEL, CSRDUMP, IDEALIZED, LOGBENDTRAJECTORY,
          DELPARTFREQ, SPTDUMPFREQ, SCSOLVEFREQ, MTSSUBSTEPS, MINSTEPFORREBIN, REBINFREQ, NUMBLOCKS,
          RECYCLEBLOCKS, NLHS, CZERO, RNGTYPE, CLOTUNEONLY, MEMORYDUMP}) {
        const Attribute& attribute = itsAttr[attr];
        if (!attribute.defaultUsed() && reported.insert(attribute.getName()).second) {
            *gmsg << level1 << "* OPTION " << attribute.getName() << " has no effect in OPALX."
                  << endl;
        }
    }

    // Set message flags.
    FileStream::setEcho(echo);

    if (Attributes::getBool(itsAttr[TELL])) {
        *gmsg << "\nCurrent settings of options:\n" << *this << endl;
    }

    Option* main = dynamic_cast<Option*>(OpalData::getInstance()->find("OPTION"));
    if (main) {
        main->update(itsAttr);
    }
}

void Option::update(const std::vector<Attribute>& othersAttributes) {
    for (int i = 0; i < SIZE; ++i) {
        itsAttr[i] = othersAttributes[i];
    }
}
