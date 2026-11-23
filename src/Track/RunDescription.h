//
// RunDescription: the values one RUN command tracks with, copied out of the parsed objects.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#ifndef OPALX_RUN_DESCRIPTION_H
#define OPALX_RUN_DESCRIPTION_H

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "Algorithms/PartData.h"

namespace opalx::run {

    /// An evaluated attribute value. Attribute types without a plain value (range, reference,
    /// token list) are kept as their text image.
    using AttributeValue = std::variant<
            double, bool, std::string, std::vector<double>, std::vector<bool>,
            std::vector<std::string>>;

    /**
     * @brief One attribute of a parsed object, evaluated when the RUN command executes.
     */
    struct AttributeEntry {
        std::string name;          ///< Attribute name, as in the input.
        std::string type;          ///< Attribute type name ("real", "string", "real array", ...).
        bool hasValue    = false;  ///< False if the attribute has neither a value nor a default.
        bool defaultUsed = true;   ///< True if the input did not set it.
        AttributeValue value;      ///< The evaluated value, if hasValue.
        std::string error;         ///< Non-empty if evaluating the attribute failed.
    };

    /**
     * @brief A parsed input object: its keyword, its name and all its attributes.
     */
    struct ObjectDescription {
        std::string keyword;  ///< Keyword the object was made from, e.g. "QUADRUPOLE".
        std::string name;     ///< Name given in the input.
        std::vector<AttributeEntry> attributes;
    };

    /**
     * @brief One EMISSIONSOURCE of a beam and the distribution it names.
     */
    struct SourceDescription {
        std::string name;
        std::string distribution;        ///< DISTRIBUTION attribute of the source.
        bool distributionFound = false;  ///< False if no DISTRIBUTION has that name.
        std::string distributionType;    ///< TYPE of the distribution, upper case.
        bool emitted = false;            ///< EMITTED of the distribution.
        double t0    = 0;                ///< T0 of the source (s).
    };

    struct BeamDescription {
        ObjectDescription object;
        /// The reference particle (OPALX units: eV, eV/c, elementary charges), from
        /// Beam::getReference(). The writer writes its mass, charge, momentum and gamma.
        PartData reference;
        double chargePerParticle = 0;  ///< Macro charge per particle (C).
        double massPerParticle   = 0;  ///< Macro mass per particle (GeV/c^2).
        std::string emissionSourceList;

        // Values read through Beam's own functions, for validate(). The text writer leaves them
        // out: the attributes above already hold the input they come from.
        std::string particleName;
        bool isPhoton = false;
        std::string daughterBeam;  ///< DAUGHTERBEAM, empty if not set.
        std::vector<std::string> globalProcesses;
        bool hasExplicitEnergy  = false;  ///< GAMMA, ENERGY or PC given.
        bool hasEmissionSources = false;  ///< SOURCES given.
        bool hasPolarization    = false;
        std::optional<double> zStop;             ///< BEAM ZSTOP.
        bool untilMaxSteps = false;              ///< BEAM UNTILMAXSTEPS.
        std::string emissionSourceListError;     ///< Why the SOURCES name could not be read.
        bool emissionSourceListFound = false;    ///< An EMISSIONSOURCELIST of that name exists.
        std::vector<SourceDescription> sources;  ///< The sources of that list, in order.

        // Read through Beam's own functions for buildTracking(); not written either.
        double charge        = 0;                ///< Beam::getCharge() (elementary charges).
        double mass          = 0;                ///< Beam::getMass() (GeV/c^2).
        std::size_t numAlloc = 0;                ///< NALLOC.
        std::vector<double> polarization;        ///< Beam::getPolarization(), {0, 0, 0} if not set.
        bool hasReferenceDirection = false;      ///< REFDIRECTION given.
        std::vector<double> referenceDirection;  ///< REFDIRECTION as a unit vector.
        std::string banner;  ///< The text Beam::print() writes, for the RUN log.
    };

    /**
     * @brief One element occurrence in the tracked line, in tracking order.
     */
    struct ElementDescription {
        ObjectDescription object;        ///< Attributes of the input element this occurrence uses.
        std::string elementType;         ///< Runtime element type, e.g. "MULTIPOLE".
        std::size_t occurrence = 0;      ///< Position in the tracked line, starting at 0.
        bool reflected         = false;  ///< Member written as -NAME.
    };

    struct LineDescription {
        std::string name;
        std::vector<ElementDescription> elements;
        std::array<double, 3> origin{0.0, 0.0, 0.0};          ///< Start position (m). Not written.
        std::array<double, 4> direction{1.0, 0.0, 0.0, 0.0};  ///< Start direction. Not written.
        /// The time dependences and end-field models of the input, by name; elements name them.
        std::vector<ObjectDescription> models;
    };

    /**
     * @brief The TRACK settings and the RUN command attributes.
     */
    struct TrackDescription {
        std::string line;  ///< LINE of the TRACK command. Not written; see LineDescription.
        std::vector<std::string> beamNames;
        std::vector<double> dT;
        std::vector<unsigned long long> maxSteps;
        std::vector<double> zStop;
        bool zStopGiven = false;
        double zStart   = 0;
        double t0       = 0;
        ObjectDescription run;        ///< Attributes of the RUN command.
        std::string fieldSolverType;  ///< TYPE of the RUN FIELDSOLVER, if found. Not written.
    };

    /**
     * @brief The OPTION values in effect when RUN executes (a copy of the Options globals).
     */
    struct OptionValues {
        bool echo = false, info = false, mtrace = false, warn = false;
        int seed       = 0;
        int psDumpFreq = 0, statDumpFreq = 0, checkpointFreq = 0, stepInfoFreq = 0;
        int printRankDistrFreq = 0;
        bool useQMAttributes   = false;
        int repartFreq = 0, minBinEmitted = 0;
        bool rhoDump = false, ebDump = false, rankDump = false;
        int autoPhase   = 0;
        bool enableHDF5 = false, enableVTK = false;
        std::array<unsigned, 3> elementMappingGrid{};
        double elementMappingReach = 0;
        bool asciidump             = false;
        double boundpDestroy = 0, beamHaloBoundary = 0;
        int version             = 0;
        double haloShift        = 0;
        bool computePercentiles = false, aggressiveStateSync = false;
        double loadBalancingThreshold = 0, loadBalancingPairWeight = 0;
    };

    /**
     * @brief One DUMPEMFIELDS command: the file and the Cartesian space-time grid on which the
     *        external fields of the tracked line are written.
     */
    struct EMFieldDumpDescription {
        /// FILE_NAME; a relative name is placed in the auxiliary output directory.
        std::string fileName;
        std::array<double, 4> start{};    ///< X_START, Y_START, Z_START (m) and T_START (ns).
        std::array<double, 4> spacing{};  ///< DX, DY, DZ (m) and DT (ns).
        std::array<int, 4> steps{};       ///< X_STEPS, Y_STEPS, Z_STEPS, T_STEPS: grid points.
    };

    /**
     * @brief Values of the run that do not belong to one input object.
     */
    struct RunContext {
        int runIndex = 0;  ///< 1 for the first RUN of the input file, 2 for the second, ...
        std::string inputBasename;
        std::string inputFileName;  ///< The input file as given on the command line.
        bool restart = false;
        std::string restartFileName;
        std::string openMode;
        double p0 = 0;  ///< The input variable P0, used to scale magnet strengths.
    };

    /**
     * @brief Everything one RUN command tracks with, as plain values.
     *
     * Filled by freeze() when RUN executes. It holds no pointers into the parsed objects.
     */
    struct RunDescription {
        RunContext context;
        OptionValues options;
        TrackDescription track;
        LineDescription line;
        std::vector<BeamDescription> beams;
        std::vector<ObjectDescription> emissionSources;
        std::vector<ObjectDescription> distributions;
        std::optional<ObjectDescription> fieldSolver;
        std::optional<ObjectDescription> binning;
        /// The DUMPEMFIELDS commands executed before the RUN, in input order.
        std::vector<EMFieldDumpDescription> emFieldDumps;
    };

}  // namespace opalx::run

#endif  // OPALX_RUN_DESCRIPTION_H
