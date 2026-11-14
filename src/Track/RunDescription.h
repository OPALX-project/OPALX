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

namespace opalx::run {

    /// An evaluated attribute value. Attribute types without a plain value (place, range,
    /// reference, token list, table line) are kept as their text image.
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
     * @brief The reference particle a BEAM defines (OPALX units: eV, eV/c, elementary charges).
     */
    struct ReferenceParticle {
        double mass     = 0;
        double charge   = 0;
        double momentum = 0;
        double gamma    = 0;
        double beta     = 0;  ///< Not written by the text writer.
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
        ReferenceParticle reference;
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
        std::vector<double> tuneInitial;         ///< TUNEINITIAL.
        std::string emissionSourceListError;     ///< Why the SOURCES name could not be read.
        bool emissionSourceListFound = false;    ///< An EMISSIONSOURCELIST of that name exists.
        std::vector<SourceDescription> sources;  ///< The sources of that list, in order.
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
        std::string topology;  ///< "LINE" or "RING".
        std::vector<ElementDescription> elements;
        /// Circumference of a RING: the sum of the arc lengths of the built elements
        /// (opalx::run::ringLength()). Not written.
        double length = 0;
        std::array<double, 3> origin{0.0, 0.0, 0.0};          ///< Start position (m). Not written.
        std::array<double, 4> direction{1.0, 0.0, 0.0, 0.0};  ///< Start direction. Not written.
        /// The time dependences and end-field models of the input, by name; elements name them.
        std::vector<ObjectDescription> models;
    };

    /**
     * @brief The TRACK settings (from Track::block) and the RUN command attributes.
     */
    struct TrackDescription {
        std::vector<std::string> beamNames;
        std::vector<double> dT;
        std::vector<unsigned long long> maxSteps;
        std::vector<double> zStop;
        bool zStopGiven             = false;
        double zStart               = 0;
        double t0                   = 0;
        int stepsPerTurn            = 0;
        double kineticEnergyStopGeV = 0;
        bool hasInitialOrbit        = false;
        ObjectDescription run;        ///< Attributes of the RUN command.
        std::string fieldSolverType;  ///< TYPE of the RUN FIELDSOLVER, if found. Not written.
    };

    /**
     * @brief The OPTION values in effect when RUN executes (a copy of the Options globals).
     */
    struct OptionValues {
        bool echo = false, info = false, mtrace = false, warn = false;
        int fft2d5BatchSize = 0, seed = 0;
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
     * @brief Values of the run that do not belong to one input object.
     */
    struct RunContext {
        int runIndex = 0;  ///< 1 for the first RUN of the input file, 2 for the second, ...
        std::string inputBasename;
        bool restart = false;
        std::string restartFileName;
        int restartStep = 0;
        std::string openMode;
        double p0               = 0;  ///< The input variable P0, used to scale magnet strengths.
        double globalPhaseShift = 0;
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
    };

}  // namespace opalx::run

#endif  // OPALX_RUN_DESCRIPTION_H
