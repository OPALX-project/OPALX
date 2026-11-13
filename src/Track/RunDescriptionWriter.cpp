//
// Text dump of a RunDescription, to compare what two OPALX versions read from the same input.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#include "Track/RunDescriptionWriter.h"

#include <cstdio>
#include <ostream>
#include <string>
#include <type_traits>

namespace opalx::run {

    namespace {

        std::string number(double value) {
            char buffer[32];
            std::snprintf(buffer, sizeof buffer, "%.17g", value);
            return buffer;
        }

        std::string text(const AttributeValue& value) {
            return std::visit(
                    [](const auto& v) -> std::string {
                        using T = std::decay_t<decltype(v)>;
                        if constexpr (std::is_same_v<T, double>) {
                            return number(v);
                        } else if constexpr (std::is_same_v<T, bool>) {
                            return v ? "true" : "false";
                        } else if constexpr (std::is_same_v<T, std::string>) {
                            return "\"" + v + "\"";
                        } else {
                            std::string out = "[";
                            for (std::size_t i = 0; i < v.size(); ++i) {
                                if (i > 0) out += ", ";
                                using E = std::decay_t<decltype(v[i])>;
                                if constexpr (std::is_same_v<E, std::string>) {
                                    out += "\"" + v[i] + "\"";
                                } else if constexpr (std::is_same_v<E, double>) {
                                    out += number(v[i]);
                                } else {
                                    out += v[i] ? "true" : "false";
                                }
                            }
                            return out + "]";
                        }
                    },
                    value);
        }

        void writeObject(const ObjectDescription& object, std::ostream& os) {
            os << "keyword = " << object.keyword << "\n";
            os << "name = " << object.name << "\n";
            for (const AttributeEntry& a : object.attributes) {
                os << "attr " << a.name << " (" << a.type << ", "
                   << (a.defaultUsed ? "default" : "set") << ") = ";
                if (!a.error.empty()) {
                    os << "<error: " << a.error << ">";
                } else if (!a.hasValue) {
                    os << "<none>";
                } else {
                    os << text(a.value);
                }
                os << "\n";
            }
        }

        template <class T>
        std::string list(const std::vector<T>& values) {
            std::string out;
            for (std::size_t i = 0; i < values.size(); ++i) {
                if (i > 0) out += ", ";
                if constexpr (std::is_same_v<T, double>) {
                    out += number(values[i]);
                } else if constexpr (std::is_same_v<T, std::string>) {
                    out += values[i];
                } else {
                    out += std::to_string(values[i]);
                }
            }
            return out;
        }

        void writeOptions(const OptionValues& o, std::ostream& os) {
            const auto flag = [](bool v) {
                return v ? "true" : "false";
            };
            os << "echo = " << flag(o.echo) << "\n"
               << "info = " << flag(o.info) << "\n"
               << "mtrace = " << flag(o.mtrace) << "\n"
               << "warn = " << flag(o.warn) << "\n"
               << "fft2d5BatchSize = " << o.fft2d5BatchSize << "\n"
               << "seed = " << o.seed << "\n"
               << "psDumpFreq = " << o.psDumpFreq << "\n"
               << "statDumpFreq = " << o.statDumpFreq << "\n"
               << "checkpointFreq = " << o.checkpointFreq << "\n"
               << "stepInfoFreq = " << o.stepInfoFreq << "\n"
               << "printRankDistrFreq = " << o.printRankDistrFreq << "\n"
               << "useQMAttributes = " << flag(o.useQMAttributes) << "\n"
               << "repartFreq = " << o.repartFreq << "\n"
               << "minBinEmitted = " << o.minBinEmitted << "\n"
               << "rhoDump = " << flag(o.rhoDump) << "\n"
               << "ebDump = " << flag(o.ebDump) << "\n"
               << "rankDump = " << flag(o.rankDump) << "\n"
               << "autoPhase = " << o.autoPhase << "\n"
               << "enableHDF5 = " << flag(o.enableHDF5) << "\n"
               << "enableVTK = " << flag(o.enableVTK) << "\n"
               << "elementMappingGrid = " << o.elementMappingGrid[0] << ", "
               << o.elementMappingGrid[1] << ", " << o.elementMappingGrid[2] << "\n"
               << "elementMappingReach = " << number(o.elementMappingReach) << "\n"
               << "asciidump = " << flag(o.asciidump) << "\n"
               << "boundpDestroy = " << number(o.boundpDestroy) << "\n"
               << "beamHaloBoundary = " << number(o.beamHaloBoundary) << "\n"
               << "version = " << o.version << "\n"
               << "haloShift = " << number(o.haloShift) << "\n"
               << "computePercentiles = " << flag(o.computePercentiles) << "\n"
               << "aggressiveStateSync = " << flag(o.aggressiveStateSync) << "\n"
               << "loadBalancingThreshold = " << number(o.loadBalancingThreshold) << "\n"
               << "loadBalancingPairWeight = " << number(o.loadBalancingPairWeight) << "\n";
        }

    }  // namespace

    void write(const RunDescription& d, std::ostream& os) {
        const auto flag = [](bool v) {
            return v ? "true" : "false";
        };
        os << "# OPALX run description, RUN " << d.context.runIndex << "\n";

        os << "\n[context]\n"
           << "inputBasename = " << d.context.inputBasename << "\n"
           << "restart = " << flag(d.context.restart) << "\n"
           << "restartFileName = " << d.context.restartFileName << "\n"
           << "restartStep = " << d.context.restartStep << "\n"
           << "openMode = " << d.context.openMode << "\n"
           << "p0 = " << number(d.context.p0) << "\n"
           << "globalPhaseShift = " << number(d.context.globalPhaseShift) << "\n";

        os << "\n[options]\n";
        writeOptions(d.options, os);

        const TrackDescription& t = d.track;
        os << "\n[track]\n"
           << "beams = " << list(t.beamNames) << "\n"
           << "dT = " << list(t.dT) << "\n"
           << "maxSteps = " << list(t.maxSteps) << "\n"
           << "zStop = " << list(t.zStop) << "\n"
           << "zStopGiven = " << flag(t.zStopGiven) << "\n"
           << "zStart = " << number(t.zStart) << "\n"
           << "t0 = " << number(t.t0) << "\n"
           << "stepsPerTurn = " << t.stepsPerTurn << "\n"
           << "kineticEnergyStopGeV = " << number(t.kineticEnergyStopGeV) << "\n"
           << "hasInitialOrbit = " << flag(t.hasInitialOrbit) << "\n";

        os << "\n[run]\n";
        writeObject(t.run, os);

        for (const BeamDescription& beam : d.beams) {
            os << "\n[beam " << beam.object.name << "]\n";
            os << "reference.mass = " << number(beam.reference.mass) << "\n"
               << "reference.charge = " << number(beam.reference.charge) << "\n"
               << "reference.momentum = " << number(beam.reference.momentum) << "\n"
               << "reference.gamma = " << number(beam.reference.gamma) << "\n"
               << "chargePerParticle = " << number(beam.chargePerParticle) << "\n"
               << "massPerParticle = " << number(beam.massPerParticle) << "\n"
               << "emissionSourceList = " << beam.emissionSourceList << "\n";
            writeObject(beam.object, os);
        }
        for (const ObjectDescription& source : d.emissionSources) {
            os << "\n[emissionsource " << source.name << "]\n";
            writeObject(source, os);
        }
        for (const ObjectDescription& distribution : d.distributions) {
            os << "\n[distribution " << distribution.name << "]\n";
            writeObject(distribution, os);
        }
        if (d.fieldSolver) {
            os << "\n[fieldsolver " << d.fieldSolver->name << "]\n";
            writeObject(*d.fieldSolver, os);
        }
        if (d.binning) {
            os << "\n[binning " << d.binning->name << "]\n";
            writeObject(*d.binning, os);
        }
        if (d.boundaryGeometry) {
            os << "\n[boundarygeometry " << d.boundaryGeometry->name << "]\n";
            writeObject(*d.boundaryGeometry, os);
        }

        os << "\n[line " << d.line.name << "]\n"
           << "topology = " << d.line.topology << "\n"
           << "elements = " << d.line.elements.size() << "\n";
        for (const ElementDescription& element : d.line.elements) {
            os << "\n[element " << element.occurrence << " " << element.object.name << "]\n"
               << "elementType = " << element.elementType << "\n"
               << "reflected = " << flag(element.reflected) << "\n";
            writeObject(element.object, os);
        }
    }

    std::string dumpFileName(const RunDescription& description) {
        return description.context.inputBasename + "_run"
               + std::to_string(description.context.runIndex) + ".txt";
    }

}  // namespace opalx::run
