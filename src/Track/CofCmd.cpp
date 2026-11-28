// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Track/CofCmd.h"
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <locale>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <vector>
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/ClosedOrbitSolver.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/LinearMapEigenAnalysis.h"
#include "Algorithms/PartData.h"
#include "Algorithms/TpsClosedOrbitMap.h"
#include "Attributes/Attributes.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "Elements/OpalBeamline.h"
#include "PartBunch/CartesianDomainConfig.h"
#include "Structure/Beam.h"
#include "Utilities/OpalException.h"
#include "Utility/Inform.h"

extern Inform* gmsg;
namespace {
    enum Outer {
        LINE,
        BEAM,
        DT,
        MAXSTEPS,
        T0,
        TIMEINTEGRATOR,
        MAXPATH,
        SECTION,
        GEOMTOL,
        ANGLETOL,
        METHOD,
        DIMENSION,
        X,
        PX,
        Y,
        PY,
        MAXIT,
        XTOL,
        PTOL,
        FDSTEP,
        SCALES,
        DAMPING,
        JACOBIAN,
        OUTPUT,
        SIZE
    };
    void require(bool ok, const std::string& message) {
        if (!ok) throw OpalException("COF", message);
    }
    // Small streaming writer for the COF export: JSON escaping, finite numbers,
    // classic-locale round-trip precision and two-space indentation.
    class CofJsonWriter {
        std::ostringstream out;
        std::vector<bool> first;

        void quoted(const std::string& text) {
            constexpr char hex[] = "0123456789abcdef";
            out << '"';
            for (unsigned char c : text) {
                if (c == '"' || c == '\\')
                    out << '\\' << c;
                else if (c < 0x20)
                    out << "\\u00" << hex[c >> 4] << hex[c & 15];
                else
                    out << c;
            }
            out << '"';
        }
        void prefix(const char* key) {
            if (!first.empty()) {
                if (!first.back()) out << ',';
                out << '\n' << std::string(2 * first.size(), ' ');
                first.back() = false;
            }
            if (key) {
                quoted(key);
                out << ": ";
            }
        }
        void finite(double value) {
            require(std::isfinite(value), "Cannot serialize a nonfinite COF result.");
            out << value;
        }
        template <class F>
        void container(const char* key, char open, char close, F contents) {
            prefix(key);
            out << open;
            first.push_back(true);
            contents();
            const bool empty = first.back();
            first.pop_back();
            if (!empty) out << '\n' << std::string(2 * first.size(), ' ');
            out << close;
        }

    public:
        CofJsonWriter() {
            out.imbue(std::locale::classic());
            out << std::setprecision(std::numeric_limits<double>::max_digits10);
        }
        template <class F>
        void object(const char* key, F contents) {
            container(key, '{', '}', contents);
        }
        template <class F>
        void array(const char* key, F contents) {
            container(key, '[', ']', contents);
        }
        void text(const char* key, const std::string& value) {
            prefix(key);
            quoted(value);
        }
        void number(const char* key, double value) {
            prefix(key);
            finite(value);
        }
        void boolean(const char* key, bool value) {
            prefix(key);
            out << (value ? "true" : "false");
        }
        void null(const char* key) {
            prefix(key);
            out << "null";
        }
        template <class Values>
        void numbers(const char* key, const Values& values) {
            prefix(key);
            out << '[';
            bool comma = false;
            for (double value : values) {
                if (comma) out << ", ";
                finite(value);
                comma = true;
            }
            out << ']';
        }
        std::string str() const { return out.str() + '\n'; }
    };
    std::string lowerCase(std::string text) {
        std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) {
            return std::tolower(c);
        });
        return text;
    }
    // Protect the top-level deck and every disk input opened during this run.
    // Canonical paths cover spelling/symlink aliases; equivalent also covers hard links.
    void validateOutput(const std::string& output) {
        namespace fs = std::filesystem;
        std::error_code error;
        const auto destination = fs::weakly_canonical(output, error);
        require(!error, "Cannot resolve COF OUTPUT " + output + ": " + error.message());
        const bool destinationExists = fs::exists(output, error);
        require(!error, "Cannot inspect COF OUTPUT " + output + ": " + error.message());
        auto inputs          = OpalData::getInstance()->getInputFiles();
        const auto mainInput = OpalData::getInstance()->getInputFn();
        if (!mainInput.empty()) inputs.push_back(mainInput);
        for (const auto& input : inputs) {
            const auto source = fs::weakly_canonical(input, error);
            require(!error, "Cannot resolve COF input " + input + ": " + error.message());
            const auto message = "COF OUTPUT resolves to an input file: " + output;
            require(destination != source, message);
            if (destinationExists) {
                const bool sourceExists = fs::exists(input, error);
                require(!error, "Cannot inspect COF input " + input + ": " + error.message());
                if (sourceExists) {
                    const bool same = fs::equivalent(output, input, error);
                    require(!error, "Cannot compare COF OUTPUT " + output + " with input " + input
                                            + ": " + error.message());
                    require(!same, message);
                }
            }
        }
    }
    // Keep the previous result intact until a complete replacement is closed.
    // Exclusive creation applies the umask without changing process-global state;
    // the short sibling name does not consume the destination basename's NAME_MAX.
    void writeResult(const std::string& output, const std::string& json) {
        namespace fs    = std::filesystem;
        const auto fail = [&](const char* operation, int error) {
            throw OpalException(
                    "COF",
                    std::string(operation) + " COF OUTPUT " + output + ": " + std::strerror(error));
        };
        const fs::path directory = fs::path(output).parent_path();
        std::random_device random;
        fs::path temporary;
        int fd = -1;
        for (unsigned attempt = 0; attempt < 64; ++attempt) {
            std::ostringstream name;
            name << ".cof-" << std::hex << std::setfill('0') << std::setw(8)
                 << static_cast<std::uint32_t>(random());
            temporary = directory / name.str();
            fd        = ::open(temporary.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0666);
            if (fd >= 0) break;
            if (errno != EEXIST) fail("Cannot create temporary", errno);
        }
        require(fd >= 0, "Cannot create a unique temporary COF OUTPUT for " + output);
        struct Cleanup {
            const fs::path& path;
            int descriptor;
            bool remove = true;
            ~Cleanup() {
                if (descriptor >= 0) ::close(descriptor);
                std::error_code ignored;
                if (remove) fs::remove(path, ignored);
            }
        } cleanup{temporary, fd};
        struct stat created;
        if (::fstat(fd, &created) != 0) fail("Cannot inspect temporary", errno);
        // Write privately, then restore the creation mode or existing destination mode.
        if (::fchmod(fd, 0600) != 0) fail("Cannot secure temporary", errno);
        std::unique_ptr<FILE, decltype(&std::fclose)> file(::fdopen(fd, "w"), &std::fclose);
        if (!file) fail("Cannot open temporary", errno);
        cleanup.descriptor = -1;  // The FILE now owns the descriptor.
        if (std::fwrite(json.data(), 1, json.size(), file.get()) != json.size())
            fail("Error writing", errno);
        if (std::fflush(file.get()) != 0) fail("Error flushing", errno);
        mode_t mode = created.st_mode & 0777;
        struct stat existing;
        if (::stat(output.c_str(), &existing) == 0) {
            // chown can clear permission bits, so restore the mode afterwards.
            if (existing.st_gid != created.st_gid
                && ::fchown(fd, static_cast<uid_t>(-1), existing.st_gid) != 0)
                fail("Cannot preserve group for", errno);
            mode = existing.st_mode & 0777;
        } else if (errno != ENOENT)
            fail("Cannot inspect existing", errno);
        if (::fchmod(fd, mode) != 0) fail("Cannot set permissions for", errno);
        if (std::fclose(file.release()) != 0) fail("Error closing", errno);
        validateOutput(output);  // Recheck aliases immediately before publication.
        std::error_code error;
        fs::rename(temporary, output, error);
        require(!error, "Cannot replace COF OUTPUT " + output + ": " + error.message());
        cleanup.remove = false;
    }
    struct Report {
        std::string details;
        std::string optics;
    };
    /** Squared distance to a finite design centreline, independent of aperture.
     * Straight/RBEND bodies use their box axis; SBEND bodies use the circular
     * design arc including both endpoints. Selecting the nearest segment avoids
     * applying the infinite longitudinal slab of an opposite ring arm. This is
     * a local-orbit neighbourhood assignment, not a material/solid model for
     * intersecting beam pipes. Equal-distance segments are all checked.
     */
    double bodyDistance2(const ElementBase& element, const Vector_t<double, 3>& r) {
        const auto& geometry = element.getGeometry();
        const double length  = geometry.getElementLength();
        const double h       = geometry.getCurvature();
        if (geometry.kind() != GeometryKind::SBend || std::abs(h) <= 1e-15) {
            const double dz = r(2) - std::clamp(r(2), 0.0, length);
            return r(0) * r(0) + r(1) * r(1) + dz * dz;
        }
        const double angle = h * length;
        const double ex = (std::cos(angle) - 1) / h, ez = std::sin(angle) / h;
        double distance = std::min(
                r(0) * r(0) + r(2) * r(2), (r(0) - ex) * (r(0) - ex) + (r(2) - ez) * (r(2) - ez));
        double phase = std::atan2(h * r(2), 1 + h * r(0));
        if (phase * h < 0) phase += std::copysign(2 * std::acos(-1.0), h);
        if (phase / h <= length) {
            const double radial = std::hypot(r(0) + 1 / h, r(2)) - std::abs(1 / h);
            distance            = std::min(distance, radial * radial);
        }
        return distance + r(1) * r(1);
    }
    unsigned count(double value, const char* name) {
        require(std::isfinite(value) && value >= 1 && value <= 1000000000
                        && std::floor(value) == value,
                std::string(name) + " must be a positive integer <= 1e9.");
        return static_cast<unsigned>(value);
    }
    void positive(double value, const char* name) {
        require(std::isfinite(value) && value > 0,
                std::string(name) + " must be finite and positive.");
    }
    struct Controls {
        OneTurnMap::Settings tracking;
        ExternalFieldRayTracker::IntegrationMethod integrator;
        CoordinateSystemTrafo section;
        double geometryTolerance, angleTolerance;
    };
    struct RunSettings {
        ClosedOrbitSolver::Coordinates initial{};
        ClosedOrbitSolver::Settings solver;
        std::string output;
        std::string commandName;
    };

    class LatticeVisitor : public DefaultVisitor {
        OpalBeamline& lattice;
        PartBunch_t& bunch;
        std::string owner;
        void applyDefault(const ElementBase& element) override {
            if (!element.getName().empty() && element.getName().front() == '#') return;
            const auto type = element.getType();
            // Diagnostics are passive for scalar optics. Keep PROBE outside
            // magnetic placement and do not construct sinks or histories while
            // COF tracks trial rays. The same ring can then serve COF and TRACK
            // without trial crossings contaminating or overwriting real output.
            if (type == ElementType::MONITOR || type == ElementType::PROBE) return;
            static constexpr std::array acceptedTypes{
                    ElementType::DRIFT,           ElementType::MARKER, ElementType::MULTIPOLE,
                    ElementType::SBEND,           ElementType::RBEND,  ElementType::CYCLOTRONSECTOR,
                    ElementType::SCALINGFFAMAGNET};
            require(std::find(acceptedTypes.begin(), acceptedTypes.end(), type)
                            != acceptedTypes.end(),
                    "Unsupported element in static magnetic COF: " + element.getName());
            require(element.getBeamlineMembership().ownerName == owner
                            && element.getBeamlineTopology() == BeamlineTopology::RING,
                    "Ambiguous ring membership for " + element.getName());
            const auto before = lattice.getElements();
            lattice.visit(element, *this, bunch);
            for (const auto& added : lattice.getElements())
                if (!before.count(added) && type != ElementType::MARKER) ordered.push_back(added);
        }

    public:
        std::vector<std::shared_ptr<ElementBase>> ordered;
        LatticeVisitor(
                const Beamline& line, OpalBeamline& runtime, PartBunch_t& reference,
                std::string name)
            : DefaultVisitor(line, false, false),
              lattice(runtime),
              bunch(reference),
              owner(std::move(name)) {}
        void visitFlaggedElmPtr(const FlaggedElmPtr& pointer) override {
            require(!pointer.getReflectionFlag(), "Reflected COF members are unsupported.");
            DefaultVisitor::visitFlaggedElmPtr(pointer);
        }
    };

    // Aperture ownership is independent of field-candidate selection. Conventional
    // bodies use the nearest finite centreline; circular FFA/cyclotron elements
    // own their angular windows. Every selected element uses its native loss check.
    void checkApertures(
            OpalBeamline& lattice, const std::set<std::shared_ptr<ElementBase>>& elements,
            const ExternalFieldRayTracker::State& ray) {
        double nearest = std::numeric_limits<double>::infinity();
        for (const auto& element : elements) {
            if (element->getType() == ElementType::CYCLOTRONSECTOR
                || element->getType() == ElementType::SCALINGFFAMAGNET
                || element->getGeometry().getElementLength() <= 0)
                continue;
            nearest = std::min(
                    nearest,
                    bodyDistance2(*element, lattice.transformToLocalCS(element, ray.position)));
        }
        for (const auto& element : elements) {
            if (element->getType() != ElementType::CYCLOTRONSECTOR
                && element->getType() != ElementType::SCALINGFFAMAGNET) {
                if (element->getGeometry().getElementLength() <= 0) continue;
                const double distance =
                        bodyDistance2(*element, lattice.transformToLocalCS(element, ray.position));
                // Floating-point tie allowance in squared metres;
                // not an aperture or Newton convergence tolerance.
                if (distance > nearest
                                       + 64 * std::numeric_limits<double>::epsilon()
                                                 * std::max(1.0, nearest))
                    continue;
            }
            require(!element->isReferenceParticleLost(
                            lattice.transformToLocalCS(element, ray.position),
                            lattice.rotateToLocalCS(element, ray.momentum), ray.time),
                    "Aperture/material loss in " + element->getName());
        }
    }
    void validateClosure(
            OpalBeamline& lattice, const std::vector<std::shared_ptr<ElementBase>>& ordered,
            const Controls& c, const std::string& name, const std::string& beamName,
            std::ostringstream& out) {
        const auto first = lattice.getNominalEntryTransform(ordered.front());
        const auto last  = lattice.getNominalExitTransform(ordered.back());
        const Vector_t<double, 3> zero(0.0);
        const Vector_t<double, 3> positionDifference =
                first.transformFrom(zero) - last.transformFrom(zero);
        const double positionError = euclidean_norm(positionDifference);
        double angleError          = 0;
        for (unsigned d = 0; d < 3; ++d) {
            Vector_t<double, 3> axis(0.0);
            axis[d]      = 1;
            const auto a = first.rotateFrom(axis), b = last.rotateFrom(axis);
            const Vector_t<double, 3> axisCross = cross(a, b);
            angleError = std::max(angleError, std::atan2(euclidean_norm(axisCross), dot(a, b)));
        }
        out << std::setprecision(17) << "COF ring=" << name << " beam=" << beamName
            << "\nNominal closure position_m=" << positionError << " angle_rad=" << angleError
            << '\n';
        require(positionError <= c.geometryTolerance && angleError <= c.angleTolerance,
                out.str() + "Nominal ring geometry does not close within GEOMTOL/ANGLETOL.");
    }
    void writeSolverReport(
            std::ostream& out, const Controls& c, const RunSettings& s,
            const ClosedOrbitSolver::Result& result) {
        out << "Integrator=" << ExternalFieldRayTracker::integrationMethodName(c.integrator)
            << " DT_s=" << c.tracking.dt << " p0_mc=" << c.tracking.momentum
            << " MAXSTEPS=" << c.tracking.maxSteps << " MAXPATH_m=" << c.tracking.maxPath
            << " T0_s=" << c.tracking.time << "\nXTOL_m=" << s.solver.positionTolerance
            << " PTOL_mc=" << s.solver.momentumTolerance << " MAXIT=" << s.solver.maxIterations
            << " DAMPING=" << s.solver.damping << "\nStatus=" << int(result.status) << ' '
            << result.message << " evaluations=" << result.evaluations << '\n';
        out << "FDSTEP=";
        for (double v : s.solver.finiteDifferenceSteps)
            out << v << ' ';
        out << "\nSCALES=";
        for (double v : s.solver.scales)
            out << v << ' ';
        out << '\n';
        for (const auto& iteration : result.iterations) {
            out << "iteration damping=" << iteration.damping
                << " pivot_ratio=" << iteration.pivotRatio << " residual=";
            for (double v : iteration.residual)
                out << v << ' ';
            out << " correction=";
            for (double v : iteration.correction)
                out << v << ' ';
            out << '\n';
        }
        out << "Coordinates (x,px,y,py): ";
        for (double v : result.coordinates)
            out << v << ' ';
        out << "\nResidual: ";
        for (double v : result.residual)
            out << v << ' ';
        out << '\n';
    }
    void setLaunchState(
            ClosedOrbitInitialState& state, const Controls& c,
            const ClosedOrbitSolver::Coordinates& u) {
        state.position = c.section.transformFrom(Vector_t<double, 3>(u[0], u[2], 0));
        state.momentum = c.section.rotateFrom(
                Vector_t<double, 3>(
                        u[1], u[3],
                        std::sqrt(
                                c.tracking.momentum * c.tracking.momentum - u[1] * u[1]
                                - u[3] * u[3])));
    }
    std::string serializeResult(
            const ClosedOrbitInitialState& state, const Controls& c, const RunSettings& s,
            const std::string& beamName, const ClosedOrbitSolver::Result& result,
            const PartData& reference, const LinearMapEigenAnalysis::Result& spectrum,
            const OneTurnMap::FullTunes& fullTunes, const OneTurnMap::Result& returned,
            double drift) {
        const std::string speciesId = lowerCase(state.species);
        const char* stability       = "MARGINAL";
        using Stability             = LinearMapEigenAnalysis::Stability;
        if (spectrum.stability == Stability::Stable) stability = "STABLE";
        if (spectrum.stability == Stability::Unstable) stability = "UNSTABLE";
        if (spectrum.stability == Stability::NonUnitCircle) stability = "NON_UNIT_CIRCLE";
        CofJsonWriter json;
        json.object(nullptr, [&] {
            json.text("schema", "opalx.closed_orbits");
            json.number("schema_version", 1);
            json.array("species", [&] {
                json.object(nullptr, [&] {
                    json.text("id", speciesId);
                    json.text("name", state.species);
                    json.number("rest_energy_eV", state.massEV);
                    json.number("charge_e", state.chargeE);
                });
            });
            json.array("orbits", [&] {
                json.object(nullptr, [&] {
                    json.text("id", "orbit_" + lowerCase(beamName));
                    json.text("beam_id", beamName);
                    json.text("cof_id", s.commandName);
                    json.text("line_id", state.line);
                    json.object("reference", [&] {
                        json.text("species_id", speciesId);
                        json.number(
                                "kinetic_energy_MeV", (reference.getE() - reference.getM()) / 1e6);
                    });
                    json.object("section", [&] {
                        json.numbers("origin_m", c.section.transformFrom(Vector_t<double, 3>(0.0)));
                        json.numbers("lab_to_section_quaternion_wxyz", c.section.getRotation());
                        json.array("coordinate_order", [&] {
                            for (const auto* name : {"x", "px", "y", "py"})
                                json.text(nullptr, name);
                        });
                        json.array("coordinate_units", [&] {
                            for (const auto* unit : {"m", "p/(mc)", "m", "p/(mc)"})
                                json.text(nullptr, unit);
                        });
                    });
                    json.object("solution", [&] {
                        json.numbers("section_coordinates", result.coordinates);
                        json.numbers("closure_residual", result.residual);
                        json.object("launch", [&] {
                            json.text("frame", "laboratory");
                            json.number("time_s", state.time);
                            json.numbers("position_m", state.position);
                            json.numbers("momentum_mc", state.momentum);
                        });
                    });
                    json.object("solver", [&] {
                        json.text("status", "converged");
                        json.text("method", "NEWTON");
                        json.text(
                                "integrator",
                                ExternalFieldRayTracker::integrationMethodName(c.integrator));
                        json.number("dt_s", c.tracking.dt);
                        json.object("settings", [&] {
                            json.number("position_tolerance_m", s.solver.positionTolerance);
                            json.number("momentum_tolerance_mc", s.solver.momentumTolerance);
                            json.number("max_iterations", s.solver.maxIterations);
                            json.number("max_steps", c.tracking.maxSteps);
                            json.number("max_path_m", c.tracking.maxPath);
                            json.numbers("initial_coordinates", s.initial);
                            json.numbers("finite_difference_steps", s.solver.finiteDifferenceSteps);
                            json.numbers("coordinate_scales", s.solver.scales);
                            json.boolean("damping", s.solver.damping);
                            json.number("rank_tolerance", s.solver.rankTolerance);
                            json.number("max_backtracks", s.solver.maxBacktracks);
                            json.number("geometry_tolerance_m", c.geometryTolerance);
                            json.number("angle_tolerance_rad", c.angleTolerance);
                        });
                        json.number("iterations", result.iterations.size());
                        json.number("map_evaluations", result.evaluations);
                        json.text("diagnostic", result.message);
                    });
                    json.object("optics", [&] {
                        json.array("one_turn_matrix", [&] {
                            for (const auto& row : result.matrix)
                                json.numbers(nullptr, row);
                        });
                        json.text("stability", stability);
                        json.boolean("near_integer", spectrum.nearInteger);
                        json.array("eigenvalues", [&] {
                            for (const auto& value : spectrum.eigenvalues)
                                json.object(nullptr, [&] {
                                    json.number("real", value.real());
                                    json.number("imag", value.imag());
                                });
                        });
                        json.number(
                                "maximum_relative_eigenpair_residual", spectrum.maximumResidual);
                        json.number(
                                "basis_reciprocal_condition", spectrum.basisReciprocalCondition);
                        json.array("fractional_modes", [&] {
                            for (const auto& mode : spectrum.modes)
                                json.object(nullptr, [&] {
                                    if (mode.tune)
                                        json.number("tune", *mode.tune);
                                    else
                                        json.null("tune");
                                    if (mode.complementaryTune)
                                        json.number("complement", *mode.complementaryTune);
                                    else
                                        json.null("complement");
                                });
                        });
                        json.object("full_tunes", [&] {
                            json.text("status", fullTunes.available ? "available" : "unavailable");
                            json.text("method", "uncoupled_planar_phase_winding");
                            if (fullTunes.available) {
                                json.number("horizontal_radial", fullTunes.tunes[0]);
                                json.number("vertical", fullTunes.tunes[1]);
                                json.number("phase_samples", fullTunes.samples);
                                json.number("maximum_phase_step_rad", fullTunes.maximumPhaseStep);
                            } else {
                                json.null("horizontal_radial");
                                json.null("vertical");
                                json.null("phase_samples");
                                json.null("maximum_phase_step_rad");
                            }
                            json.text("diagnostic", fullTunes.diagnostic);
                        });
                    });
                    json.object("diagnostics", [&] {
                        json.number("relative_energy_drift", drift);
                        json.number(
                                "return_duration_s",
                                returned.ray.time - state.time - returned.ray.timeCorrection);
                        json.number(
                                "path_length_m",
                                returned.ray.pathLength - returned.ray.pathLengthCorrection);
                        json.number("section_residual_m", returned.sectionResidual);
                    });
                });
            });
        });
        return json.str();
    }
    OneTurnMap makeReturnMap(
            OpalBeamline& lattice, ExternalFieldRayTracker& tracker, const Controls& c) {
        const auto elements = lattice.getElements();
        return OneTurnMap(
                [&lattice, &tracker, elements](
                        const auto& ray, double dt, std::vector<OneTurnMap::Step>* accepted) {
                    checkApertures(lattice, elements, ray);
                    std::vector<OneTurnMap::Step> local;
                    auto* intervals = accepted ? accepted : &local;
                    const auto end  = tracker.advance(ray, dt, intervals);
                    for (const auto& interval : *intervals) {
                        checkApertures(lattice, elements, interval.midpoint);
                        checkApertures(lattice, elements, interval.end);
                    }
                    return end;
                },
                c.section, c.tracking, tracker.maximumTimeStep());
    }
    Report solveAndReport(
            BeamSequence& sequence, Beam& beam, OpalBeamline& lattice,
            const std::vector<std::shared_ptr<ElementBase>>& ordered, const Controls& c,
            const RunSettings& s, ClosedOrbitInitialState& state) {
        require(!ordered.empty(), "RING must contain magnetic/transport elements.");
        if (!s.output.empty()) validateOutput(s.output);
        std::ostringstream out;
        validateClosure(lattice, ordered, c, sequence.getOpalName(), beam.getOpalName(), out);
        ExternalFieldRayTracker tracker(lattice, beam.getReference(), c.integrator);
        auto map          = makeReturnMap(lattice, tracker, c);
        const auto result = ClosedOrbitSolver::solve(map, s.initial, s.solver);
        writeSolverReport(out, c, s, result);
        if (result.status != ClosedOrbitSolver::Status::Converged) {
            throw OpalException("COF", out.str());
        }
        const auto returned   = map(result.coordinates);
        const auto& reference = beam.getReference();
        const double drift    = (std::sqrt(1 + dot(returned.ray.momentum, returned.ray.momentum))
                              - reference.getGamma())
                             / (reference.getGamma() - 1);
        out << "Return time_s=" << returned.ray.time - c.tracking.time
            << " path_m=" << returned.ray.pathLength << " relative_energy_drift=" << drift << '\n';
        LinearMapEigenAnalysis::Settings ev;
        ev.scales = s.solver.scales;
        const auto opticsMap =
                opalx::maps::toTpsClosedOrbitMap(result, c.section, c.tracking.momentum);
        const auto spectrum = LinearMapEigenAnalysis::analyze(opticsMap, ev);
        std::ostringstream optics;
        optics << std::setprecision(17) << "COF matrix and tunes: ring=" << sequence.getOpalName()
               << " beam=" << beam.getOpalName()
               << " Near-integer warning: " << (spectrum.nearInteger ? "yes" : "no")
               << "\nOne-turn-matrix (x, px, y, py):\n";
        const char* coordinates[] = {"x", "px", "y", "py"};
        std::ostringstream matrix;
        matrix << "     ";
        for (const auto* coordinate : coordinates)
            matrix << std::setw(16) << coordinate;
        matrix << '\n' << std::scientific << std::setprecision(8);
        for (unsigned i = 0; i < 4; ++i) {
            matrix << std::setw(3) << coordinates[i] << " [";
            for (double v : result.matrix[i])
                matrix << std::setw(16) << v;
            matrix << " ]\n";
        }
        optics << matrix.str();
        std::ostringstream eigenanalysis;
        LinearMapEigenAnalysis::writeReport(eigenanalysis, spectrum);
        std::istringstream lines(eigenanalysis.str());
        for (std::string line; std::getline(lines, line);) {
            if (line.rfind("Near-integer warning:", 0) == 0
                || line.rfind("Matrix eigenanalysis alone:", 0) == 0)
                continue;
            // Keep the shared report API unchanged; route its eigenvalue rows
            // to the detailed COF report rather than the level-4 optics summary.
            (line.rfind("lambda[", 0) == 0 ? out : optics) << line << '\n';
        }
        OneTurnMap::FullTunes fullTunes;
        if (spectrum.stability == LinearMapEigenAnalysis::Stability::Stable)
            fullTunes = map.fullTunes(
                    result.coordinates, s.solver.finiteDifferenceSteps, result.matrix);
        else
            fullTunes.diagnostic = "One-turn eigenanalysis is not numerically stable.";
        if (fullTunes.available)
            optics << "Full tunes (uncoupled planar):\n"
                   << "  " << std::left << std::setw(22) << "horizontal/radial"
                   << "nur=" << fullTunes.tunes[0] << '\n'
                   << "  " << std::setw(22) << "vertical" << "nuy=" << fullTunes.tunes[1] << '\n';
        else
            optics << "Full tunes unavailable: " << fullTunes.diagnostic << '\n';
        out << "Phase samples=" << fullTunes.samples
            << " maximum phase step_rad=" << fullTunes.maximumPhaseStep << '\n';
        setLaunchState(state, c, result.coordinates);
        state.json = serializeResult(
                state, c, s, beam.getOpalName(), result, reference, spectrum, fullTunes, returned,
                drift);
        if (!s.output.empty()) writeResult(s.output, state.json);
        return {out.str(), optics.str()};
    }

    // All ranks build the empty reference container once; scalar ray work and files
    // are rank-zero only. Broadcast diagnostic text before throwing on solver failure.
    ClosedOrbitInitialState calculate(
            BeamSequence& sequence, Beam& beam, const Controls& c, const RunSettings& s) {
        sequence.prepareForTracking();
        // COF needs only empty particle storage for element initialization, not a solver.
        const opalx::spacecharge::CartesianDomainConfig3D domain;
        PartBunch_t bunch({beam.getCharge()}, {beam.getMass()}, {&beam}, {0}, 1, "LF2", domain);
        OpalBeamline lattice(
                sequence.fetchLine()->getOrigin3D(), sequence.fetchLine()->getInitialDirection());
        LatticeVisitor visitor(*sequence.fetchLine(), lattice, bunch, sequence.getOpalName());
        visitor.execute();
        lattice.prepareSections();
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        ClosedOrbitInitialState state;
        state.section = c.section;
        state.time    = c.tracking.time;
        state.massEV  = beam.getReference().getM();
        state.chargeE = beam.getReference().getQ();
        state.line    = sequence.getOpalName();
        state.species = beam.getParticleName();
        std::string error;
        Report report;
        if (rank == 0) try {
                report = solveAndReport(sequence, beam, lattice, visitor.ordered, c, s, state);
            } catch (const OpalException& ex) {
                error = ex.what();
            } catch (const std::exception& ex) {
                error = ex.what();
            }
        auto broadcast = [&](std::string& value) {
            int length = static_cast<int>(value.size());
            MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);
            value.resize(length);
            MPI_Bcast(value.data(), length, MPI_CHAR, 0, MPI_COMM_WORLD);
        };
        broadcast(error);
        require(error.empty(), error);
        broadcast(state.json);
        double values[6]{};
        if (rank == 0)
            for (unsigned d = 0; d < 3; ++d) {
                values[d]     = state.position[d];
                values[d + 3] = state.momentum[d];
            }
        MPI_Bcast(values, 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (unsigned d = 0; d < 3; ++d) {
            state.position[d] = values[d];
            state.momentum[d] = values[d + 3];
        }
        if (rank == 0) {
            *gmsg << level5 << report.details << endl;
            *gmsg << level4 << report.optics << endl;
        }
        return state;
    }
}  // namespace

CofCmd::CofCmd()
    : Action(SIZE, "COF", "Find a fixed-energy closed orbit of a static magnetic RING.") {
    itsAttr[LINE] = Attributes::makeString("LINE", "Required RING name.");
    itsAttr[BEAM] =
            Attributes::makeString("BEAM", "Required BEAM with explicit PC, ENERGY or GAMMA.");
    itsAttr[DT]             = Attributes::makeReal("DT", "Positive ray timestep [s].", 1e-12);
    itsAttr[MAXSTEPS]       = Attributes::makeReal("MAXSTEPS", "One-return step budget.", 200000);
    itsAttr[T0]             = Attributes::makeReal("T0", "Initial time [s].", 0);
    itsAttr[TIMEINTEGRATOR] = Attributes::makePredefinedString(
            "TIMEINTEGRATOR", "Shared ray integrator.", {"BORIS", "LF2", "RK4", "DOP853"}, "RK4");
    itsAttr[MAXPATH] =
            Attributes::makeReal("MAXPATH", "Required positive return search bound [m].");
    itsAttr[SECTION] = Attributes::makeRealArray(
            "SECTION", "Optional global X,Y,Z,THETA,PHI,PSI [m,rad]; default RING frame.");
    itsAttr[GEOMTOL] =
            Attributes::makeReal("GEOMTOL", "Nominal end-to-start position tolerance [m].", 1e-9);
    itsAttr[ANGLETOL] =
            Attributes::makeReal("ANGLETOL", "Nominal end-to-start axis tolerance [rad].", 1e-10);
    itsAttr[METHOD] = Attributes::makePredefinedString("METHOD", "Solver.", {"NEWTON"}, "NEWTON");
    itsAttr[DIMENSION]  = Attributes::makeReal("DIMENSION", "Only 4 is supported.", 4);
    const char* names[] = {"X", "PX", "Y", "PY"};
    for (unsigned i = 0; i < 4; ++i)
        itsAttr[X + i] = Attributes::makeReal(names[i], "Section coordinate [m or p/(mc)].", 0);
    itsAttr[MAXIT] = Attributes::makeReal("MAXIT", "Maximum Newton iterations.", 20);
    itsAttr[XTOL] =
            Attributes::makeReal("XTOL", "Position residual and correction tolerance [m].", 1e-10);
    itsAttr[PTOL] = Attributes::makeReal(
            "PTOL", "Momentum residual and correction tolerance [p/(mc)].", 1e-10);
    itsAttr[FDSTEP] = Attributes::makeRealArray(
            "FDSTEP", "Four absolute central-difference steps; default all 1e-6.");
    itsAttr[SCALES] = Attributes::makeRealArray("SCALES", "Four coordinate scales; default all 1.");
    itsAttr[DAMPING]  = Attributes::makeBool("DAMPING", "Use Newton backtracking.", true);
    itsAttr[JACOBIAN] = Attributes::makePredefinedString(
            "JACOBIAN", "Differentiation.", {"CENTRAL"}, "CENTRAL");
    itsAttr[OUTPUT] = Attributes::makeString(
            "OUTPUT", "Optional exact JSON filename (replaced after a successful solve).");
    registerOwnership(AttributeHandler::COMMAND);
}
CofCmd::CofCmd(const std::string& name, CofCmd* parent) : Action(name, parent) {}
CofCmd* CofCmd::clone(const std::string& name) { return new CofCmd(name, this); }
void CofCmd::execute() {
    result_m.reset();
    OpalData::getInstance()->hasCofRun = true;
    require(!OpalData::getInstance()->inRestartRun(), "COF does not support restart mode.");
    auto* sequence = BeamSequence::find(Attributes::getString(itsAttr[LINE]));
    require(sequence->fetchLine()->getBeamlineTopology() == BeamlineTopology::RING,
            "LINE must name a RING.");
    auto* beam = Beam::find(Attributes::getString(itsAttr[BEAM]));
    require(beam->hasExplicitEnergy(), "BEAM requires explicit PC, ENERGY (total GeV) or GAMMA.");
    require(beam->getGlobalProcessNames().empty(), "COF does not support BEAM global processes.");
    const auto& reference = beam->getReference();
    positive(reference.getM(), "Rest mass");
    require(std::isfinite(reference.getQ()) && reference.getQ() != 0,
            "COF requires nonzero finite charge.");
    Controls c;
    c.tracking.momentum = reference.getP() / reference.getM();
    positive(c.tracking.momentum, "Reference momentum");
    c.tracking.dt = Attributes::getReal(itsAttr[DT]);
    positive(c.tracking.dt, "DT");
    c.tracking.maxSteps = count(Attributes::getReal(itsAttr[MAXSTEPS]), "MAXSTEPS");
    require(!itsAttr[MAXPATH].defaultUsed(), "MAXPATH must specify the return search bound [m].");
    c.tracking.maxPath = Attributes::getReal(itsAttr[MAXPATH]);
    positive(c.tracking.maxPath, "MAXPATH");
    c.tracking.time = Attributes::getReal(itsAttr[T0]);
    require(std::isfinite(c.tracking.time), "T0 must be finite.");
    c.geometryTolerance = Attributes::getReal(itsAttr[GEOMTOL]);
    positive(c.geometryTolerance, "GEOMTOL");
    c.angleTolerance = Attributes::getReal(itsAttr[ANGLETOL]);
    positive(c.angleTolerance, "ANGLETOL");
    c.integrator = ExternalFieldRayTracker::parseIntegrationMethod(
            Attributes::getString(itsAttr[TIMEINTEGRATOR]));
    c.section = CoordinateSystemTrafo(
            sequence->fetchLine()->getOrigin3D(), sequence->fetchLine()->getInitialDirection());
    if (!itsAttr[SECTION].defaultUsed()) {
        const auto v = Attributes::getRealArray(itsAttr[SECTION]);
        require(v.size() == 6, "SECTION requires X,Y,Z,THETA,PHI,PSI.");
        for (double x : v)
            require(std::isfinite(x), "SECTION must be finite.");
        Quaternion theta(std::cos(v[3] / 2), 0, std::sin(v[3] / 2), 0);
        Quaternion phi(std::cos(v[4] / 2), std::sin(v[4] / 2), 0, 0);
        Quaternion psi(std::cos(v[5] / 2), 0, 0, std::sin(v[5] / 2));
        c.section = CoordinateSystemTrafo(
                Vector_t<double, 3>(v[0], v[1], v[2]), (theta * (phi * psi)).conjugate());
    }
    require(Attributes::getReal(itsAttr[DIMENSION]) == 4, "Only DIMENSION=4 is supported.");
    RunSettings s;
    s.commandName = getOpalName();
    for (unsigned i = 0; i < 4; ++i)
        s.initial[i] = Attributes::getReal(itsAttr[X + i]);
    s.solver.maxIterations     = count(Attributes::getReal(itsAttr[MAXIT]), "MAXIT");
    s.solver.positionTolerance = Attributes::getReal(itsAttr[XTOL]);
    s.solver.momentumTolerance = Attributes::getReal(itsAttr[PTOL]);
    s.solver.damping           = Attributes::getBool(itsAttr[DAMPING]);
    const std::array arrayAttributes{
            std::make_pair(FDSTEP, &s.solver.finiteDifferenceSteps),
            std::make_pair(SCALES, &s.solver.scales)};
    for (const auto& item : arrayAttributes) {
        if (!itsAttr[item.first].defaultUsed()) {
            const auto values = Attributes::getRealArray(itsAttr[item.first]);
            require(values.size() == 4, "FDSTEP and SCALES require exactly four values.");
            std::copy(values.begin(), values.end(), item.second->begin());
        }
    }
    if (!itsAttr[OUTPUT].defaultUsed()) {
        s.output = Attributes::getString(itsAttr[OUTPUT]);
        require(!s.output.empty(), "Explicit OUTPUT must not be empty.");
    }
    result_m = calculate(*sequence, *beam, c, s);
}

const ClosedOrbitInitialState& CofCmd::findResult(const std::string& name) {
    auto* command = dynamic_cast<CofCmd*>(OpalData::getInstance()->find(name));
    require(command && command->result_m.has_value(),
            "INITIALORBIT must name a successfully executed COF: " + name);
    return *command->result_m;
}
