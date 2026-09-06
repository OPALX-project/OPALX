// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Track/CofCmd.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Directory.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/ClosedOrbitSolver.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/LinearMapEigenAnalysis.h"
#include "Attributes/Attributes.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "Elements/OpalBeamline.h"
#include "OpalParser/OpalParser.h"
#include "OpalParser/Statement.h"
#include "Structure/Beam.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolverCmd.h"
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
        SIZE
    };
    void require(bool ok, const std::string& message) {
        if (!ok) throw OpalException("COF", message);
    }
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
    };

    // The local parser owns these actions; clones retain callbacks only for the
    // duration of this block. No static command context or globally stored action.
    class CofRun : public Action {
        enum {
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
            N
        };
        std::function<void(const RunSettings&)> run;
        CofRun(const std::string& name, CofRun* parent) : Action(name, parent), run(parent->run) {}

    public:
        explicit CofRun(std::function<void(const RunSettings&)> callback)
            : Action(N, "RUN", "Solve the static fixed-energy 4D closed orbit."),
              run(std::move(callback)) {
            itsAttr[METHOD] =
                    Attributes::makePredefinedString("METHOD", "Solver.", {"NEWTON"}, "NEWTON");
            itsAttr[DIMENSION]  = Attributes::makeReal("DIMENSION", "Only 4 is supported.", 4);
            const char* names[] = {"X", "PX", "Y", "PY"};
            for (unsigned i = 0; i < 4; ++i)
                itsAttr[X + i] =
                        Attributes::makeReal(names[i], "Section coordinate [m or p/(mc)].", 0);
            itsAttr[MAXIT] = Attributes::makeReal("MAXIT", "Maximum Newton iterations.", 20);
            itsAttr[XTOL]  = Attributes::makeReal(
                    "XTOL", "Position residual and correction tolerance [m].", 1e-10);
            itsAttr[PTOL] = Attributes::makeReal(
                    "PTOL", "Momentum residual and correction tolerance [p/(mc)].", 1e-10);
            itsAttr[FDSTEP] = Attributes::makeRealArray(
                    "FDSTEP", "Four absolute central-difference steps; default all 1e-6.");
            itsAttr[SCALES] =
                    Attributes::makeRealArray("SCALES", "Four coordinate scales; default all 1.");
            itsAttr[DAMPING]  = Attributes::makeBool("DAMPING", "Use Newton backtracking.", true);
            itsAttr[JACOBIAN] = Attributes::makePredefinedString(
                    "JACOBIAN", "Differentiation.", {"CENTRAL"}, "CENTRAL");
            itsAttr[OUTPUT] =
                    Attributes::makeString("OUTPUT", "Output prefix (existing files rejected).");
        }
        CofRun* clone(const std::string& name) override { return new CofRun(name, this); }
        void execute() override {
            require(Attributes::getReal(itsAttr[DIMENSION]) == 4, "Only DIMENSION=4 is supported.");
            RunSettings s;
            for (unsigned i = 0; i < 4; ++i)
                s.initial[i] = Attributes::getReal(itsAttr[X + i]);
            s.solver.maxIterations     = count(Attributes::getReal(itsAttr[MAXIT]), "MAXIT");
            s.solver.positionTolerance = Attributes::getReal(itsAttr[XTOL]);
            s.solver.momentumTolerance = Attributes::getReal(itsAttr[PTOL]);
            s.solver.damping           = Attributes::getBool(itsAttr[DAMPING]);
            for (auto item :
                 {std::make_pair(FDSTEP, &s.solver.finiteDifferenceSteps),
                  std::make_pair(SCALES, &s.solver.scales)}) {
                if (!itsAttr[item.first].defaultUsed()) {
                    const auto values = Attributes::getRealArray(itsAttr[item.first]);
                    require(values.size() == 4, "FDSTEP and SCALES require exactly four values.");
                    std::copy(values.begin(), values.end(), item.second->begin());
                }
            }
            s.output = itsAttr[OUTPUT].defaultUsed()
                               ? OpalData::getInstance()->getInputBasename() + "_cof"
                               : Attributes::getString(itsAttr[OUTPUT]);
            require(!s.output.empty(), "OUTPUT must not be empty.");
            run(s);
        }
    };
    class CofEnd : public Action {
        std::function<void()> end;
        CofEnd(const std::string& name, CofEnd* parent) : Action(name, parent), end(parent->end) {}

    public:
        explicit CofEnd(std::function<void()> callback)
            : Action(0, "ENDCOF", "Close the COF block."), end(std::move(callback)) {}
        CofEnd* clone(const std::string& name) override { return new CofEnd(name, this); }
        void execute() override { end(); }
    };
    class CofParser : public OpalParser {
        Directory directory;

    public:
        bool ended = false, ran = false;
        explicit CofParser(std::function<void(const RunSettings&)> callback) {
            directory.insert("RUN", new CofRun([this, callback](const RunSettings& s) {
                                 require(!ran, "Use one RUN per COF block.");
                                 callback(s);
                                 ran = true;
                             }));
            directory.insert("ENDCOF", new CofEnd([this]() {
                                 ended = true;
                                 stop();
                             }));
        }
        Object* find(const std::string& name) const override { return directory.find(name); }
        void parse(Statement& stat) const override {
            stat.start();
            require(stat.keyword("RUN") || stat.keyword("ENDCOF"),
                    "Only RUN and ENDCOF are allowed inside COF.");
            stat.start();
            // Bypass definition/assignment dispatch: context-owning actions cannot escape.
            parseAction(stat);
        }
    };

    class EmptyFieldSolver : public FieldSolverCmd {
    public:
        EmptyFieldSolver() {
            Attributes::setPredefinedString(itsAttr[FIELDSOLVER::TYPE], "NONE");
            for (auto index : {FIELDSOLVER::PARFFTX, FIELDSOLVER::PARFFTY, FIELDSOLVER::PARFFTZ})
                Attributes::setBool(itsAttr[index], true);
            setNX(8);
            setNY(8);
            setNZ(8);
            setFieldSolverCmdType();
            setDomainDecomposition();
        }
    };

    class LatticeVisitor : public DefaultVisitor {
        OpalBeamline& lattice;
        PartBunch_t& bunch;
        std::string owner;
        void applyDefault(const ElementBase& element) override {
            if (!element.getName().empty() && element.getName().front() == '#') return;
            const auto type = element.getType();
            // A monitor is passive for scalar optics. Do not initialise its
            // diagnostic sink (which could overwrite a TRACK monitor file).
            if (type == ElementType::MONITOR) return;
            require(type == ElementType::DRIFT || type == ElementType::MARKER
                            || type == ElementType::MULTIPOLE || type == ElementType::SBEND
                            || type == ElementType::RBEND || type == ElementType::CYCLOTRONSECTOR,
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

    // All ranks build the empty reference container once; scalar ray work and files
    // are rank-zero only. Broadcast diagnostic text before throwing on solver failure.
    void calculate(BeamSequence& sequence, Beam& beam, const Controls& c, const RunSettings& s) {
        sequence.prepareForTracking();
        EmptyFieldSolver fields;
        DataSink sink(std::vector<H5PartWrapper*>{}, false, 0, s.output);
        PartBunch_t bunch(
                {beam.getCharge()}, {beam.getMass()}, {&beam}, {0}, 1, "LF2", &fields, &sink);
        OpalBeamline lattice(
                sequence.fetchLine()->getOrigin3D(), sequence.fetchLine()->getInitialDirection());
        LatticeVisitor visitor(*sequence.fetchLine(), lattice, bunch, sequence.getOpalName());
        visitor.execute();
        lattice.prepareSections();
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::string error, report;
        if (rank == 0) try {
                require(!visitor.ordered.empty(), "RING must contain magnetic/transport elements.");
                const auto first = lattice.getNominalEntryTransform(visitor.ordered.front());
                const auto last  = lattice.getNominalExitTransform(visitor.ordered.back());
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
                    angleError =
                            std::max(angleError, std::atan2(euclidean_norm(axisCross), dot(a, b)));
                }
                std::ostringstream out;
                out << std::setprecision(17) << "COF ring=" << sequence.getOpalName()
                    << " beam=" << beam.getOpalName()
                    << "\nNominal closure position_m=" << positionError
                    << " angle_rad=" << angleError << '\n';
                require(positionError <= c.geometryTolerance && angleError <= c.angleTolerance,
                        out.str()
                                + "Nominal ring geometry does not close within GEOMTOL/ANGLETOL.");
                for (const auto& suffix : {".txt", ".json", "-orbit.csv"})
                    require(!std::filesystem::exists(s.output + suffix),
                            "COF output exists: " + s.output + suffix);
                ExternalFieldRayTracker tracker(lattice, beam.getReference(), c.integrator);
                // Field selection excludes points outside apertures. Select the
                // nearest finite conventional body independently of aperture,
                // then ask its native loss predicate. Checking every body would
                // falsely lose rays in the longitudinal slab of the opposite arm.
                // Stateless selection is also valid for bisection and shadow rays.
                const auto elements       = lattice.getElements();
                const auto checkApertures = [&](const ExternalFieldRayTracker::State& ray) {
                    double nearest = std::numeric_limits<double>::infinity();
                    for (const auto& element : elements) {
                        if (element->getType() == ElementType::CYCLOTRONSECTOR
                            || element->getGeometry().getElementLength() <= 0)
                            continue;
                        nearest = std::min(
                                nearest, bodyDistance2(
                                                 *element, lattice.transformToLocalCS(
                                                                   element, ray.position)));
                    }
                    for (const auto& element : elements) {
                        if (element->getType() != ElementType::CYCLOTRONSECTOR) {
                            if (element->getGeometry().getElementLength() <= 0) continue;
                            const double distance = bodyDistance2(
                                    *element, lattice.transformToLocalCS(element, ray.position));
                            // Floating-point tie allowance in squared metres;
                            // not an aperture or Newton convergence tolerance.
                            if (distance > nearest
                                                   + 64 * std::numeric_limits<double>::epsilon()
                                                             * std::max(1.0, nearest))
                                continue;
                        }
                        Vector_t<double, 3> electric(0.0), magnetic(0.0);
                        require(!element->applyToReferenceParticle(
                                        lattice.transformToLocalCS(element, ray.position),
                                        lattice.rotateToLocalCS(element, ray.momentum), ray.time,
                                        electric, magnetic),
                                "Aperture/material loss in " + element->getName());
                    }
                };
                OneTurnMap map(
                        [&](const auto& ray, double dt, std::vector<OneTurnMap::Step>* accepted) {
                            checkApertures(ray);
                            std::vector<OneTurnMap::Step> local;
                            auto* intervals = accepted ? accepted : &local;
                            const auto end  = tracker.advance(ray, dt, intervals);
                            for (const auto& interval : *intervals) {
                                checkApertures(interval.midpoint);
                                checkApertures(interval.end);
                            }
                            return end;
                        },
                        c.section, c.tracking);
                const auto result = ClosedOrbitSolver::solve(map, s.initial, s.solver);
                out << "Integrator=" << ExternalFieldRayTracker::integrationMethodName(c.integrator)
                    << " DT_s=" << c.tracking.dt << " p0_mc=" << c.tracking.momentum
                    << " MAXSTEPS=" << c.tracking.maxSteps << " MAXPATH_m=" << c.tracking.maxPath
                    << " T0_s=" << c.tracking.time << "\nXTOL_m=" << s.solver.positionTolerance
                    << " PTOL_mc=" << s.solver.momentumTolerance
                    << " MAXIT=" << s.solver.maxIterations << " DAMPING=" << s.solver.damping
                    << "\nStatus=" << int(result.status) << ' ' << result.message
                    << " evaluations=" << result.evaluations << '\n';
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
                out << "\nMatrix:\n";
                for (const auto& row : result.matrix) {
                    for (double v : row)
                        out << v << ' ';
                    out << '\n';
                }
                if (result.status != ClosedOrbitSolver::Status::Converged) {
                    std::ofstream file(s.output + ".txt");
                    file << out.str();
                    throw OpalException("COF", out.str());
                }
                const auto returned   = map(result.coordinates);
                const auto& reference = beam.getReference();
                const double drift =
                        (std::sqrt(1 + dot(returned.ray.momentum, returned.ray.momentum))
                         - reference.getGamma())
                        / (reference.getGamma() - 1);
                out << "Return time_s=" << returned.ray.time - c.tracking.time
                    << " path_m=" << returned.ray.pathLength << " relative_energy_drift=" << drift
                    << '\n';
                LinearMapEigenAnalysis::Settings ev;
                ev.scales           = s.solver.scales;
                const auto spectrum = LinearMapEigenAnalysis::analyze(result.matrix, ev);
                LinearMapEigenAnalysis::writeReport(out, spectrum);
                out << "Stability is a numerical matrix diagnostic; validate DT and FDSTEP "
                       "convergence.\n";
                std::ofstream json(s.output + ".json"), orbit(s.output + "-orbit.csv"),
                        file(s.output + ".txt");
                require(bool(json) && bool(orbit) && bool(file), "Cannot open COF output files.");
                json << std::setprecision(17) << "{\"converged\":true,\"dt_s\":" << c.tracking.dt
                     << ",\"energy_MeV\":" << (reference.getE() - reference.getM()) / 1e6
                     << ",\"coordinates\":[";
                auto array = [&](const auto& values) {
                    bool comma = false;
                    for (double v : values) {
                        if (comma) json << ',';
                        json << v;
                        comma = true;
                    }
                };
                array(result.coordinates);
                json << "],\"residual\":[";
                array(result.residual);
                json << "],\"matrix\":[";
                for (unsigned i = 0; i < 4; ++i) {
                    if (i) json << ',';
                    json << '[';
                    array(result.matrix[i]);
                    json << ']';
                }
                json << "],\"eigenvalues\":[";
                for (unsigned i = 0; i < 4; ++i) {
                    if (i) json << ',';
                    json << '[' << spectrum.eigenvalues[i].real() << ','
                         << spectrum.eigenvalues[i].imag() << ']';
                }
                const char* stability = "MARGINAL";
                using Stability       = LinearMapEigenAnalysis::Stability;
                if (spectrum.stability == Stability::Stable) stability = "STABLE";
                if (spectrum.stability == Stability::Unstable) stability = "UNSTABLE";
                if (spectrum.stability == Stability::NonUnitCircle) stability = "NON_UNIT_CIRCLE";
                json << "],\"stability\":\"" << stability
                     << "\",\"near_integer\":" << (spectrum.nearInteger ? "true" : "false")
                     << ",\"fractional_modes\":[";
                for (size_t i = 0; i < spectrum.modes.size(); ++i) {
                    if (i) json << ',';
                    const auto& mode = spectrum.modes[i];
                    json << "{\"tune\":";
                    if (mode.tune)
                        json << *mode.tune;
                    else
                        json << "null";
                    json << ",\"complement\":";
                    if (mode.complementaryTune)
                        json << *mode.complementaryTune;
                    else
                        json << "null";
                    json << '}';
                }
                json << "],\"relative_energy_drift\":" << drift << "}\n";
                orbit << std::setprecision(17) << "s_m,time_s,X_m,Y_m,Z_m,PX_mc,PY_mc,PZ_mc\n";
                ExternalFieldRayTracker::State ray;
                const auto& u = result.coordinates;
                ray.position  = c.section.transformFrom(Vector_t<double, 3>(u[0], u[2], 0));
                ray.momentum  = c.section.rotateFrom(
                        Vector_t<double, 3>(
                                u[1], u[3],
                                std::sqrt(
                                        c.tracking.momentum * c.tracking.momentum - u[1] * u[1]
                                        - u[3] * u[3])));
                ray.time   = c.tracking.time;
                auto write = [&](const auto& r) {
                    orbit << r.pathLength << ',' << r.time;
                    for (unsigned d = 0; d < 3; ++d)
                        orbit << ',' << r.position[d];
                    for (unsigned d = 0; d < 3; ++d)
                        orbit << ',' << r.momentum[d];
                    orbit << '\n';
                };
                write(ray);
                for (unsigned i = 1; i < returned.steps; ++i) {
                    ray = tracker.advance(ray, c.tracking.dt);
                    write(ray);
                }
                write(returned.ray);
                file << out.str();
                file.flush();
                json.flush();
                orbit.flush();
                require(bool(file) && bool(json) && bool(orbit), "Error writing COF output.");
                report = out.str();
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
        broadcast(report);
        require(error.empty(), error);
        if (rank == 0) *gmsg << report << endl;
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
    registerOwnership(AttributeHandler::COMMAND);
    AttributeHandler::addAttributeOwner("COF", AttributeHandler::COMMAND, "ENDCOF");
}
CofCmd::CofCmd(const std::string& name, CofCmd* parent) : Action(name, parent) {}
CofCmd* CofCmd::clone(const std::string& name) { return new CofCmd(name, this); }
void CofCmd::execute() {
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
    CofParser parser([&](const RunSettings& s) {
        calculate(*sequence, *beam, c, s);
    });
    parser.run();
    require(parser.ended && parser.ran, "COF requires one RUN followed by ENDCOF (not EOF/STOP).");
}
