// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Track/CofCmd.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/ClosedOrbitSolver.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/LinearMapEigenAnalysis.h"
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
                if (!s.output.empty())
                    require(!std::filesystem::exists(s.output), "COF output exists: " + s.output);
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
                std::ostringstream json;
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
                const auto& u  = result.coordinates;
                state.position = c.section.transformFrom(Vector_t<double, 3>(u[0], u[2], 0));
                state.momentum = c.section.rotateFrom(
                        Vector_t<double, 3>(
                                u[1], u[3],
                                std::sqrt(
                                        c.tracking.momentum * c.tracking.momentum - u[1] * u[1]
                                        - u[3] * u[3])));
                json << "],\"relative_energy_drift\":" << drift << ",\"time_s\":" << state.time
                     << ",\"position_m\":[";
                array(state.position);
                json << "],\"momentum_mc\":[";
                array(state.momentum);
                json << "],\"line\":" << std::quoted(state.line)
                     << ",\"species\":" << std::quoted(state.species)
                     << ",\"mass_eV\":" << state.massEV << ",\"charge_e\":" << state.chargeE
                     << ",\"section_origin_m\":[";
                array(c.section.transformFrom(Vector_t<double, 3>(0.0)));
                json << "],\"section_rotation_wxyz\":[";
                array(c.section.getRotation());
                json << "]}\n";
                state.json = json.str();
                if (!s.output.empty()) {
                    std::ofstream file(s.output);
                    file << state.json;
                    file.flush();
                    require(bool(file), "Error writing COF OUTPUT: " + s.output);
                }
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
        if (rank == 0) *gmsg << report << endl;
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
            "OUTPUT", "Optional exact JSON filename (existing files rejected).");
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
    for (unsigned i = 0; i < 4; ++i)
        s.initial[i] = Attributes::getReal(itsAttr[X + i]);
    s.solver.maxIterations     = count(Attributes::getReal(itsAttr[MAXIT]), "MAXIT");
    s.solver.positionTolerance = Attributes::getReal(itsAttr[XTOL]);
    s.solver.momentumTolerance = Attributes::getReal(itsAttr[PTOL]);
    s.solver.damping           = Attributes::getBool(itsAttr[DAMPING]);
    const std::array arrayAttributes{
            std::make_pair(FDSTEP, &s.solver.finiteDifferenceSteps),
            std::make_pair(SCALES, &s.solver.scales)};
    for (const auto item : arrayAttributes) {
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
