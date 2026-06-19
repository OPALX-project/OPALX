#ifndef OPAL_POISSON_BACKENDS_H
#define OPAL_POISSON_BACKENDS_H

#include <variant>

#include "Manager/FieldSolverBase.h"
#include "PartBunch/Solvers/PoissonBackend.hpp"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

template <typename T, unsigned Dim, typename Solver>
class IPPLPoissonBackendBase : public PoissonBackend<T, Dim> {
protected:
    Solver_t<T, Dim>& solverVariant_m;
    Field_t<Dim>* rho_m;
    VField_t<T, Dim>* E_m;
    Field_t<Dim>* phi_m;

    IPPLPoissonBackendBase(
            Solver_t<T, Dim>& solverVariant, Field_t<Dim>* rho, VField_t<T, Dim>* E,
            Field_t<Dim>* phi)
        : solverVariant_m(solverVariant), rho_m(rho), E_m(E), phi_m(phi) {}

    Solver& solver() { return std::get<Solver>(solverVariant_m); }
    const Solver& solver() const { return std::get<Solver>(solverVariant_m); }

    void emplaceSolver(const ippl::ParameterList& sp, const std::string& context) {
        solverVariant_m.template emplace<Solver>();
        Solver& activeSolver = solver();
        activeSolver.mergeParameters(sp);

        if (!rho_m) {
            throw OpalException(context, "rho_m is not initialized.");
        }
        if (!E_m) {
            throw OpalException(context, "E_m is not initialized.");
        }

        activeSolver.setRhs(*rho_m);
        activeSolver.setLhs(*E_m);
    }

    void rejectShiftedGreens(const SolveRequest<T, Dim>& request, const std::string& context) const {
        if (request.hasShiftedGreens()) {
            throw OpalException(
                    context, "Shifted Greens solves are not supported by " + this->name());
        }
    }

public:
    void refreshAfterLayoutChange() override { solver().setRhs(*rho_m); }
};

template <typename T, unsigned Dim>
class NullPoissonBackend final
    : public IPPLPoissonBackendBase<T, Dim, NullSolver_t<T, Dim>> {
    using Solver = NullSolver_t<T, Dim>;
    using Base   = IPPLPoissonBackendBase<T, Dim, Solver>;

public:
    NullPoissonBackend(
            Solver_t<T, Dim>& solverVariant, Field_t<Dim>* rho, VField_t<T, Dim>* E,
            Field_t<Dim>* phi)
        : Base(solverVariant, rho, E, phi) {
        ippl::ParameterList sp;
        this->emplaceSolver(sp, "NullPoissonBackend");
    }

    std::string name() const override { return "NONE"; }

    SolverCapabilities capabilities() const override { return {.isNoOp = true}; }

    void solve(const SolveRequest<T, Dim>& request) override {
        this->rejectShiftedGreens(request, "NullPoissonBackend::solve");
        this->solver().solve();
    }

    T couplingConstant() const override { return 1.0 / (4.0 * Physics::pi * Physics::epsilon_0); }
};

template <typename T, unsigned Dim>
class FFTPoissonBackend final : public IPPLPoissonBackendBase<T, Dim, FFTSolver_t<T, Dim>> {
    using Solver = FFTSolver_t<T, Dim>;
    using Base   = IPPLPoissonBackendBase<T, Dim, Solver>;

public:
    FFTPoissonBackend(
            Solver_t<T, Dim>& solverVariant, Field_t<Dim>* rho, VField_t<T, Dim>* E,
            Field_t<Dim>* phi)
        : Base(solverVariant, rho, E, phi) {
        ippl::ParameterList sp;
        sp.add("output_type", Solver::GRAD);
        sp.add("use_heffte_defaults", false);
        sp.add("use_pencils", true);
        sp.add("use_reorder", false);
        sp.add("use_gpu_aware", true);
        sp.add("comm", ippl::p2p_pl);
        sp.add("r2c_direction", 0);
        this->emplaceSolver(sp, "FFTPoissonBackend");
    }

    std::string name() const override { return "FFT"; }

    SolverCapabilities capabilities() const override {
        return {
                .supportsShiftedGreens   = false,
                .providesPotential       = false,
                .providesGradient        = true,
                .requiresPotentialBCs    = false,
                .debugDumpRhoBeforeSolve = true,
                .debugDumpScalarAfterSolve = true,
                .debugDumpVectorAfterSolve = true};
    }

    void solve(const SolveRequest<T, Dim>& request) override {
        this->rejectShiftedGreens(request, "FFTPoissonBackend::solve");
        this->solver().solve();
    }

    T couplingConstant() const override { return 1.0 / Physics::epsilon_0; }
};

template <typename T, unsigned Dim>
class OpenPoissonBackend final : public IPPLPoissonBackendBase<T, Dim, OpenSolver_t<T, Dim>> {
    using Solver = OpenSolver_t<T, Dim>;
    using Base   = IPPLPoissonBackendBase<T, Dim, Solver>;

public:
    OpenPoissonBackend(
            Solver_t<T, Dim>& solverVariant, Field_t<Dim>* rho, VField_t<T, Dim>* E,
            Field_t<Dim>* phi)
        : Base(solverVariant, rho, E, phi) {
        ippl::ParameterList sp;
        sp.add("output_type", Solver::SOL_AND_GRAD);
        sp.add("use_heffte_defaults", false);
        sp.add("use_pencils", true);
        sp.add("use_reorder", false);
        sp.add("use_gpu_aware", true);
        sp.add("comm", ippl::p2p_pl);
        sp.add("r2c_direction", 0);
        sp.add("algorithm", Solver::HOCKNEY);
        this->emplaceSolver(sp, "OpenPoissonBackend");
    }

    std::string name() const override { return "OPEN"; }

    SolverCapabilities capabilities() const override {
        return {
                .supportsShiftedGreens   = true,
                .providesPotential       = true,
                .providesGradient        = true,
                .requiresPotentialBCs    = false,
                .debugDumpRhoBeforeSolve = true,
                .debugDumpScalarAfterSolve = true,
                .debugDumpVectorAfterSolve = true};
    }

    void solve(const SolveRequest<T, Dim>& request) override {
        if (request.greensShift.has_value()) {
            this->solver().shiftedGreensFunction(*request.greensShift);
            this->solver().solve();
            this->solver().greensFunction();
            return;
        }

        this->solver().solve();
    }

    T couplingConstant() const override { return 1.0 / Physics::epsilon_0; }
};

template <typename T, unsigned Dim>
class CGPoissonBackend final : public IPPLPoissonBackendBase<T, Dim, CGSolver_t<T, Dim>> {
    using Solver = CGSolver_t<T, Dim>;
    using Base   = IPPLPoissonBackendBase<T, Dim, Solver>;

public:
    CGPoissonBackend(
            Solver_t<T, Dim>& solverVariant, Field_t<Dim>* rho, VField_t<T, Dim>* E,
            Field_t<Dim>* phi)
        : Base(solverVariant, rho, E, phi) {
        ippl::ParameterList sp;
        sp.add("output_type", Solver::GRAD);
        sp.add("tolerance", 1e-12);

        solverVariant.template emplace<Solver>();
        Solver& activeSolver = this->solver();
        activeSolver.mergeParameters(sp);
        if (!rho) {
            throw OpalException("CGPoissonBackend", "rho_m is not initialized.");
        }
        activeSolver.setRhs(*rho);

        throw OpalException("CGPoissonBackend", "Cannot use CGSolver yet, not fully implemented.");
    }

    std::string name() const override { return "CG"; }

    SolverCapabilities capabilities() const override {
        return {
                .supportsShiftedGreens   = false,
                .providesPotential       = true,
                .providesGradient        = true,
                .requiresPotentialBCs    = true,
                .usesSeparatePotentialField = true,
                .debugDumpRhoBeforeSolve = true,
                .debugDumpScalarAfterSolve = true,
                .debugDumpVectorAfterSolve = false};
    }

    void solve(const SolveRequest<T, Dim>& request) override {
        this->rejectShiftedGreens(request, "CGPoissonBackend::solve");
        this->solver().solve();
    }

    T couplingConstant() const override { return 1.0 / Physics::epsilon_0; }
};

#endif  // OPAL_POISSON_BACKENDS_H
