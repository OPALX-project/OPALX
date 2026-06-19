#include "PartBunch/Solvers/PoissonBackendRegistry.hpp"

#include <algorithm>

#include "PartBunch/Solvers/PoissonBackends.hpp"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

namespace {
    SolverCapabilities noneCapabilities() { return {.isNoOp = true}; }

    SolverCapabilities fftCapabilities() {
        return {
                .supportsShiftedGreens   = false,
                .providesPotential       = false,
                .providesGradient        = true,
                .requiresPotentialBCs    = false,
                .debugDumpRhoBeforeSolve = true,
                .debugDumpScalarAfterSolve = true,
                .debugDumpVectorAfterSolve = true};
    }

    SolverCapabilities openCapabilities() {
        return {
                .supportsShiftedGreens   = true,
                .providesPotential       = true,
                .providesGradient        = true,
                .requiresPotentialBCs    = false,
                .debugDumpRhoBeforeSolve = true,
                .debugDumpScalarAfterSolve = true,
                .debugDumpVectorAfterSolve = true};
    }

    SolverCapabilities cgCapabilities() {
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

    [[noreturn]] void throwUnknownSolver(const std::string& context, const std::string& solverName) {
        throw OpalException(context, "No known solver matches the argument: " + solverName);
    }
}  // namespace

template <>
std::unique_ptr<PoissonBackend<double, 3>> PoissonBackendRegistry<double, 3>::create(
        const std::string& solverName, Solver_t<double, 3>& solverVariant, Field_t<3>* rho,
        VField_t<double, 3>* E, Field_t<3>* phi) {
    if (solverName == "NONE") {
        return std::make_unique<NullPoissonBackend<double, 3>>(solverVariant, rho, E, phi);
    }
    if (solverName == "FFT") {
        return std::make_unique<FFTPoissonBackend<double, 3>>(solverVariant, rho, E, phi);
    }
    if (solverName == "OPEN") {
        return std::make_unique<OpenPoissonBackend<double, 3>>(solverVariant, rho, E, phi);
    }
    if (solverName == "CG") {
        return std::make_unique<CGPoissonBackend<double, 3>>(solverVariant, rho, E, phi);
    }

    throwUnknownSolver("PoissonBackendRegistry::create", solverName);
}

template <>
SolverCapabilities PoissonBackendRegistry<double, 3>::capabilitiesFor(
        const std::string& solverName) {
    if (solverName == "NONE") {
        return noneCapabilities();
    }
    if (solverName == "FFT") {
        return fftCapabilities();
    }
    if (solverName == "OPEN") {
        return openCapabilities();
    }
    if (solverName == "CG") {
        return cgCapabilities();
    }

    throwUnknownSolver("PoissonBackendRegistry::capabilitiesFor", solverName);
}

template <>
double PoissonBackendRegistry<double, 3>::couplingConstantFor(const std::string& solverName) {
    if (solverName == "OPEN" || solverName == "FFT" || solverName == "CG") {
        return 1.0 / Physics::epsilon_0;
    }
    if (solverName == "NONE") {
        return 1.0 / (4.0 * Physics::pi * Physics::epsilon_0);
    }

    throwUnknownSolver("PoissonBackendRegistry::couplingConstantFor", solverName);
}

template <>
std::vector<std::string> PoissonBackendRegistry<double, 3>::solverNames() {
    return {"NONE", "FFT", "OPEN", "CG"};
}

template <>
bool PoissonBackendRegistry<double, 3>::supports(const std::string& solverName) {
    const auto names = solverNames();
    return std::find(names.begin(), names.end(), solverName) != names.end();
}
