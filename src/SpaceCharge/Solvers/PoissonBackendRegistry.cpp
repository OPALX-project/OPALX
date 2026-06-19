#include "SpaceCharge/Solvers/PoissonBackendRegistry.hpp"

#include <array>

#include "SpaceCharge/Solvers/PoissonBackends.hpp"
#include "Utilities/OpalException.h"

namespace {
    using BackendFactory = std::unique_ptr<PoissonBackend<double, 3>> (*)(
            Solver_t<double, 3>&, Field_t<3>*, VField_t<double, 3>*, Field_t<3>*);

    struct BackendSpec {
        const char* name;
        SolverCapabilities capabilities;
        double couplingConstant;
        BackendFactory create;
    };

    template <typename Backend>
    std::unique_ptr<PoissonBackend<double, 3>> makeBackend(
            Solver_t<double, 3>& solverVariant, Field_t<3>* rho, VField_t<double, 3>* E,
            Field_t<3>* phi) {
        return std::make_unique<Backend>(solverVariant, rho, E, phi);
    }

    const std::array<BackendSpec, 4>& backendSpecs() {
        static const std::array<BackendSpec, 4> specs = {{
                {NullPoissonBackend<double, 3>::solverName(),
                 NullPoissonBackend<double, 3>::backendCapabilities(),
                 NullPoissonBackend<double, 3>::backendCouplingConstant(),
                 &makeBackend<NullPoissonBackend<double, 3>>},
                {FFTPoissonBackend<double, 3>::solverName(),
                 FFTPoissonBackend<double, 3>::backendCapabilities(),
                 FFTPoissonBackend<double, 3>::backendCouplingConstant(),
                 &makeBackend<FFTPoissonBackend<double, 3>>},
                {OpenPoissonBackend<double, 3>::solverName(),
                 OpenPoissonBackend<double, 3>::backendCapabilities(),
                 OpenPoissonBackend<double, 3>::backendCouplingConstant(),
                 &makeBackend<OpenPoissonBackend<double, 3>>},
                {CGPoissonBackend<double, 3>::solverName(),
                 CGPoissonBackend<double, 3>::backendCapabilities(),
                 CGPoissonBackend<double, 3>::backendCouplingConstant(),
                 &makeBackend<CGPoissonBackend<double, 3>>},
        }};
        return specs;
    }

    const BackendSpec* findBackendSpec(const std::string& solverName) {
        for (const BackendSpec& spec : backendSpecs()) {
            if (solverName == spec.name) {
                return &spec;
            }
        }
        return nullptr;
    }

    [[noreturn]] void throwUnknownSolver(const std::string& context, const std::string& solverName) {
        throw OpalException(context, "No known solver matches the argument: " + solverName);
    }
}  // namespace

template <>
std::unique_ptr<PoissonBackend<double, 3>> PoissonBackendRegistry<double, 3>::create(
        const std::string& solverName, Solver_t<double, 3>& solverVariant, Field_t<3>* rho,
        VField_t<double, 3>* E, Field_t<3>* phi) {
    if (const BackendSpec* spec = findBackendSpec(solverName)) {
        return spec->create(solverVariant, rho, E, phi);
    }

    throwUnknownSolver("PoissonBackendRegistry::create", solverName);
}

template <>
SolverCapabilities PoissonBackendRegistry<double, 3>::capabilitiesFor(
        const std::string& solverName) {
    if (const BackendSpec* spec = findBackendSpec(solverName)) {
        return spec->capabilities;
    }

    throwUnknownSolver("PoissonBackendRegistry::capabilitiesFor", solverName);
}

template <>
double PoissonBackendRegistry<double, 3>::couplingConstantFor(const std::string& solverName) {
    if (const BackendSpec* spec = findBackendSpec(solverName)) {
        return spec->couplingConstant;
    }

    throwUnknownSolver("PoissonBackendRegistry::couplingConstantFor", solverName);
}

template <>
std::vector<std::string> PoissonBackendRegistry<double, 3>::solverNames() {
    std::vector<std::string> names;
    names.reserve(backendSpecs().size());
    for (const BackendSpec& spec : backendSpecs()) {
        names.emplace_back(spec.name);
    }
    return names;
}

template <>
bool PoissonBackendRegistry<double, 3>::supports(const std::string& solverName) {
    return findBackendSpec(solverName) != nullptr;
}
