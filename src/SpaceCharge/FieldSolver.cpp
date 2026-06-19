#include "SpaceCharge/FieldSolver.hpp"

#include "Utilities/OpalException.h"

extern Inform* gmsg;

template <>
void FieldSolver<double, 3>::dumpVectField(std::string what) {
    diagnostics_m.dumpVectorField(what, this->getE());
}

template <>
void FieldSolver<double, 3>::dumpScalField(std::string what) {
    diagnostics_m.dumpScalarField(what, this->getRho(), this->getPhi(), currentCapabilities());
}

template <>
void FieldSolver<double, 3>::setPotentialBCs() {
    Inform m("FieldSolver::setPotentialBCs");
    if (!hasValidBCHandler()) {
        throw OpalException("FieldSolver::setPotentialBCs", "BC Handler not set or invalid.");
    }

    const SolverCapabilities capabilities = currentCapabilities();
    if (!capabilities.requiresPotentialBCs) {
        m << level3 << "Potential BCs only need to be set for solvers that require them. Current "
          << "solver type: " << this->getStype() << endl;
        return;
    }

    // Need to do it like that, because for some reason IPPL wants a reference,
    // therefore cannot simply say "setFieldBC(...toIPPLBConds())".
    typedef ippl::BConds<Field_t<3>, 3> bc_type;
    bc_type bc_container = getBCHandler()->toIPPLBConds<Field_t<3>>();
    phi_m->setFieldBC(bc_container);
    m << level4 << "Potential BCs in FieldSolver updated using BCHandler." << endl;
}

template <>
void FieldSolver<double, 3>::initSolver() {
    Inform m("FieldSolver::initSolver");
    m << level3 << "Initializing backend for solver type: " << this->getStype() << endl;

    backend_m = PoissonBackendRegistry<double, 3>::create(
            this->getStype(), this->getSolver(), rho_m, E_m, phi_m);
    setPotentialBCs();
    diagnostics_m.resetCallCounter();
}

template <>
void FieldSolver<double, 3>::runSolver(
        const SolveRequest<double, 3>& request, bool force_skip_field_dump) {
    Inform m("FieldSolver::runSolver");
    m << level3 << "Running solver with type: " << this->getStype()
      << ". Force skip field dump: " << force_skip_field_dump << endl;

    if (!backend_m) {
        initSolver();
    }

    [[maybe_unused]] const SolverCapabilities capabilities = backend_m->capabilities();
    if (request.hasShiftedGreens() && !capabilities.supportsShiftedGreens) {
        throw OpalException(
                "FieldSolver::runSolver",
                "SHIFTED_GREENS_FUNCTION requires a field solver backend with shifted "
                "Green's-function support (got '" + this->getStype() + "').");
    }

#ifdef OPALX_FIELD_DEBUG
    if (!force_skip_field_dump && capabilities.debugDumpRhoBeforeSolve) {
        this->dumpScalField("rho");
    }
#endif

    backend_m->solve(request);

#ifdef OPALX_FIELD_DEBUG
    if (!force_skip_field_dump && capabilities.debugDumpScalarAfterSolve) {
        this->dumpScalField("phi");
    }
    if (!force_skip_field_dump && capabilities.debugDumpVectorAfterSolve) {
        this->dumpVectField("ef");
    }
#endif

    diagnostics_m.incrementCallCounter();
}

template <>
void FieldSolver<double, 3>::runSolver(bool force_skip_field_dump) {
    runSolver(SolveRequest<double, 3>{}, force_skip_field_dump);
}

template <>
double FieldSolver<double, 3>::getCouplingConstant() const {
    if (backend_m) {
        return backend_m->couplingConstant();
    }
    return PoissonBackendRegistry<double, 3>::couplingConstantFor(this->getStype());
}

template <>
void FieldSolver<double, 3>::refreshBackendAfterLayoutChange() {
    if (backend_m) {
        backend_m->refreshAfterLayoutChange();
    }
}
