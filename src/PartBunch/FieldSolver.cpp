#include "PartBunch/FieldSolver.hpp"

#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>

#include <filesystem>

#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include "SpaceCharge/SelfFieldDiagnostics.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

extern Inform* gmsg;

template <>
const opalx::spacecharge::IpplPoissonCapabilities& FieldSolver<double, 3>::getBackendCapabilities()
        const;

template <>
void FieldSolver<double, 3>::setPotentialBCs();

template <>
opalx::spacecharge::PoissonBackendKind FieldSolver<double, 3>::backendKindFromName(
        const std::string& solver) {
    const std::string value = Util::toUpper(solver);
    if (value == "NONE") {
        return opalx::spacecharge::PoissonBackendKind::None;
    }
    if (value == "FFT") {
        return opalx::spacecharge::PoissonBackendKind::FftPeriodic;
    }
    if (value == "OPEN") {
        return opalx::spacecharge::PoissonBackendKind::Open;
    }
    if (value == "CG") {
        return opalx::spacecharge::PoissonBackendKind::ConjugateGradient;
    }
    throw OpalException(
            "FieldSolver::backendKindFromName", "No known solver matches the argument: " + solver);
}

template <>
opalx::spacecharge::GreenFunctionKind FieldSolver<double, 3>::greenFunctionKindFromName(
        const std::string& greenFunction) {
    const std::string value = Util::toUpper(greenFunction);
    if (value == "STANDARD") {
        return opalx::spacecharge::GreenFunctionKind::Standard;
    }
    if (value == "INTEGRATED") {
        return opalx::spacecharge::GreenFunctionKind::Integrated;
    }
    throw OpalException(
            "FieldSolver::greenFunctionKindFromName",
            "Unknown GREENSF value \"" + greenFunction
                    + "\". Supported values are STANDARD and INTEGRATED.");
}

template <>
opalx::spacecharge::IpplPoissonFields FieldSolver<double, 3>::fields() {
    return {&workspace_m->chargeDensity(), &workspace_m->electricField(),
            &workspace_m->potential()};
}

template <>
void FieldSolver<double, 3>::dumpVectField(std::string what, std::size_t solveIndex) {
    /*
      what == ef
     */

    Inform m("FieldSolver::dumpVectorField");

    //    std::variant<Field_t<3>*, VField_t<double, 3>* > field;

    if (ippl::Comm->size() > 1 || solveIndex < 2) {
        m << level5 << "Skipping vector field dump for multiple ranks or first call." << endl;
        return;
    }

    /* Save the files in the output directory of the simulation. The file
     * name of vector fields is
     *
     * 'basename'-'name'_field-'******'.dat
     *
     * and of scalar fields
     *
     * 'basename'-'name'_scalar-'******'.dat
     *
     * with
     *   'basename': OPAL input file name (*.in)
     *   'name':     field name (input argument of function)
     *   '******':   solveIndex padded with zeros to 6 digits
     */

    std::string dirname = "data/";

    std::string type;
    std::string unit;

    if (Util::toUpper(what) == "EF") {
        type = "vector";
        unit = "";
    }

    VField_t<double, 3>* field = this->getE();

    auto localIdx = field->getOwned();
    auto mesh_mp  = &(field->get_mesh());
    auto spacing  = mesh_mp->getMeshSpacing();
    auto origin   = mesh_mp->getOrigin();
    int nghost    = field->getNghost();  // ghosts are excluded in getLocalNDIndex()

    auto fieldV      = field->getView();
    auto field_hostV = field->getHostMirror();
    Kokkos::deep_copy(field_hostV, fieldV);

    std::filesystem::path file(dirname);
    std::string basename = OpalData::getInstance()->getInputBasename();
    std::ostringstream oss;
    oss << basename << "-" << (what + std::string("_") + type) << "-" << std::setfill('0')
        << std::setw(6) << solveIndex << ".dat";
    std::string filename = oss.str();
    file /= filename;
    std::ofstream fout(file.string(), std::ios::out);

    fout << std::setprecision(9);

    fout << "# " << Util::toUpper(what) << " " << type << " data on grid" << std::endl
         << "# origin= " << std::fixed << origin << " h= " << std::fixed << spacing << std::endl
         << "#" << std::setw(4) << "i" << std::setw(5) << "j" << std::setw(5) << "k"
         << std::setw(17) << "x [m]" << std::setw(17) << "y [m]" << std::setw(17) << "z [m]";

    fout << std::setw(10) << what << "x [" << unit << "]" << std::setw(10) << what << "y [" << unit
         << "]" << std::setw(10) << what << "z [" << unit << "]";

    fout << std::endl;

    for (int i = localIdx[0].first() + nghost; i <= localIdx[0].last() + nghost; i++) {
        for (int j = localIdx[1].first() + nghost; j <= localIdx[1].last() + nghost; j++) {
            for (int k = localIdx[2].first() + nghost; k <= localIdx[2].last() + nghost; k++) {
                // define the physical points (cell-centered)
                double x = (i - nghost) * spacing[0] + origin[0];
                double y = (j - nghost) * spacing[1] + origin[1];
                double z = (k - nghost) * spacing[2] + origin[2];

                fout << std::setw(5) << i << std::setw(5) << j << std::setw(5) << k << std::setw(17)
                     << x << std::setw(17) << y << std::setw(17) << z << std::scientific << "\t"
                     << field_hostV(i, j, k)[0] << "\t" << field_hostV(i, j, k)[1] << "\t"
                     << field_hostV(i, j, k)[2] << std::endl;
            }
        }
    }
    fout.close();
    m << level5 << "*** FINISHED DUMPING " + Util::toUpper(what) + " FIELD *** to " << file.string()
      << endl;
}

template <>
void FieldSolver<double, 3>::dumpScalField(std::string what, std::size_t solveIndex) {
    /*
      what == phi | rho
     */

    Inform m("FS::dumpScalField() ");
    m << level5 << "Dumping scalar field: " << what << endl;

    if (ippl::Comm->size() > 1) {
        m << level5 << "Skipping scalar field dump for multiple ranks or first call." << endl;
        return;
    }

    /* Save the files in the output directory of the simulation. The file
     * name of vector fields is
     *
     * 'basename'-'name'_field-'******'.dat
     *
     * and of scalar fields
     *
     * 'basename'-'name'_scalar-'******'.dat
     *
     * with
     *   'basename': OPAL input file name (*.in)
     *   'name':     field name (input argument of function)
     *   '******':   solveIndex padded with zeros to 6 digits
     */

    // Needs to be empty...?
    std::string dirname = "data/";

    std::string type;
    std::string unit;
    bool isVectorField = false;

    if (Util::toUpper(what) == "RHO") {
        type = "scalar";
        unit = "Cb/m^3";
    } else if (Util::toUpper(what) == "PHI") {
        type = "scalar";
        unit = "V";
    }

    const bool separatePotential = getBackendCapabilities().usesSeparatePotentialField;
    Field_t<3>* field =
            (separatePotential && Util::toUpper(what) == "PHI") ? this->getPhi() : this->getRho();

    // auto localIdx = field->getOwned();
    ippl::NDIndex<3> localIdx = field->getLayout().getLocalNDIndex();
    int nghost = field->getNghost();  // ghosts are excluded in getLocalNDIndex(), but we still need
                                      // to shift indices
    auto mesh_mp = &(field->get_mesh());
    auto spacing = mesh_mp->getMeshSpacing();
    auto origin  = mesh_mp->getOrigin();

    Field_t<3>::view_type fieldV             = field->getView();
    Field_t<3>::host_mirror_type field_hostV = field->getHostMirror();
    Kokkos::deep_copy(field_hostV, fieldV);

    std::filesystem::path file(dirname);
    std::string basename = OpalData::getInstance()->getInputBasename();
    std::ostringstream oss;
    oss << basename << "-" << (what + std::string("_") + type) << "-" << std::setfill('0')
        << std::setw(6) << solveIndex << ".dat";
    std::string filename = oss.str();
    file /= filename;
    std::ofstream fout(file.string(), std::ios::out);

    fout << std::setprecision(9);

    fout << "# " << Util::toUpper(what) << " " << type << " data on grid" << std::endl
         << "# origin= " << std::fixed << origin << " h= " << std::fixed << spacing
         << " nghosts=" << nghost << std::endl
         << "#" << std::setw(4) << "i" << std::setw(5) << "j" << std::setw(5) << "k"
         << std::setw(17) << "x [m]" << std::setw(17) << "y [m]" << std::setw(17) << "z [m]";

    if (isVectorField) {
        fout << std::setw(10) << what << "x [" << unit << "]" << std::setw(10) << what << "y ["
             << unit << "]" << std::setw(10) << what << "z [" << unit << "]";
    } else {
        fout << std::setw(13) << what << " [" << unit << "]";
    }

    fout << std::endl;

    if (Util::toUpper(what) == "RHO") {
        for (int i = localIdx[0].first() + nghost; i <= localIdx[0].last() + nghost; i++) {
            for (int j = localIdx[1].first() + nghost; j <= localIdx[1].last() + nghost; j++) {
                for (int k = localIdx[2].first() + nghost; k <= localIdx[2].last() + nghost; k++) {
                    // define the physical points (cell-centered)
                    double x = (i - nghost) * spacing[0] + origin[0];
                    double y = (j - nghost) * spacing[1] + origin[1];
                    double z = (k - nghost) * spacing[2] + origin[2];

                    fout << std::setw(5) << i << std::setw(5) << j << std::setw(5) << k
                         << std::setw(17) << x << std::setw(17) << y << std::setw(17) << z
                         << std::scientific << "\t" << field_hostV(i, j, k) << std::endl;
                }
            }
        }
    } else {
        for (int i = localIdx[0].first() + nghost; i <= localIdx[0].last() + nghost; i++) {
            for (int j = localIdx[1].first() + nghost; j <= localIdx[1].last() + nghost; j++) {
                for (int k = localIdx[2].first() + nghost; k <= localIdx[2].last() + nghost; k++) {
                    // define the physical points (cell-centered)
                    double x = (i - nghost) * spacing[0] + origin[0];
                    double y = (j - nghost) * spacing[1] + origin[1];
                    double z = (k - nghost) * spacing[2] + origin[2];

                    // "+ 1" matches OPAL indexing in the output
                    fout << std::setw(5) << i << std::setw(5) << j << std::setw(5) << k
                         << std::setw(17) << x << std::setw(17) << y << std::setw(17) << z
                         << std::scientific << "\t" << field_hostV(i, j, k) << std::endl;
                }
            }
        }
    }
    fout.close();
    m << level5 << "*** FINISHED DUMPING " + Util::toUpper(what) + " FIELD *** to " << file.string()
      << endl;
}

template <>
void FieldSolver<double, 3>::initSolver() {
    Inform m("FieldSolver::initSolver");
    backend_m = opalx::spacecharge::IpplPoissonAdapter::create(
            backendKind_m, greenFunctionKind_m, this->getSolver(), fields());
    setPotentialBCs();
    m << level3 << "Initialized typed IPPL backend for solver: " << this->getStype() << endl;
}

template <>
void FieldSolver<double, 3>::refreshAfterFieldLayoutChange() {
    Inform m("FieldSolver::refreshAfterFieldLayoutChange");
    m << level4
      << "Refreshing existing solver backend for field layout change: " << this->getStype() << endl;

    if (backend_m == nullptr) {
        throw OpalException(
                "FieldSolver::refreshAfterFieldLayoutChange",
                "The IPPL backend must be initialized before it is refreshed.");
    }
    backend_m->refresh(fields());
}

template <>
void FieldSolver<double, 3>::runSolverImpl(
        const opalx::spacecharge::IpplPoissonSolveRequest& request, bool force_skip_field_dump,
        opalx::spacecharge::SelfFieldDiagnostics* diagnostics);

template <>
void FieldSolver<double, 3>::runSolver(bool force_skip_field_dump) {
    runSolverImpl({}, force_skip_field_dump, nullptr);
}

template <>
void FieldSolver<double, 3>::runSolver(
        bool force_skip_field_dump, opalx::spacecharge::SelfFieldDiagnostics& diagnostics) {
    runSolverImpl({}, force_skip_field_dump, &diagnostics);
}

template <>
void FieldSolver<double, 3>::runSolverImpl(
        const opalx::spacecharge::IpplPoissonSolveRequest& request, bool force_skip_field_dump,
        opalx::spacecharge::SelfFieldDiagnostics* diagnostics) {
    Inform m("FieldSolver::runSolver");
    m << level3 << "Running solver with type: " << this->getStype()
      << ". Force skip field dump: " << force_skip_field_dump << endl;

    if (backend_m == nullptr) {
        throw OpalException("FieldSolver::runSolver", "The typed IPPL backend is not initialized.");
    }

    const std::size_t solveIndex = diagnostics != nullptr ? diagnostics->backendSolveCount() : 0;
    static_cast<void>(solveIndex);
    std::optional<opalx::spacecharge::SelfFieldDiagnostics::ScopedEvent> backendEvent;
    if (diagnostics != nullptr) {
        backendEvent.emplace(
                *diagnostics, opalx::spacecharge::SelfFieldEventKind::BackendSolve,
                request.hasShiftedGreenFunction() ? "OPEN-shifted" : "backend");
    }

    [[maybe_unused]] const auto& capabilities = backend_m->capabilities();
#ifdef OPALX_FIELD_DEBUG
    if (!force_skip_field_dump && capabilities.debugDumpChargeBeforeSolve) {
        this->dumpScalField("rho", solveIndex);
    }
#endif

    backend_m->solve(request);

#ifdef OPALX_FIELD_DEBUG
    if (!force_skip_field_dump && capabilities.debugDumpScalarAfterSolve) {
        this->dumpScalField("phi", solveIndex);
    }
    if (!force_skip_field_dump && capabilities.debugDumpVectorAfterSolve) {
        this->dumpVectField("ef", solveIndex);
    }
#endif

    if (diagnostics != nullptr) {
        diagnostics->completeBackendSolve();
    }
}

template <>
void FieldSolver<double, 3>::runShiftedOpenSolver(
        const ippl::Vector<double, 3>& shift,
        opalx::spacecharge::SelfFieldDiagnostics& diagnostics) {
    if (!getBackendCapabilities().supportsShiftedGreenFunction) {
        throw OpalException(
                "FieldSolver::runShiftedOpenSolver",
                "SHIFTED_GREENS_FUNCTION requires FIELDSOLVER type OPEN (got '" + this->getStype()
                        + "').");
    }

    Inform m("FieldSolver::runShiftedOpenSolver");
    m << level4 << "Running shifted open solver with shift = " << shift << endl;

    opalx::spacecharge::IpplPoissonSolveRequest request;
    request.greenFunctionShift = shift;
    runSolverImpl(request, false, &diagnostics);
}

template <>
double FieldSolver<double, 3>::getCouplingConstant() const {
    return backend_m != nullptr
                   ? backend_m->couplingConstant()
                   : opalx::spacecharge::IpplPoissonAdapter::couplingConstantFor(backendKind_m);
}

template <>
const opalx::spacecharge::IpplPoissonCapabilities& FieldSolver<double, 3>::getBackendCapabilities()
        const {
    return backend_m != nullptr
                   ? backend_m->capabilities()
                   : opalx::spacecharge::IpplPoissonAdapter::capabilitiesFor(backendKind_m);
}

template <>
void FieldSolver<double, 3>::setPotentialBCs() {
    Inform m("FieldSolver::setPotentialBCs");
    // Check if BC handler is set
    if (!hasValidBCHandler()) {
        throw OpalException("FieldSolver::setPotentialBCs", "BC Handler not set or invalid.");
    }

    if (!getBackendCapabilities().requiresPotentialBoundaryConditions) {
        m << level3 << "Potential BCs only need to be set for CG solver. Current solver type: "
          << this->getStype() << endl;
        return;
    }

    // Need to do it like that, because for some reason IPPL wants a reference,
    // therefore cannot simply say "setFieldBC(...toIPPLBConds())".
    typedef ippl::BConds<Field_t<Dim>, Dim> bc_type;
    bc_type bc_container = getBCHandler()->toIPPLBConds<Field_t<Dim>>();
    getPhi()->setFieldBC(bc_container);
    m << level4 << "Potential BCs in FieldSolver updated using BCHandler." << endl;
}
