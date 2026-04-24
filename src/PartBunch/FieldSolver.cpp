#include "PartBunch/FieldSolver.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>

#include <filesystem>

#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include "Utilities/Util.h"

extern Inform* gmsg;

template <>
template <typename Solver>
void FieldSolver<double, 3>::initSolverWithParams(const ippl::ParameterList& sp) {
    Inform m("FieldSolver::initSolverWithParams");
    this->getSolver().template emplace<Solver>();
    Solver& solver = std::get<Solver>(this->getSolver());
    solver.mergeParameters(sp);
    m << level3 << "Initialized solver with params: " << this->getStype() << endl;

    // test if rho_m exists (just in case)
    if (!rho_m) {
        throw OpalException("FieldSolver::initSolverWithParams", "rho_m is not initialized.");
    }

    solver.setRhs(*rho_m);
    m << level5 << "Set solver RHS." << endl;

    if constexpr ((std::is_same_v<Solver, CGSolver_t<T, Dim>>) /*|| 
                  (std::is_same_v<Solver, FEMSolver_t<T, Dim>>) || 
                  (std::is_same_v<Solver, FEMPreconSolver_t<T, Dim>>)*/) {
        /// \todo for now, don't use the CG solver!
        throw OpalException(
                "FieldSolver::initSolverWithParams",
                "Cannot use CGSolver yet, not fully implemented.");
        // The CG solver and FEMPoissonSolver compute the potential
        // directly and use this to get the electric field
        solver.setLhs(*phi_m);
        solver.setGradient(*E_m);
        m << level5 << "Set gradient for CG." << endl;
    } else {
        // The periodic Poisson solver, Open boundaries solver,
        // and the TG solver compute the electric field directly
        solver.setLhs(*E_m);
    }

    m << level4 << "Set solver RHS+LHS." << endl;
    call_counter_m = 0;
}

template <>
void FieldSolver<double, 3>::dumpVectField(std::string what) {
    /*
      what == ef
     */

    Inform m("FieldSolver::dumpVectorField");

    //    std::variant<Field_t<3>*, VField_t<double, 3>* > field;

    if (ippl::Comm->size() > 1 || call_counter_m < 2) {
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
     *   '******':   call_counter_m padded with zeros to 6 digits
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
        << std::setw(6) << call_counter_m << ".dat";
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
void FieldSolver<double, 3>::dumpScalField(std::string what) {
    /*
      what == phi | rho
     */

    Inform m("FS::dumpScalField() ");
    m << level5 << "Dumping scalar field: " << what << endl;

    if (ippl::Comm->size() > 1 /*|| call_counter_m<2*/) {
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
     *   '******':   call_counter_m padded with zeros to 6 digits
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

    Field_t<3>* field = (this->getStype() == "CG" && Util::toUpper(what) == "PHI")
                                ? this->getPhi()
                                : this->getRho();  // both rho and phi are in the same variable (in
                                                   // place computation)

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
        << std::setw(6) << call_counter_m << ".dat";
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
void FieldSolver<double, 3>::initOpenSolver() {
    ippl::ParameterList sp;
    sp.add("output_type", OpenSolver_t<double, 3>::SOL_AND_GRAD);
    sp.add("use_heffte_defaults", false);
    sp.add("use_pencils", true);
    sp.add("use_reorder", false);
    sp.add("use_gpu_aware", true);
    sp.add("comm", ippl::p2p_pl);
    sp.add("r2c_direction", 0);
    sp.add("algorithm", OpenSolver_t<double, 3>::HOCKNEY);
    initSolverWithParams<OpenSolver_t<double, 3>>(sp);
}

template <>
void FieldSolver<double, 3>::initFFTSolver() {
    if constexpr (Dim == 2 || Dim == 3) {
        ippl::ParameterList sp;
        // Needs sol for phi_m testing/output and GRAD for E_m computation
        /// \todo don't print phi_m to file if FFT and GRAD is selected!
        sp.add("output_type", FFTSolver_t<double, 3>::GRAD);
        // sp.add("output_type", OpenSolver_t<double, 3>::SOL_AND_GRAD);
        sp.add("use_heffte_defaults", false);
        sp.add("use_pencils", true);
        sp.add("use_reorder", false);
        sp.add("use_gpu_aware", true);
        sp.add("comm", ippl::p2p_pl);
        sp.add("r2c_direction", 0);
        initSolverWithParams<FFTSolver_t<double, 3>>(sp);
    } else {
        throw OpalException(
                "FieldSolver<double,3>::initFFTSolver",
                "FFTSolver_t is only implemented for 2D and 3D fields.");
    }
}

template <>
void FieldSolver<double, 3>::initCGSolver() {
    ippl::ParameterList sp;
    sp.add("output_type", CGSolver_t<double, 3>::GRAD);
    /// \todo should probably at some point be passed from the input file
    sp.add("tolerance", 1e-12);

    initSolverWithParams<CGSolver_t<double, 3>>(sp);
}

template <>
void FieldSolver<double, 3>::initNullSolver() {
    ippl::ParameterList sp;
    if constexpr (Dim == 2 || Dim == 3) {
        initSolverWithParams<NullSolver_t<T, Dim>>(sp);
    } else {
        throw OpalException(
                "FieldSolver<double,3>::initNullSolver",
                "NullSolver_t is only implemented for 2D and 3D fields.");
    }
}

template <>
void FieldSolver<double, 3>::initSolver() {
    Inform m;
    if (this->getStype() == "FFT") {
        initFFTSolver();
    } else if (this->getStype() == "OPEN") {
        initOpenSolver();
    } else if (this->getStype() == "CG") {
        initCGSolver();
    } else if (this->getStype() == "NONE") {
        initNullSolver();
    } else {
        throw OpalException(
                "FieldSolver::initSolver",
                "No known solver matches the argument: " + this->getStype());
    }
}

/*template <>
void FieldSolver<double,3>::setPotentialBCs() {
        // CG requires explicit periodic boundary conditions while the periodic Poisson solver
        // simply assumes them
        if (this->getStype() == "CG") {
            typedef ippl::BConds<Field<double, 3>, 3> bc_type;
            bc_type allPeriodic;
            for (unsigned int i = 0; i < 2 * 3; ++i) {
                allPeriodic[i] = std::make_shared<ippl::PeriodicFace<Field<double, 3>>>(i);
            }
            phi_m->setFieldBC(allPeriodic);
        }
    }*/

template <>
void FieldSolver<double, 3>::runSolver(bool force_skip_field_dump) {
    constexpr int Dim = 3;
    Inform m("FieldSolver::runSolver");
    /*
    Add this output such that there is no possible unused variable warning for
    force_skip_field_dump, which is only used in the debug output functions
    for now.
    */
    m << level3 << "Running solver with type: " << this->getStype()
      << ". Force skip field dump: " << force_skip_field_dump << endl;

    if (this->getStype() == "CG") {
        CGSolver_t<double, 3>& solver = std::get<CGSolver_t<double, 3>>(this->getSolver());
#ifdef OPALX_FIELD_DEBUG
        if (!force_skip_field_dump) this->dumpScalField("rho");
#endif
        // std::get<CGSolver_t<double, 3>>(this->getSolver()).solve();
        solver.solve();
#ifdef OPALX_FIELD_DEBUG
        if (!force_skip_field_dump) this->dumpScalField("phi");
#endif
        /*
        if (ippl::Comm->rank() == 0) {
            std::stringstream fname;
            fname << "data/CG_";
            fname << ippl::Comm->size();
            fname << ".csv";

            Inform log(NULL, fname.str().c_str(), Inform::APPEND);
            int iterations = solver.getIterationCount();
            // Assume the dummy solve is the first call
            if (iterations == 0) {
                log << "residue,iterations" << endl;
            }
            // Don't print the dummy solve
            if (iterations > 0) {
                log << solver.getResidue() << "," << iterations << endl;
            }
        }
        ippl::Comm->barrier();*/
        int iterations = solver.getIterationCount();
        int residue    = solver.getResidue();
        m << level4 << "CG solver finished. Iterations: " << iterations << ", Residue: " << residue
          << endl;
    } else if (this->getStype() == "FFT") {
        if constexpr (Dim == 2 || Dim == 3) {
#ifdef OPALX_FIELD_DEBUG
            if (!force_skip_field_dump) this->dumpScalField("rho");
#endif

            std::get<FFTSolver_t<double, 3>>(this->getSolver()).solve();
#ifdef OPALX_FIELD_DEBUG
            /// \todo do to not print phi/E depenging on output_type!
            if (!force_skip_field_dump) this->dumpScalField("phi");
            if (!force_skip_field_dump) this->dumpVectField("ef");
#endif
        }
    } else if (this->getStype() == "P3M") {
        if constexpr (Dim == 3) {
            std::get<FFTTruncatedGreenSolver_t<double, 3>>(this->getSolver()).solve();
        }
    } else if (this->getStype() == "OPEN") {
        if constexpr (Dim == 3) {
#ifdef OPALX_FIELD_DEBUG
            if (!force_skip_field_dump) this->dumpScalField("rho");
#endif
            std::get<OpenSolver_t<double, 3>>(this->getSolver()).solve();
#ifdef OPALX_FIELD_DEBUG
            if (!force_skip_field_dump) this->dumpScalField("phi");
            if (!force_skip_field_dump) this->dumpVectField("ef");
#endif
        }
    } else if (this->getStype() == "NONE") {
        std::get<NullSolver_t<T, Dim>>(this->getSolver()).solve();
    } else {
        throw OpalException(
                "FieldSolver::runSolver",
                "No known solver matches the argument: " + this->getStype());
    }

    call_counter_m++;
}

template <>
void FieldSolver<double, 3>::runShiftedOpenSolver(const ippl::Vector<double, 3>& shift) {
    if (this->getStype() != "OPEN") {
        throw OpalException(
                "FieldSolver::runShiftedOpenSolver",
                "SHIFTED_GREENS_FUNCTION requires FIELDSOLVER type OPEN (got '" + this->getStype()
                        + "').");
    }

    Inform m("FieldSolver::runShiftedOpenSolver");
    m << level4 << "Running shifted open solver with shift = " << shift << endl;

    auto& openSolver = std::get<OpenSolver_t<double, 3>>(this->getSolver());

    // Install the translated kernel (overwrites grntr_m inside IPPL).
    openSolver.shiftedGreensFunction(shift);

#ifdef OPALX_FIELD_DEBUG
    this->dumpScalField("rho");
#endif

    openSolver.solve();

#ifdef OPALX_FIELD_DEBUG
    this->dumpScalField("phi");
    this->dumpVectField("ef");
#endif

    // Restore the standard Green's function so subsequent primary solves in
    // later bins do not silently reuse the shifted kernel. See the header
    // doc for the defensive-guard rationale.
    openSolver.greensFunction();

    call_counter_m++;
}

// Implement getCouplingConstant
template <>
double FieldSolver<double, 3>::getCouplingConstant() const {
    /*
    In SI units, the coupling constant for the electric field is
    1/(4*pi*epsilon_0). However, some solvers seem to use different conventions
    (likely due to different Green's function conventions or FFT normalizations).

    As I can tell, the IPPL solvers all need 1/eps0. However, just in case,
    leave it like this for now. Perhaps later this can be changed to simply
    something like `return 1.0 / Physics::epsilon_0;`.
    */

    /// \todo Verify this before activating a new solver!
    const std::string stype = this->getStype();
    if (stype == "OPEN") {
        return 1.0 / Physics::epsilon_0;
    } else if (stype == "FFT") {
        return 1.0 / Physics::epsilon_0;
    } else if (stype == "CG") {
        return 1.0 / Physics::epsilon_0;
    }

    // Standard coupling constant (from before)
    return 1.0 / (4.0 * Physics::pi * Physics::epsilon_0);
}

template <>
void FieldSolver<double, 3>::setPotentialBCs() {
    Inform m("FieldSolver::setPotentialBCs");
    // Check if BC handler is set
    if (!hasValidBCHandler()) {
        throw OpalException("FieldSolver::setPotentialBCs", "BC Handler not set or invalid.");
    }

    if (this->getStype() != "CG") {
        m << level3 << "Potential BCs only need to be set for CG solver. Current solver type: "
          << this->getStype() << endl;
        return;
    }

    // Need to do it like that, because for some reason IPPL wants a reference,
    // therefore cannot simply say "setFieldBC(...toIPPLBConds())".
    typedef ippl::BConds<Field_t<Dim>, Dim> bc_type;
    bc_type bc_container = getBCHandler()->toIPPLBConds<Field_t<Dim>>();
    phi_m->setFieldBC(bc_container);
    // rho_m->setFieldBC(bc_container);
    m << level4 << "Potential BCs in FieldSolver updated using BCHandler." << endl;
}

/*
/// \todo to be implemented...
template<>
void FieldSolver<double,3>::initP3MSolver() {
    //        if constexpr (Dim == 3) {
    ippl::ParameterList sp;
    sp.add("output_type", P3MSolver_t<double, 3>::GRAD);
    sp.add("use_heffte_defaults", false);
    sp.add("use_pencils", true);
    sp.add("use_reorder", false);
    sp.add("use_gpu_aware", true);
    sp.add("comm", ippl::p2p_pl);
    sp.add("r2c_direction", 0);

    initSolverWithParams<P3MSolver_t<double, 3>>(sp);
    //  } else {
    // throw std::runtime_error("Unsupported dimensionality for P3M solver");
    // }
}

*/
