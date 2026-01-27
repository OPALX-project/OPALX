#include "PartBunch/FieldSolver.hpp"

#include <iomanip>
#include <fstream>
#include <sstream>

#include <filesystem>

#include "Utilities/Util.h"
#include "AbstractObjects/OpalData.h"

extern Inform* gmsg;

template <>
template <typename Solver>
void  FieldSolver<double,3>::initSolverWithParams(const ippl::ParameterList& sp) {
    Inform m ("initSolverWithParams");
    this->getSolver().template emplace<Solver>();
    Solver& solver = std::get<Solver>(this->getSolver());
    solver.mergeParameters(sp);

    // test if rho_m exists (just in case)
    if (!rho_m) {
        throw OpalException("FieldSolver<double,3>::initSolverWithParams", "rho_m is null pointer.");
    }

    m << "Init solver with params: " << solver.getStype() << endl;
    solver.setRhs(*rho_m);
    m << "Set solver RHS" << endl;

    if constexpr ((std::is_same_v<Solver, CGSolver_t<T, Dim>>) || 
                  (std::is_same_v<Solver, FEMSolver_t<T, Dim>>) || 
                  (std::is_same_v<Solver, FEMPreconSolver_t<T, Dim>>)) {
        // The CG solver and FEMPoissonSolver compute the potential 
        // directly and use this to get the electric field
        solver.setLhs(*phi_m);
        solver.setGradient(*E_m);
        m << "Set gradient for CG/FEM solver" << endl;
    } else {
        // The periodic Poisson solver, Open boundaries solver,
        // and the TG solver compute the electric field directly
        solver.setLhs(*E_m);
    }
    m << "Set solver LHS" << endl;
    
    /// \todo This can potentially be implemented much easier, see https://github.com/IPPL-framework/ippl/blob/f4c4102a3cb76dd2cd911eef7314adb77aacd676/alpine/FieldSolver.hpp#L148
    /*if constexpr (std::is_same_v<Solver, CGSolver_t<T, Dim>>) {
        // The CG solver computes the potential directly and
        // uses this to get the electric field
        throw OpalException("FieldSolver<double,3>::initSolverWithParams", "Cannot use CGSolver yet, not implemented.");

        /// \todo implement properly
        m << "CG solver used " << endl;
        solver.setLhs(*phi_m);
        solver.setGradient(*E_m);
    } else if constexpr (std::is_same_v<Solver, OpenSolver_t<T, Dim>>) {
        // The periodic Poisson solver, Open boundaries solver,
        // and the P3M solver compute the electric field directly
        m << "OpenSolver used" << endl;

        solver.setLhs(*E_m);
        m << "Set LHS for OpenSolver" << endl;
    } else if constexpr (std::is_same_v<Solver, FFTSolver_t<T, Dim>>) {
        // The periodic Poisson solver
        m << "FFTSolver used" << endl;
        
        solver.setLhs(*E_m);
        m << "Set LHS for FFTSolver" << endl;
    } else if constexpr (std::is_same_v<Solver, NullSolver_t<T, Dim>>) {
        m << "NullSolver used" << endl;
        solver.setLhs(*E_m);
    }*/
    
    call_counter_m = 0;
}


template <>
void FieldSolver<double,3>::dumpVectField(std::string what) {
    /*
      what == ef
     */

    Inform m("FS::dumpVectorField() ");

    //    std::variant<Field_t<3>*, VField_t<double, 3>* > field;

    if (ippl::Comm->size() > 1 || call_counter_m<2) {
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
    int nghost    = field->getNghost(); // ghosts are excluded in getLocalNDIndex()

    auto fieldV      = field->getView();
    auto field_hostV = field->getHostMirror();
    Kokkos::deep_copy(field_hostV, fieldV);     

    std::filesystem::path file(dirname);
    std::string basename = OpalData::getInstance()->getInputBasename();
    std::ostringstream oss;
    oss << basename << "-" << (what + std::string("_") + type) << "-" 
        << std::setfill('0') << std::setw(6) << call_counter_m << ".dat";
    std::string filename = oss.str();
    file /= filename;
    std::ofstream fout(file.string(), std::ios::out);

    fout << std::setprecision(9);

    fout << "# " << Util::toUpper(what) << " " << type << " data on grid" << std::endl
         << "# origin= " << std::fixed << origin << " h= " << std::fixed << spacing << std::endl 
         << "#"
         << std::setw(4)  << "i"
         << std::setw(5)  << "j"
         << std::setw(5)  << "k"
         << std::setw(17) << "x [m]"
         << std::setw(17) << "y [m]"
         << std::setw(17) << "z [m]";

    fout << std::setw(10) << what << "x [" << unit << "]"
         << std::setw(10) << what << "y [" << unit << "]"
         << std::setw(10) << what << "z [" << unit << "]";

    fout << std::endl;

    for (int i = localIdx[0].first() + nghost; i <= localIdx[0].last() + nghost; i++) {
            for (int j = localIdx[1].first() + nghost; j <= localIdx[1].last() + nghost ; j++) {
                for (int k = localIdx[2].first() + nghost; k <= localIdx[2].last() + nghost; k++) {
                    
                    // define the physical points (cell-centered)
                    double x = (i - nghost) * spacing[0] + origin[0];        
                    double y = (j - nghost) * spacing[1] + origin[1];        
                    double z = (k - nghost) * spacing[2] + origin[2];     
                
                    fout << std::setw(5) << i 
                         << std::setw(5) << j
                         << std::setw(5) << k
                         << std::setw(17) << x
                         << std::setw(17) << y
                         << std::setw(17) << z
                         << std::scientific 
                         << "\t" << field_hostV(i,j,k)[0]
                         << "\t" << field_hostV(i,j,k)[1]
                         << "\t" << field_hostV(i,j,k)[2]                      
                         << std::endl;
                }
            }
        }
    fout.close();
    m << "*** FINISHED DUMPING " + Util::toUpper(what) + " FIELD *** to " << file.string() << endl;
}

template <>
void FieldSolver<double,3>::dumpScalField(std::string what) {

    /*
      what == phi | rho
     */

    Inform m("FS::dumpScalField() ");
    m << "Dumping scalar field: " << what << endl;

    if (ippl::Comm->size() > 1 /*|| call_counter_m<2*/) {
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

    
    Field_t<3>* field = this->getRho();   // both rho and phi are in the same variable (in place computation)
    
    // auto localIdx = field->getOwned();
    ippl::NDIndex<3> localIdx = field->getLayout().getLocalNDIndex();
    int nghost    = field->getNghost(); // ghosts are excluded in getLocalNDIndex(), but we still need to shift indices
    auto mesh_mp  = &(field->get_mesh());
    auto spacing  = mesh_mp->getMeshSpacing();
    auto origin   = mesh_mp->getOrigin();

    Field_t<3>::view_type fieldV       = field->getView();
    Field_t<3>::HostMirror field_hostV = field->getHostMirror();
    Kokkos::deep_copy(field_hostV, fieldV);     

    std::filesystem::path file(dirname);
    std::string basename = OpalData::getInstance()->getInputBasename();
    std::ostringstream oss;
    oss << basename << "-" << (what + std::string("_") + type) << "-" 
        << std::setfill('0') << std::setw(6) << call_counter_m << ".dat";
    std::string filename = oss.str();
    file /= filename;
    std::ofstream fout(file.string(), std::ios::out);

    fout << std::setprecision(9);

    fout << "# " << Util::toUpper(what) << " " << type << " data on grid" << std::endl
         << "# origin= " << std::fixed << origin 
         << " h= " << std::fixed << spacing 
         << " nghosts=" << nghost << std::endl 
         << "#"
         << std::setw(4)  << "i"
         << std::setw(5)  << "j"
         << std::setw(5)  << "k"
         << std::setw(17) << "x [m]"
         << std::setw(17) << "y [m]"
         << std::setw(17) << "z [m]";

    if (isVectorField) {
        fout << std::setw(10) << what << "x [" << unit << "]"
             << std::setw(10) << what << "y [" << unit << "]"
             << std::setw(10) << what << "z [" << unit << "]";
    } else {
        fout << std::setw(13) << what << " [" << unit << "]";
    }

    fout << std::endl;

    if (Util::toUpper(what) == "RHO") {
        for (int i = localIdx[0].first() + nghost; i <= localIdx[0].last() + nghost; i++) {
            for (int j = localIdx[1].first() + nghost; j <= localIdx[1].last() + nghost; j++) {
                for (int k = localIdx[2].first() + nghost; k <= localIdx[2].last() + nghost; k++) {
                    
                    // define the physical points (cell-centered)
                    double x = (i-nghost) * spacing[0] + origin[0];        
                    double y = (j-nghost) * spacing[1] + origin[1];        
                    double z = (k-nghost) * spacing[2] + origin[2];     
                
                    fout << std::setw(5) << i
                         << std::setw(5) << j
                         << std::setw(5) << k
                         << std::setw(17) << x
                         << std::setw(17) << y
                         << std::setw(17) << z
                         << std::scientific << "\t" << field_hostV(i,j,k)                             
                         << std::endl;
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
                    fout << std::setw(5) << i
                         << std::setw(5) << j
                         << std::setw(5) << k
                         << std::setw(17) << x
                         << std::setw(17) << y
                         << std::setw(17) << z
                         << std::scientific << "\t" << field_hostV(i,j,k)                             
                         << std::endl;
                }
            }
        }
    }
    fout.close();
    m << "*** FINISHED DUMPING " + Util::toUpper(what) + " FIELD *** to " << file.string() << endl;
}

template <>
void FieldSolver<double,3>::initOpenSolver() {
    ippl::ParameterList sp;
    sp.add("output_type", OpenSolver_t<double, 3>::SOL_AND_GRAD); // see https://github.com/IPPL-framework/ippl/blob/f4c4102a3cb76dd2cd911eef7314adb77aacd676/alpine/FieldSolver.hpp#L319
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
void FieldSolver<double,3>::initFFTSolver() {
    if constexpr (Dim == 2 || Dim == 3) {
        ippl::ParameterList sp;
        sp.add("output_type", FFTSolver_t<double, 3>::GRAD);
        //sp.add("output_type", OpenSolver_t<double, 3>::SOL_AND_GRAD);
        sp.add("use_heffte_defaults", false);
        sp.add("use_pencils", true);
        sp.add("use_reorder", false);
        sp.add("use_gpu_aware", true);
        sp.add("comm", ippl::p2p_pl);
        sp.add("r2c_direction", 0);
        initSolverWithParams<FFTSolver_t<double, 3>>(sp);
    } else {
        // TODO: add exception here
        // throw std::runtime_error("Unsupported dimensionality for FFT solver");
    }
}

template <>
void FieldSolver<double,3>::initCGSolver() {
    ippl::ParameterList sp;
    sp.add("output_type", CGSolver_t<double, 3>::GRAD);
    // Increase tolerance in the 1D case
    sp.add("tolerance", 1e-4);
    
    initSolverWithParams<CGSolver_t<double, 3>>(sp);
}

template<>
void FieldSolver<double,3>::initNullSolver() {
    ippl::ParameterList sp;
    if constexpr (Dim == 2 || Dim == 3) {
        initSolverWithParams<NullSolver_t<T, Dim>>(sp);
    } else {
        throw std::runtime_error("Unsupported dimensionality for Null solver");
    }
}

template <>
void FieldSolver<double,3>::initSolver() {
    Inform m;
    if (this->getStype() == "FFT") {
        initFFTSolver();    
    } else if (this->getStype() == "OPEN") {
        initOpenSolver();    
    } else if (this->getStype() == "CG") {
        initCGSolver();
    } else if (this->getStype() == "NONE") {
        initNullSolver();
    }
    else {
        m << "No solver matches the argument: " << this->getStype() << endl;
        throw std::runtime_error("No solver match");
    }
}

template <>
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
    }

template<>
void FieldSolver<double,3>::runSolver(bool force_skip_field_dump) {
    constexpr int Dim = 3;

    if (this->getStype() == "CG") {
            // CGSolver_t<double, 3>& solver = std::get<CGSolver_t<double, 3>>(this->getSolver());
            // solver.solve();
            std::get<CGSolver_t<double, 3>>(this->getSolver()).solve();
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
    } else if (this->getStype() == "FFT") {
        if constexpr (Dim == 2 || Dim == 3) {
#ifdef OPALX_FIELD_DEBUG
            if (!force_skip_field_dump) this->dumpScalField("rho");
#endif

            std::get<FFTSolver_t<double, 3>>(this->getSolver()).solve();
#ifdef OPALX_FIELD_DEBUG
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
        throw std::runtime_error("Unknown solver type");
    }

    call_counter_m++; // maybe want "if (!force_skip_field_dump)" here?
}

// Implement getCouplingConstant
template<>
double FieldSolver<double, 3>::getCouplingConstant() const {
    /*
    In SI units, the coupling constant for the electric field is 
    1/(4*pi*epsilon_0). However, some solvers seem to use different conventions
    (likely due to different Green's function conventions or FFT normalizations). 
    */

    /// \todo Verify this before activating a new solver!
    const std::string stype = this->getStype();
    if (stype == "OPEN") {
        return 1.0 / Physics::epsilon_0; 
    } 

    // Standard coupling constant (from before)
    return 1.0 / (4.0 * Physics::pi * Physics::epsilon_0);
}  

/*
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
