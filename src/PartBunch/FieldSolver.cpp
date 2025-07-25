#include "PartBunch/FieldSolver.hpp"

#include <iomanip>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "Utilities/Util.h"
#include "AbstractObjects/OpalData.h"

extern Inform* gmsg;

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

    auto fieldV      = field->getView();
    auto field_hostV = field->getHostMirror();
    Kokkos::deep_copy(field_hostV, fieldV);     

    boost::filesystem::path file(dirname);
    boost::format filename("%1%-%2%-%|3$06|.dat");
    std::string basename = OpalData::getInstance()->getInputBasename();
    filename % basename % (what + std::string("_") + type) % call_counter_m;
    file /= filename.str();
    std::ofstream fout(file.string(), std::ios::out);

    fout << std::setprecision(9);

    fout << "# " << Util::toUpper(what) << " " << type << " data on grid" << std::endl
         << "# origin= " << std::fixed << origin << " h= " << std::fixed << spacing << std::endl 
         << std::setw(5)  << "i"
         << std::setw(5)  << "j"
         << std::setw(5)  << "k"
         << std::setw(17) << "x [m]"
         << std::setw(17) << "y [m]"
         << std::setw(17) << "z [m]";

    fout << std::setw(10) << what << "x [" << unit << "]"
         << std::setw(10) << what << "y [" << unit << "]"
         << std::setw(10) << what << "z [" << unit << "]";

    fout << std::endl;

    for (int i = localIdx[0].first() + 1; i <= localIdx[0].last() +1 ; i++) {
            for (int j = localIdx[1].first() + 1 ; j <= localIdx[1].last() +1  ; j++) {
                for (int k = localIdx[2].first() + 1 ; k <= localIdx[2].last() +1 ; k++) {
                    
                    // define the physical points (cell-centered)
                    double x = i * spacing[0] + origin[0];        
                    double y = j * spacing[1] + origin[1];        
                    double z = k * spacing[2] + origin[2];     
                
                    fout << std::setw(5) << i-1 
                         << std::setw(5) << j-1
                         << std::setw(5) << k-1
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

    //  what == phi | rho

    Inform m("FS::dumpScalField() ");

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
    bool isVectorField = false;
    
    if (Util::toUpper(what) == "RHO") {
        type = "scalar";
        unit = "Cb/m^3";
    } else if (Util::toUpper(what) == "PHI") {
        type = "scalar";
        unit = "V";
    }

    
    Field_t<3>* field = this->getRho();   // both rho and phi are in the same variable (in place computation)
    
    auto localIdx = field->getOwned();
    auto mesh_mp  = &(field->get_mesh());
    auto spacing  = mesh_mp->getMeshSpacing();
    auto origin   = mesh_mp->getOrigin();

    auto fieldV      = field->getView();
    auto field_hostV = field->getHostMirror();
    Kokkos::deep_copy(field_hostV, fieldV);     

    boost::filesystem::path file(dirname);
    boost::format filename("%1%-%2%-%|3$06|.dat");
    std::string basename = OpalData::getInstance()->getInputBasename();
    filename % basename % (what + std::string("_") + type) % call_counter_m;
    file /= filename.str();
    std::ofstream fout(file.string(), std::ios::out);

    fout << std::setprecision(9);

    fout << "# " << Util::toUpper(what) << " " << type << " data on grid" << std::endl
         << "# origin= " << std::fixed << origin << " h= " << std::fixed << spacing << std::endl 
         << std::setw(5)  << "i"
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
        for (int i = localIdx[0].first(); i <= localIdx[0].last(); i++) {
            for (int j = localIdx[1].first(); j <= localIdx[1].last(); j++) {
                for (int k = localIdx[2].first(); k <= localIdx[2].last(); k++) {
                    
                    // define the physical points (cell-centered)
                    double x = i * spacing[0] + origin[0];        
                    double y = j * spacing[1] + origin[1];        
                    double z = k * spacing[2] + origin[2];     
                
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
    else {
        for (int i = localIdx[0].first() + 1; i <= localIdx[0].last() +1 ; i++) {
            for (int j = localIdx[1].first() + 1 ; j <= localIdx[1].last() +1  ; j++) {
                for (int k = localIdx[2].first() + 1 ; k <= localIdx[2].last() +1 ; k++) {
                    
                    // define the physical points (cell-centered)
                    double x = i * spacing[0] + origin[0];        
                    double y = j * spacing[1] + origin[1];        
                    double z = k * spacing[2] + origin[2];     
                
                    fout << std::setw(5) << i-1 
                         << std::setw(5) << j-1
                         << std::setw(5) << k-1
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

template <typename T, unsigned Dim>
template <typename Solver>
void FieldSolver<T, Dim>::initSolverWithParams(const ippl::ParameterList& sp) {

    this->getSolver().template emplace<Solver>();

    if (std::holds_alternative<Solver>(this->getSolver())) {
        *gmsg << "* Solver was correctly emplaced." << endl;
    } else {
        throw std::runtime_error("Failed to emplace Solver in std::variant.");
    }

    Solver& solver = std::get<Solver>(this->getSolver());

    solver.mergeParameters(sp);

    solver.setRhs(*rho_m);

        if constexpr (std::is_same_v<Solver, CGSolver_t<T, Dim>>) {
            // The CG solver computes the potential directly and
            // uses this to get the electric field
            solver.setLhs(*phi_m);
            solver.setGradient(*E_m);
        }
        else{
            // The periodic Poisson solver, Open boundaries solver,
            // and the P3M solver compute the electric field directly
            solver.setLhs(*E_m);
            if constexpr (std::is_same_v<Solver, OpenSolver_t<T, Dim>>){
                solver.setGradFD();
            }
        }
        call_counter_m = 0;
}

template <>
void FieldSolver<double,3>::initOpenSolver() {
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
void FieldSolver<double,3>::initSolver() {
    Inform m;

    if (this->getStype() == "FFT") {
        initOpenSolver();
    } else if (this->getStype() == "FFTPERIODIC") {
        initFFTSolver();
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
void FieldSolver<double,3>::runSolver() {
    constexpr int Dim = 3;

    if (this->getStype() == "CG") {
            CGSolver_t<double, 3>& solver = std::get<CGSolver_t<double, 3>>(this->getSolver());
            solver.solve();

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
            ippl::Comm->barrier();
        } else if (this->getStype() == "FFT") {
            if constexpr (Dim == 2 || Dim == 3) {
#ifdef OPALX_FIELD_DEBUG
                this->dumpScalField("rho");
                call_counter_m++;
#endif
                std::get<OpenSolver_t<double, 3>>(this->getSolver()).solve();
#ifdef OPALX_FIELD_DEBUG
                this->dumpScalField("phi");
                this->dumpVectField("ef");
                call_counter_m++;
#endif
            }
        } else if (this->getStype() == "P3M") {
            if constexpr (Dim == 3) {
                std::get<FFTTruncatedGreenSolver_t<double, 3>>(this->getSolver()).solve();
            }
        } else if (this->getStype() == "FFTPERIODIC") {
            if constexpr (Dim == 3) {
#ifdef OPALX_FIELD_DEBUG
                this->dumpScalField("rho");
                call_counter_m++;
#endif
                std::get<FFTSolver_t<T, Dim>>(this->getSolver()).solve();
#ifdef OPALX_FIELD_DEBUG
                this->dumpScalField("phi");
                this->dumpVectField("ef");
                call_counter_m++;
#endif
            }
    } else if (this->getStype() == "NONE") {
        std::get<NullSolver_t<T, Dim>>(this->getSolver()).solve();
    } else {
        throw std::runtime_error("Unknown solver type");
    }
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

/*
template <>
void FieldSolver<double,3>::initCGSolver() {
    ippl::ParameterList sp;
    sp.add("output_type", CGSolver_t<double, 3>::GRAD);
    // Increase tolerance in the 1D case
    sp.add("tolerance", 1e-10);
    
    initSolverWithParams<CGSolver_t<double, 3>>(sp);
}

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

template class FieldSolver<double, 3>;
template void FieldSolver<double, 3>::initSolverWithParams<CGSolver_t<double, 3>>(const ippl::ParameterList&);
template void FieldSolver<double, 3>::initSolverWithParams<OpenSolver_t<double, 3>>(const ippl::ParameterList&);
template void FieldSolver<double, 3>::initSolverWithParams<FFTSolver_t<double, 3>>(const ippl::ParameterList&);
