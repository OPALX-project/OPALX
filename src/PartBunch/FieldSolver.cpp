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

    Inform m("FS::dumpScalField() ");

    //    std::variant<Field_t<3>*, VField_t<double, 3>* > field;

    if (ippl::Comm->size() > 1 || call_counter<2) {
        return;
    }

    
    m << "*** START DUMPING VECTOR FIELD ***" << endl;

/*
    constexpr bool isVectorField = std::is_same<VField_t, FieldType>::value;
    std::string type = (isVectorField) ? "field" : "scalar";
    
    
*/

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
 *   '******':   step padded with zeros to 6 digits
 */
    
    std::string dirname = "";

    std::string type;
    std::string unit;
    bool isVectorField;

    if (Util::toUpper(what) == "EF") {
        type = "vector";
        unit = "";
        isVectorField = true;
        //    field = this->getE();
    }


    boost::filesystem::path file(dirname);
    boost::format filename("%1%-%2%-%|3$06|.dat");
    std::string basename = OpalData::getInstance()->getInputBasename();
    filename % basename % (what + std::string("_") + type) % call_counter;
    file /= filename.str();
    m << "*** FILE NAME " + file.string() << endl;
    std::ofstream fout(file.string(), std::ios::out);
    fout.precision(9);

    fout << "# " << Util::toUpper(what) << " " << type << " data on grid" << std::endl
         << "#"
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

    fout.close();
    m << "*** FINISHED DUMPING " + Util::toUpper(what) + " FIELD ***" << endl;
}

template <>
void FieldSolver<double,3>::dumpScalField(std::string what) {
    /*
      what == phi | rho
     */

    Inform m("FS::dumpScalField() ");

    //    std::variant<Field_t<3>*, VField_t<double, 3>* > field;

    if (ippl::Comm->size() > 1 || call_counter<2) {
        return;
    }

    m << "*** START DUMPING SCALAR FIELD ***" << endl;

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
 *   '******':   step padded with zeros to 6 digits
 */
    
    int step = 0;
    std::string dirname = "";

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
    filename % basename % (what + std::string("_") + type) % step;
    file /= filename.str();
    m << "*** FILE NAME " + file.string() << endl;
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

    for (int i = localIdx[0].first() +1; i <= localIdx[0].last() +1; i++) {
        for (int j = localIdx[1].first() +1; j <= localIdx[1].last() +1; j++) {
            for (int k = localIdx[2].first() +1; k <= localIdx[2].last() +1; k++) {

                // define the physical points (cell-centered)
                double x = i * spacing[0] + origin[0];        
                double y = j * spacing[1] + origin[1];        
                double z = k * spacing[2] + origin[2];     
                
                fout << std::setw(5) << i + 1
                     << std::setw(5) << j + 1
                     << std::setw(5) << k + 1
                     << std::setw(17) << x
                     << std::setw(17) << y
                     << std::setw(17) << z
                     << std::scientific << "\t" << field_hostV(i,j,k)                             
                     << std::endl;
            }
        }
    }
    
    fout.close();
    m << "*** FINISHED DUMPING " + Util::toUpper(what) + " FIELD ***" << endl;
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

    initOpenSolver();
    /*
    Inform m;
    if (this->getStype() == "FFT") {
    
    } else if (this->getStype() == "CG") {
        initCGSolver();
    } else if (this->getStype() == "P3M") {
        initP3MSolver();
    } else if (this->getStype() == "FFTOPEN") {
        initFFTSolver();
    } else if (this->getStype() == "NONE") {
        initNullSolver();
    }
    else {
        m << "No solver matches the argument: " << this->getStype() << endl;
    }
    */
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
                this->dumpScalField("rho");
                call_counter++;
                std::get<OpenSolver_t<double, 3>>(this->getSolver()).solve();
                this->dumpScalField("phi");
                call_counter++;
            }
        } else if (this->getStype() == "P3M") {
            if constexpr (Dim == 3) {
                std::get<P3MSolver_t<double, 3>>(this->getSolver()).solve();
            }
        } else if (this->getStype() == "FFTOPEN") {
            if constexpr (Dim == 3) {
                this->dumpScalField("rho");
                call_counter++;
                std::get<OpenSolver_t<double, 3>>(this->getSolver()).solve();
                this->dumpScalField("phi");
                call_counter++;
            }
        } else {
            throw std::runtime_error("Unknown solver type");
        }

}

/*
template<>
void FieldSolver<double,3>::initNullSolver() {
    //        if constexpr (Dim == 2 || Dim == 3) {
    //ippl::ParameterList sp;
    throw std::runtime_error("Not implemented Null solver");
    //  } else {
    //throw std::runtime_error("Unsupported dimensionality for Null solver");
    //}
}
*/

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
