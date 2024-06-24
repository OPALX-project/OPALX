#include "PartBunch/FieldSolver.hpp"

extern Inform* gmsg;

template <>
void FieldSolver<double,3>::initOpenSolver() {
    //      if constexpr (Dim == 3) {
        ippl::ParameterList sp;
        sp.add("output_type", OpenSolver_t<double, 3>::GRAD);
        sp.add("use_heffte_defaults", false);
        sp.add("use_pencils", true);
        sp.add("use_reorder", false);
        sp.add("use_gpu_aware", true);
        sp.add("comm", ippl::p2p_pl);
        sp.add("r2c_direction", 0);
        sp.add("algorithm", OpenSolver_t<double, 3>::HOCKNEY);
        initSolverWithParams<OpenSolver_t<double, 3>>(sp);
        // } else {
        //throw std::runtime_error("Unsupported dimensionality for OPEN solver");
        //}
}

template <>
void FieldSolver<double,3>::initSolver() {
    Inform m;
    if (this->getStype() == "FFT") {
        initOpenSolver();
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
                std::get<OpenSolver_t<double, 3>>(this->getSolver()).solve();
            }
        } else if (this->getStype() == "P3M") {
            if constexpr (Dim == 3) {
                std::get<P3MSolver_t<double, 3>>(this->getSolver()).solve();
            }
        } else if (this->getStype() == "FFTOPEN") {
            if constexpr (Dim == 3) {
                std::get<FFTSolver_t<double, 3>>(this->getSolver()).solve();
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