#ifndef OPAL_FIELD_SOLVER_H
#define OPAL_FIELD_SOLVER_H

#include <memory>
#include "Manager/BaseManager.h"
#include "Manager/FieldSolverBase.h"

// Define the FieldSolver class
template <typename T, unsigned Dim>
class FieldSolver : public ippl::FieldSolverBase<T, Dim> {
private:
    Field_t<Dim>* rho_m;
    VField_t<T, Dim>* E_m;
    Field_t<Dim>* phi_m;
    unsigned int call_counter_m;
public:
    FieldSolver(std::string solver,
                Field_t<Dim>* rho,
                VField_t<T, Dim>* E,
                Field<T, Dim>* phi)
        : ippl::FieldSolverBase<T, Dim>(solver),
          rho_m(rho), E_m(E), phi_m(phi), call_counter_m(0) {

        setPotentialBCs();
    }

    ~FieldSolver() {
    }

    void dumpScalField(std::string what);
    void dumpVectField(std::string what);

    Field_t<Dim>* getRho() {
        return rho_m;
    }
    void setRho(Field_t<Dim>* rho) {
        rho_m = rho;
    }

    VField_t<T, Dim>* getE() const {
        return E_m;
    }
    void setE(VField_t<T, Dim>* E) {
        E_m = E;
    }

    Field<T, Dim>* getPhi() const {
        return phi_m;
    }
    void setPhi(Field<T, Dim>* phi) {
        phi_m = phi;
    }

    void initOpenSolver();

    void initSolver() override ;

    void setPotentialBCs();

    void runSolver() override;

    template <typename Solver>
    void initSolverWithParams(const ippl::ParameterList& sp);

    void initNullSolver();
    
    void initFFTSolver() {
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
    //initSolverWithParams<OpenSolver_t<double, 3>>(sp);
    }
    
    void initCGSolver() { }

    void initP3MSolver() { }

};

// Explicit specialization declaration
template<>
void FieldSolver<double, 3>::initNullSolver();

#endif
