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

    /// Counts number of times the solver has been called
    size_t call_counter_m;
public:
    FieldSolver(std::string solver,
                Field_t<Dim>* rho,
                VField_t<T, Dim>* E,
                Field<T, Dim>* phi)
        : ippl::FieldSolverBase<T, Dim>(solver),
          rho_m(rho), E_m(E), phi_m(phi), call_counter_m(0) {

        setPotentialBCs();
    }

    ~FieldSolver() override = default;

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

    /**
     * @brief Get the solver's coupling constant.
     *
     * Returns the scalar coupling constant used by the field solver to scale
     * interactions between particles and the field. This value is applied 
     * during `ParBunch::computeSpaceCharge`. Its physical meaning and units 
     * potentially depend on the specific solver type used.
     *
     * @return The coupling constant of type T (usually double).
     */
    T getCouplingConstant() const;

    void initOpenSolver();

    void initSolver() override ;

    void setPotentialBCs();

    void runSolver() override {
        // The default runSolver should always count towards the call counter!
        runSolver(false);
    }

    /**
     * @brief Reset the solver call counter to zero.
     *
     * Sets the internal call counter (`call_counter_m`) back to 0 so that
     * subsequent calls will be counted from a clean state.
     * 
     * @note This function is necessary to exclude potential solver warm-up
     * calls from being counted towards output or logging that depends on the
     * number of solver executions. 
     */
    void resetCallCounter() { call_counter_m = 0; }

    /**
     * @brief Execute the field solver for the current simulation state.
     *
     * Performs a single solve cycle using the solver's current configuration,
     * boundary conditions and particle/mesh data. The solver updates the
     * internal field representations.
     *
     * @param force_skip_field_dump
     *     If true, suppress any field-dump output that would otherwise be
     *     produced by this call. If false, field output behavior follows the
     *     configured/normal schedule. 
     * 
     * @note This second implementation is necessary since the pure 
     * `runSolver()` routine is defined in the base class as not taking any 
     * arguments. 
     */
    void runSolver(bool force_skip_field_dump);

    template <typename Solver>
    void initSolverWithParams(const ippl::ParameterList& sp);

    void initNullSolver();
    
    void initFFTSolver();
    
    void initCGSolver() { }

    void initP3MSolver() { }

};

// Explicit specialization declaration
template<>
void FieldSolver<double, 3>::initNullSolver();

#endif
