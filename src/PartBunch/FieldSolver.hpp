#ifndef OPAL_FIELD_SOLVER_H
#define OPAL_FIELD_SOLVER_H

#include <memory>
#include "BCHandler.hpp"
#include "Manager/BaseManager.h"
#include "Manager/FieldSolverBase.h"

// Define the FieldSolver class
template <typename T, unsigned Dim>
class FieldSolver : public ippl::FieldSolverBase<T, Dim> {
private:
    Field_t<Dim>* rho_m;
    VField_t<T, Dim>* E_m;
    Field_t<Dim>* phi_m;

    using BCHandler_t = BCHandler<Dim>;
    std::shared_ptr<BCHandler_t> bcHandler_m;

    /// Counts number of times the solver has been called
    size_t call_counter_m;

public:
    FieldSolver(
            std::string solver, Field_t<Dim>* rho, VField_t<T, Dim>* E, Field_t<Dim>* phi,
            std::shared_ptr<BCHandler_t> bcHandler)
        : ippl::FieldSolverBase<T, Dim>(solver),
          rho_m(rho),
          E_m(E),
          phi_m(phi),
          bcHandler_m(bcHandler),
          call_counter_m(0) {
        setPotentialBCs();
    }

    ~FieldSolver() override = default;

    void dumpScalField(std::string what);
    void dumpVectField(std::string what);

    Field_t<Dim>* getRho() { return rho_m; }
    void setRho(Field_t<Dim>* rho) { rho_m = rho; }

    VField_t<T, Dim>* getE() const { return E_m; }
    void setE(VField_t<T, Dim>* E) { E_m = E; }

    Field<T, Dim>* getPhi() const { return phi_m; }
    void setPhi(Field<T, Dim>* phi) { phi_m = phi; }

    std::shared_ptr<BCHandler_t> getBCHandler() const { return bcHandler_m; }
    void setBCHandler(std::shared_ptr<BCHandler_t> bcHandler) {
        bcHandler_m = bcHandler;
        setPotentialBCs();
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

    void initSolver() override;

    /**
     * @brief Set boundary conditions for the electrostatic potential field.
     *
     * Converts the boundary-condition specification provided by the BC handler
     * into the IPPL boundary-condition format for Field_t<Dim> and applies the
     * resulting conditions to the internal potential field (phi_m) by calling
     * its setFieldBC method.
     *
     * @throws OpalException if the BC handler is not set or invalid.
     */
    void setPotentialBCs();

    bool hasValidBCHandler() const { return (bcHandler_m != nullptr); }

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

    size_t getCallCounter() { return call_counter_m; }

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

    /**
     * @brief Run an Open-solver solve with a shifted free-space Green's function.
     *
     * Installs a translated Green's kernel @f$G(r) = -1/(4\pi|r - \texttt{shift}|)@f$
     * on the underlying IPPL @c FFTOpenPoissonSolver via
     * @c shiftedGreensFunction(shift), runs @c solve(), then restores the
     * standard Green's function (@c greensFunction()) so subsequent solves in
     * later bins are not affected.
     *
     * The manual restore guards against two adjacent bins whose stretched mesh
     * spacings happen to coincide: the mesh-change detector inside the IPPL
     * @c solve() would then @em not recompute the kernel and the next primary
     * solve would silently reuse the shifted one. With the current binning
     * algorithm this collision is not supposed to happen, so the extra FFT
     * per shifted pass is defensive and can be removed once the invariant is
     * guaranteed by the binner.
     *
     * @param shift  Per-axis translation of the Green's function in physical units.
     *
     * @throws OpalException if the configured solver is not @c "OPEN".
     */
    void runShiftedOpenSolver(const ippl::Vector<double, Dim>& shift);

    template <typename Solver>
    void initSolverWithParams(const ippl::ParameterList& sp);

    void initNullSolver();

    void initFFTSolver();

    void initCGSolver() {}

    void initP3MSolver() {}
};

// Explicit specialization declaration
template <>
void FieldSolver<double, 3>::initNullSolver();

#endif
