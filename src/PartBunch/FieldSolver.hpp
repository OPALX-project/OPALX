#ifndef OPAL_FIELD_SOLVER_H
#define OPAL_FIELD_SOLVER_H

#include <memory>
#include <string>
#include <utility>
#include "BCHandler.hpp"
#include "Manager/BaseManager.h"
#include "Manager/FieldSolverBase.h"
#include "SpaceCharge/Ippl/IpplPoissonAdapter.h"
#include "SpaceCharge/Pic/PicWorkspace.h"
#include "Utilities/OpalException.h"

namespace opalx::spacecharge {
    class SelfFieldDiagnostics;
}

// Define the FieldSolver class
template <typename T, unsigned Dim>
class FieldSolver : public ippl::FieldSolverBase<T, Dim> {
private:
    using Workspace_t = opalx::spacecharge::PicWorkspace<T, Dim>;

    std::shared_ptr<Workspace_t> workspace_m;
    std::string greensFunction_m;
    opalx::spacecharge::PoissonBackendKind backendKind_m;
    opalx::spacecharge::GreenFunctionKind greenFunctionKind_m;
    std::unique_ptr<opalx::spacecharge::IpplPoissonAdapter> backend_m;

    using BCHandler_t = BCHandler<Dim>;
    std::shared_ptr<BCHandler_t> bcHandler_m;

public:
    FieldSolver(
            std::string solver, std::shared_ptr<Workspace_t> workspace,
            std::shared_ptr<BCHandler_t> bcHandler, std::string greensFunction = "STANDARD")
        : ippl::FieldSolverBase<T, Dim>(solver),
          workspace_m(std::move(workspace)),
          greensFunction_m(std::move(greensFunction)),
          backendKind_m(backendKindFromName(solver)),
          greenFunctionKind_m(greenFunctionKindFromName(greensFunction_m)),
          bcHandler_m(bcHandler) {
        if (!workspace_m) {
            throw OpalException("FieldSolver::FieldSolver", "PIC workspace is null.");
        }
        setPotentialBCs();
    }

    ~FieldSolver() override = default;

    void dumpScalField(std::string what, std::size_t solveIndex);
    void dumpVectField(std::string what, std::size_t solveIndex);

    Field_t<Dim>* getRho() { return &workspace_m->chargeDensity(); }
    const Field_t<Dim>* getRho() const { return &workspace_m->chargeDensity(); }
    VField_t<T, Dim>* getE() { return &workspace_m->electricField(); }
    const VField_t<T, Dim>* getE() const { return &workspace_m->electricField(); }
    Field<T, Dim>* getPhi() { return &workspace_m->potential(); }
    const Field<T, Dim>* getPhi() const { return &workspace_m->potential(); }
    [[nodiscard]] Workspace_t& getWorkspace() { return *workspace_m; }
    [[nodiscard]] const Workspace_t& getWorkspace() const { return *workspace_m; }
    [[nodiscard]] opalx::spacecharge::IpplPoissonAdapter& getBackendAdapter() {
        if (backend_m == nullptr) {
            throw OpalException(
                    "FieldSolver::getBackendAdapter",
                    "The IPPL backend must be initialized before it is borrowed.");
        }
        return *backend_m;
    }

    std::shared_ptr<BCHandler_t> getBCHandler() const { return bcHandler_m; }
    void setBCHandler(std::shared_ptr<BCHandler_t> bcHandler) {
        bcHandler_m = bcHandler;
        setPotentialBCs();
    }

    /**
     * @brief Return the OPAL `GREENSF` selection used for open-boundary FFT solves.
     *
     * The value follows OPAL input syntax (`STANDARD` or `INTEGRATED`). It is
     * mapped to IPPL's `greens_function` parameter when `TYPE=OPEN` initializes
     * the `FFTOpenPoissonSolver`.
     */
    std::string getGreensFunction() const { return greensFunction_m; }

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

    [[nodiscard]] const opalx::spacecharge::IpplPoissonCapabilities& getBackendCapabilities() const;

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

    void runSolver() override { runSolver(false); }

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

    /** @brief Execute and diagnose one runtime backend solve. */
    void runSolver(
            bool force_skip_field_dump, opalx::spacecharge::SelfFieldDiagnostics& diagnostics);

    /**
     * @brief Execute one typed backend request and retain common diagnostics and dump handling.
     *
     * Correction plans use this host-only entry point to select the standard or shifted Green
     * function without bypassing the shared backend solve counters, timers, or debug dumps.
     */
    void runSolver(
            const opalx::spacecharge::IpplPoissonSolveRequest& request, bool force_skip_field_dump,
            opalx::spacecharge::SelfFieldDiagnostics& diagnostics);

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
    void runShiftedOpenSolver(
            const ippl::Vector<double, Dim>& shift,
            opalx::spacecharge::SelfFieldDiagnostics& diagnostics);

private:
    static opalx::spacecharge::PoissonBackendKind backendKindFromName(const std::string& solver);
    static opalx::spacecharge::GreenFunctionKind greenFunctionKindFromName(
            const std::string& greenFunction);

    [[nodiscard]] opalx::spacecharge::IpplPoissonFields fields();

    void runSolverImpl(
            const opalx::spacecharge::IpplPoissonSolveRequest& request, bool force_skip_field_dump,
            opalx::spacecharge::SelfFieldDiagnostics* diagnostics);
};

#endif
