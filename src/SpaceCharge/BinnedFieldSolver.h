#ifndef OPALX_BINNED_FIELD_SOLVER_H
#define OPALX_BINNED_FIELD_SOLVER_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "SpaceCharge/FieldSolver.hpp"
#include "SpaceCharge/BinContext.hpp"
#include "SpaceCharge/Corrections/BoundaryCorrection.hpp"
#include "SpaceCharge/Corrections/ImageChargeCorrection.hpp"
#include "SpaceCharge/Corrections/ShiftedGreensCorrection.hpp"
#include "SpaceCharge/Drivers/BinnedLorentzSelfFieldDriver.hpp"
#include "SpaceCharge/Drivers/MonolithicSelfFieldDriver.hpp"
#include "SpaceCharge/Drivers/SelfFieldDriver.hpp"
#include "SpaceCharge/FieldAccumulator.hpp"
#include "SpaceCharge/FieldSolverConfig.hpp"
#include "SpaceCharge/ImageChargeScatterController.h"
#include "PartBunch/PartBunch.h"
#include "Utilities/OpalException.h"

/**
 * @brief Field solver wrapper that implements the full binned self-field algorithm.
 *
 * @tparam T   Particle numeric type (typically `double`).
 * @tparam Dim Spatial dimensionality (currently only `Dim == 3` is supported).
 *
 * Design overview:
 * - The solver owns no particle data; it borrows a `PartBunch` by reference.
 * - Runtime selection:
 *   - If adaptive bins are available (`bunch.hasBinning()` / `bunch.getBins()`),
 *     the per-bin algorithm is executed.
 *   - Otherwise, the legacy monolithic scatter/solve/gather path is used.
 * - The solver currently supports only `ChargeQ -> rho` scattering and `ElectricFieldE`
 *   gathering, but gathers both E and B fields.
 * - Physics details:
 *   - See https://github.com/aliemen/HS24-masters-thesis for details.
 *   - After bins are calculated, it solves electrostatics per bin in a quasi-static approximation.
 *   - Fields per bin are then transformed to the lab frame and accumulated into the temporary
 *     fields, this also produces the magnetic field contributions.
 *   - Finally, the accumulated fields are gathered back to the particles.
 *   - This procedure approximates full Maxwell's equations for the self-fields.
 *   - Without a bins object, it falls back to the legacy electrostatic approximation.
 */
template <typename T, unsigned Dim>
class BinnedFieldSolver : public FieldSolver<T, Dim> {
    static_assert(Dim == 3, "BinnedFieldSolver currently supports Dim == 3 only.");

    template <typename, unsigned>
    friend class BinnedLorentzSelfFieldDriver;
    template <typename, unsigned>
    friend class MonolithicSelfFieldDriver;

public:
    using PartBunch_t    = PartBunch<T, Dim>;
    using ParticleCtr_t  = typename PartBunch_t::ParticleContainer_t;
    using AdaptBins_t    = typename PartBunch_t::AdaptBins_t;
    using BCHandler_t    = BCHandler<Dim>;
    using bin_index_type = typename AdaptBins_t::bin_index_type;
    using size_type      = typename AdaptBins_t::size_type;

    using particle_position_type = typename PartBunch_t::Base::particle_position_type;

    /**
     * @brief Which particle attribute to scatter from to build the mesh charge density `rho`.
     *
     * Currently only `ChargeQ` is implemented.
     */
    enum class ScatterAttribute { ChargeQ };

    /**
     * @brief Which particle attribute to gather the accumulated electric field into.
     *
     * Currently only `ElectricFieldE` is implemented.
     */
    enum class GatherAttribute { ElectricFieldE };

    /**
     * @brief Construct a binned/legacy-compatible solver.
     *
     * @param config     Normalized field-solver command and backend capability state.
     * @param rho        Pointer to the mesh charge-density field storage.
     * @param E          Pointer to the mesh electric-field storage (solver output).
     * @param phi        Pointer to the potential-field storage (solver internal use).
     * @param bcHandler  Shared pointer to the boundary-condition handler.
     */
    BinnedFieldSolver(
            const opalx::FieldSolverConfig<Dim>& config, Field_t<Dim>* rho,
            VField_t<T, Dim>* E, Field_t<Dim>* phi, std::shared_ptr<BCHandler_t> bcHandler);

    /**
     * @brief Compute space-charge self-fields for the given particle bunch.
     *
     * If the bunch provides adaptive binning (`bunch.getBins()`), the solver executes
     * the per-bin algorithm:
     * `scatter rho corrections -> solve -> Lorentz scaling -> accumulate -> gather`.
     * Otherwise, it executes the legacy monolithic algorithm:
     * `scatter all particles -> solve once -> gather directly`.
     *
     * @param bunch Particle bunch to update. Ownership remains with the caller.
     *
     * @throws OpalException If required internal data (particle container / temp E field)
     *                        is missing, or if unsupported scatter/gather modes are selected.
     */
    void computeSelfFields(PartBunch_t& bunch);

    /**
     * @brief Name of the self-field driver selected for the bunch state.
     */
    const char* selectedSelfFieldDriverName(const PartBunch_t& bunch) const;

    /**
     * @brief Normalized solver configuration used by the runtime.
     */
    const opalx::FieldSolverConfig<Dim>& getConfig() const { return config_m; }

    /**
     * @brief Set particle scatter attribute (extensible; default is `ChargeQ`).
     * @param attr Attribute to scatter from.
     */
    void setScatterAttribute(const ScatterAttribute attr);

    /**
     * @brief Set particle gather attribute (extensible; default is `ElectricFieldE`).
     * @param attr Attribute to gather into.
     */
    void setGatherAttribute(const GatherAttribute attr);

    /// @brief Configure optional image-charge scatter pass.
    void setImageChargeConfiguration(bool enabled, double zPlane);
    bool isImageChargeEnabled() const { return imageCorrection_m.isEnabled(); }
    double getImageChargePlaneZ() const { return imageCorrection_m.planeZ(); }

    /// @brief Configure the shifted Green's function Dirichlet correction (alternative to
    /// image charges). Mutually exclusive with @c setImageChargeConfiguration(true, ...).
    /// Requires a backend with shifted-Greens support.
    void setShiftedGreensConfiguration(bool enabled, double zPlane);
    bool isShiftedGreensEnabled() const { return shiftedGreensCorrection_m.isEnabled(); }
    double getShiftedGreensPlaneZ() const { return shiftedGreensCorrection_m.planeZ(); }

    /// @brief Set the maximum number of timesteps for which image charges are active (0 =
    /// unlimited).
    void setZerofaceMaxSteps(int maxSteps);
    int getZerofaceMaxSteps() const { return zerofaceMaxSteps_m; }

    /// @brief Check whether the explicit image-charge pass should run for a given timestep.
    bool isImageChargeActiveForStep(size_t step) const;

    /// @brief Check whether the shifted Green's function correction should run for a given
    /// timestep. Reuses the same step budget (@c zerofaceMaxSteps_m) as the image-charge path.
    bool isShiftedGreensActiveForStep(size_t step) const;

    /// @brief Configure dump frequency for dirichlet-plane diagnostics (`0` disables dumps).
    void setZeroFacePlaneDumpFrequency(int frequency);
    int getZeroFacePlaneDumpFrequency() const { return zeroFacePlaneDumpFrequency_m; }

    struct BinKinematics {
        Vector_t<double, Dim> pmean = Vector_t<double, Dim>(0.0);
        double gammaBin             = 1.0;
    };

private:
    ScatterAttribute scatterAttribute_m;
    GatherAttribute gatherAttribute_m;
    opalx::FieldSolverConfig<Dim> config_m;
    int zeroFacePlaneDumpFrequency_m = 0;
    int zerofaceMaxSteps_m           = 0;
    ImageChargeScatterController<T, Dim> imageScatterController_m;
    ImageChargeCorrection<T, Dim> imageCorrection_m;
    ShiftedGreensCorrection<T, Dim> shiftedGreensCorrection_m;
    bool warnedPlaneDumpParallelUnsupported_m = false;

    FieldAccumulator<T, Dim> fieldAccumulator_m;
    std::unique_ptr<MonolithicSelfFieldDriver<T, Dim>> monolithicDriver_m;
    std::unique_ptr<BinnedLorentzSelfFieldDriver<T, Dim>> binnedDriver_m;

    /**
     * @brief Row entry for the level-3 bin statistics table.
     */
    struct BinStatsRow {
        long long binNumber;            //!< Merged bin index (or `-1` for legacy mode).
        unsigned long long nParticles;  //!< Number of particles in the (merged) bin.
        double gammaBin;                //!< Global average gamma for the (merged) bin.
    };

    /**
     * @brief Print the bin statistics table at level 3.
     *
     * The output includes columns for bin index, particle count, and `gammaBin`.
     * In binned mode, rows correspond to each merged bin. In legacy mode, a single
     * row with `binNumber = -1` is printed.
     *
     * @param tableName    Logical table name (used in the header).
     * @param nBinsOrZero  Number of bins (for header metadata). Use `0` for legacy mode.
     * @param rows         Table rows to print.
     */
    void printBinStatsTable(
            const std::string& binningCmdName, const std::vector<BinStatsRow>& rows);

private:
    /**
     * @brief Dump and report potential values interpolated on the image-charge plane.
     *
     * If diagnostics are disabled or conditions are not met for the current step,
     * this function returns without side effects.
     *
     * @param bunch   Particle bunch context for time/step and DataSink access.
     * @param solveTag Label used in output file naming (`legacy`, `binned`, ...).
     */
    void dumpDirichletPlaneDiagnosticsIfRequested(PartBunch_t& bunch, const std::string& solveTag);

    /**
     * @brief Build and prepare adaptive bins for the current step.
     *
     * This performs a full rebin to the maximum bin count, sorts the particle container by bin,
     * generates the adaptive histogram (merging), and restores bin configuration.
     *
     * @param bunch Bunch whose bins are updated.
     * @param bins  Adaptive bin structure owned/managed by the bunch.
     */
    void rebinAndPrepare(PartBunch_t& bunch, std::shared_ptr<AdaptBins_t> bins);

    const BoundaryCorrection<T, Dim>* configuredBoundaryCorrection() const;
    const BoundaryCorrection<T, Dim>* activeBoundaryCorrectionForStep(size_t step) const;
    SelfFieldDriver<T, Dim>& selectedSelfFieldDriver(const PartBunch_t& bunch);
    const SelfFieldDriver<T, Dim>& selectedSelfFieldDriver(const PartBunch_t& bunch) const;
    ImageScatterMode primaryScatterModeForStep(
            const BoundaryCorrection<T, Dim>* activeCorrection) const;
    void applyBoundaryCorrectionPass(
            PartBunch_t& bunch, std::shared_ptr<AdaptBins_t> bins,
            const BinContext<T, Dim>& context, const BoundaryCorrectionPass<T, Dim>& pass,
            bool& dumpedDirichletPlaneThisStep);

public:
    /**
     * @brief Compute per-bin global mean momentum and gamma.
     *
     * The function computes the global mean momentum vector `pmean` across all
     * particles in the merged bin and derives:
     * `gammaBin = mean(sqrt(1 + dot(p_i, p_i)))`.
     *
     * @param bunch        Bunch providing particle data.
     * @param bins         Bins providing the merged-bin iteration policy and indexing.
     * @param binIndex     Bin index in the merged-bin space.
     * @param nPartGlobal  Global particle count for this merged bin.
     *
     * @return Per-bin kinematics bundle (`pmean`, `gammaBin`).
     */
    BinKinematics computeGammaBinGlobal(
            PartBunch_t& bunch, std::shared_ptr<AdaptBins_t> bins, const bin_index_type binIndex,
            const size_type nPartGlobal) const;

    /**
     * @brief Build mesh `rho` for a specific merged bin and apply all corrections.
     *
     * Steps mirror the legacy ordering from the existing OPAL-X code paths and include:
     * - dt scaling for charge scattering,
     * - config-driven mesh normalization and optional background subtraction,
     * - Lorentz rest-frame scaling: `rho /= gammaBin`,
     * - scaling by the solver coupling constant.
     *
     * @param bunch        Bunch providing geometry and charge data.
     * @param bins         Adaptive bins providing bin iteration and hash indexing.
     * @param binIndex     Merged bin index.
     * @param nPartGlobal  Global number of particles in that merged bin.
     * @param gammaBin     Global average gamma for that merged bin.
     */
    void prepareRhoForBin(
            PartBunch_t& bunch, std::shared_ptr<AdaptBins_t> bins,
            const BinContext<T, Dim>& context,
            ImageScatterMode mode = ImageScatterMode::PrimaryAndImage);

    /**
     * @brief Gather the accumulated lab-frame fields from temporaries back to particles.
     */
    void gatherFromTempToParticles(PartBunch_t& bunch);
};

// Reduce compile-time churn: instantiate the only supported concrete solver in one TU.
extern template class BinnedFieldSolver<double, 3>;

#endif  // OPALX_BINNED_FIELD_SOLVER_H
