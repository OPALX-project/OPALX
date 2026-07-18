#ifndef OPALX_BINNED_FIELD_SOLVER_H
#define OPALX_BINNED_FIELD_SOLVER_H

#include <cstdint>
#include <iomanip>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "PartBunch/FieldSolver.hpp"
#include "PartBunch/PartBunch.h"
#include "SpaceCharge/Pic/CorrectionPlan.h"
#include "SpaceCharge/Pic/FieldComposer.h"
#include "SpaceCharge/Pic/IterationPlan.h"
#include "SpaceCharge/Pic/PicScatterGather.h"
#include "SpaceCharge/Pic/PicWorkspace.h"
#include "SpaceCharge/SelfFieldDiagnostics.h"
#include "Utilities/OpalException.h"

/**
 * @brief Field solver wrapper that implements the full binned self-field algorithm.
 *
 * @tparam T   Particle numeric type (typically `double`).
 * @tparam Dim Spatial dimensionality (currently only `Dim == 3` is supported).
 *
 * Design overview:
 * - The solver owns no particle data; it borrows a `PartBunch` by reference.
 * - Runtime selection comes from a solver-owned IterationPlan. BinningPlan emits one unit at a
 *   time; WholeBunchPlan emits a single monolithic unit.
 * - The solver currently supports only `ChargeQ -> rho` scattering and `ElectricFieldE`
 *   gathering, but gathers both E and B fields.
 * - Physics details:
 *   - See https://github.com/aliemen/HS24-masters-thesis for details.
 *   - After bins are calculated, it solves electrostatics per bin in a quasi-static approximation.
 *   - Fields per bin are then transformed to the lab frame and accumulated into the temporary
 *     fields, this also produces the magnetic field contributions.
 *   - Finally, the accumulated fields are gathered back to the particles.
 *   - This procedure approximates full Maxwell's equations for the self-fields.
 *   - A whole-bunch unit uses the legacy electrostatic approximation.
 */
template <typename T, unsigned Dim>
class BinnedFieldSolver : public FieldSolver<T, Dim> {
    static_assert(Dim == 3, "BinnedFieldSolver currently supports Dim == 3 only.");

public:
    using PartBunch_t          = PartBunch<T, Dim>;
    using ParticleCtr_t        = typename PartBunch_t::ParticleContainer_t;
    using Workspace_t          = opalx::spacecharge::PicWorkspace<T, Dim>;
    using ScatterGather_t      = opalx::spacecharge::PicScatterGather<T, Dim>;
    using FieldComposer_t      = opalx::spacecharge::FieldComposer<T, Dim>;
    using BCHandler_t          = BCHandler<Dim>;
    using IterationPlan_t      = opalx::spacecharge::IterationPlan<T, Dim>;
    using PreparedIteration_t  = opalx::spacecharge::PreparedIteration;
    using PreparedCorrection_t = opalx::spacecharge::PreparedCorrection<T, Dim>;
    using SolveUnit_t          = opalx::spacecharge::SolveUnit<T, Dim>;
    using SolvePass_t          = opalx::spacecharge::SolvePass<T, Dim>;

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
     * @param solver     Concrete solver name (e.g. `FFT`, `OPEN`, `CG`, `NONE`).
     * @param workspace  Shared persistent mesh, field, and scratch storage.
     * @param bcHandler  Shared pointer to the boundary-condition handler.
     * @param greensFunction  OPAL `GREENSF` selection (`STANDARD` or `INTEGRATED`) forwarded to
     *                        IPPL's open-boundary FFT solver.
     */
    BinnedFieldSolver(
            std::string solver, std::shared_ptr<Workspace_t> workspace,
            std::shared_ptr<BCHandler_t> bcHandler, std::string greensFunction = "STANDARD");

    /**
     * @brief Compute space-charge self-fields for the given particle bunch.
     *
     * If @p iterationPlan is a BinningPlan, the solver executes the per-bin algorithm:
     * `scatter rho corrections -> solve -> Lorentz scaling -> accumulate -> gather`.
     * Otherwise, it executes the legacy monolithic algorithm:
     * `scatter all particles -> solve once -> gather directly`.
     *
     * @param bunch Particle bunch to update. Ownership remains with the caller.
     *
     * @throws OpalException If required internal data (particle container / temp E field)
     *                        is missing, or if unsupported scatter/gather modes are selected.
     */
    void computeSelfFields(
            PartBunch_t& bunch, IterationPlan_t& iterationPlan, std::uint64_t particleGeneration,
            const PreparedCorrection_t& correction,
            opalx::spacecharge::BinConfigurationObserver* binObserver,
            opalx::spacecharge::SelfFieldDiagnostics& diagnostics);

    /** @brief Initialize the temporary `.stat` bin-count compatibility cache. */
    void configureIterationMetadata(bool binningConfigured, std::size_t maximumBinCount);

    /** @brief Preserve the legacy StatWriter rule for the reported `nBins` value. */
    [[nodiscard]] int legacyReportedBinCount() const;

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

private:
    ScatterAttribute scatterAttribute_m;
    GatherAttribute gatherAttribute_m;
    ScatterGather_t scatterGather_m;
    FieldComposer_t fieldComposer_m;
    bool warnedPlaneDumpParallelUnsupported_m = false;

    bool binningConfigured_m      = false;
    std::size_t currentBinCount_m = 1;
    std::size_t maximumBinCount_m = 0;

    /**
     * @brief Row entry for the level-3 bin statistics table.
     */
    struct BinStatsRow {
        long long binNumber;            //!< Merged bin index (or `-1` for legacy mode).
        unsigned long long nParticles;  //!< Number of particles in the (merged) bin.
        double gammaBin;                //!< Global average gamma for the (merged) bin.
    };

    std::vector<BinStatsRow> binStats_m;

    /**
     * @brief Print the bin statistics table at level 3.
     *
     * The output includes columns for bin index, particle count, and `gammaBin`.
     * In binned mode, rows correspond to each merged bin. In legacy mode, a single
     * row with `binNumber = -1` is printed.
     *
     * @param binningCmdName Logical binning-command name used in the Inform label.
     * @param rows           Table rows to print.
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
    void dumpDirichletPlaneDiagnosticsIfRequested(
            PartBunch_t& bunch, const std::string& solveTag, double planeZ,
            opalx::spacecharge::SelfFieldDiagnostics& diagnostics);

    /**
     * @brief Execute every solve unit emitted by one prepared iteration.
     *
     * Whole-bunch and binned plans share this host loop. Unit field mode selects the
     * monolithic electrostatic path or the per-bin Lorentz-transformed path.
     */
    void executeIterationPlan(
            PartBunch_t& bunch, IterationPlan_t& iterationPlan, const PreparedIteration_t& prepared,
            std::uint64_t particleGeneration, const PreparedCorrection_t& correction,
            opalx::spacecharge::SelfFieldDiagnostics& diagnostics);

    /**
     * @brief Execute one prepared pass for one whole-bunch or binned solve unit.
     */
    void executeSolvePass(
            PartBunch_t& bunch, const SolveUnit_t& unit, const SolvePass_t& pass,
            bool& dumpedDirichletPlaneThisStep,
            opalx::spacecharge::SelfFieldDiagnostics& diagnostics);

    /**
     * @brief Build mesh `rho` for a specific merged bin and apply all corrections.
     *
     * Steps mirror the legacy ordering from the existing OPAL-X code paths and include:
     * - dt scaling for charge scattering,
     * - mesh normalization and background subtraction (non-OPEN),
     * - Lorentz rest-frame scaling: `rho /= gammaBin`,
     * - scaling by the solver coupling constant.
     *
     * @param bunch        Bunch providing geometry and charge data.
     * @param unit         Prepared particle selection and kinematics for one solve unit.
     */
    void prepareRhoForUnit(PartBunch_t& bunch, const SolveUnit_t& unit, const SolvePass_t& pass);
};

// Reduce compile-time churn: instantiate the only supported concrete solver in one TU.
extern template class BinnedFieldSolver<double, 3>;

#endif  // OPALX_BINNED_FIELD_SOLVER_H
