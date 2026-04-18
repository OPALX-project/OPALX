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

#include "PartBunch/FieldSolver.hpp"
#include "PartBunch/ImageChargeScatterController.h"
#include "PartBunch/PartBunch.h"
#include "Utilities/OpalException.h"

/**
 * @brief Field solver wrapper that implements the full binned self-field algorithm.
 *
 * @tparam T   Particle numeric type (typically `double`).
 * @tparam Dim Spatial dimensionality (currently only `Dim == 3` is supported).
 *
 * Design overview:
 * - The solver owns no particle data; it consumes a `PartBunch` via `std::shared_ptr`.
 * - Runtime selection:
 *   - If adaptive bins are available (`bunch->hasBinning()` / `bunch->getBins()`),
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

public:
    using PartBunch_t   = PartBunch<T, Dim>;
    using ParticleCtr_t = typename PartBunch_t::ParticleContainer_t;
    using AdaptBins_t   = typename PartBunch_t::AdaptBins_t;
    using BCHandler_t   = BCHandler<Dim>;
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
     * @brief Controls which charges are scattered during rho preparation.
     *
     * `PrimaryAndImage` scatters both (legacy combined behavior).
     * `PrimaryOnly` scatters only the real bunch charges.
     * `ImageOnly` scatters only the mirrored image charges.
     */
    enum class ImageScatterMode { PrimaryAndImage, PrimaryOnly, ImageOnly };

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
     * @param rho        Pointer to the mesh charge-density field storage.
     * @param E          Pointer to the mesh electric-field storage (solver output).
     * @param phi        Pointer to the potential-field storage (solver internal use).
     * @param bcHandler  Shared pointer to the boundary-condition handler.
     * @param tablePrintFrequency Global timestep frequency of printing the bin stats table
     *                             to console in binned mode. If `0`, printing is disabled.
     */
    BinnedFieldSolver(std::string solver,
                       Field_t<Dim>* rho,
                       VField_t<T, Dim>* E,
                       Field_t<Dim>* phi,
                       std::shared_ptr<BCHandler_t> bcHandler,
                       int tablePrintFrequency);

    /**
     * @brief Compute space-charge self-fields for the given particle bunch.
     *
     * If the bunch provides adaptive binning (`bunch->getBins()`), the solver executes
     * the per-bin algorithm:
     * `scatter rho corrections -> solve -> Lorentz scaling -> accumulate -> gather`.
     * Otherwise, it executes the legacy monolithic algorithm:
     * `scatter all particles -> solve once -> gather directly`.
     *
     * @param bunch Pointer to the particle bunch (must not be null).
     *
     * @throws OpalException If required internal data (particle container / temp E field)
     *                        is missing, or if unsupported scatter/gather modes are selected.
     */
    void computeSelfFields(std::shared_ptr<PartBunch_t> bunch);

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
    bool isImageChargeEnabled() const { return imageScatterController_m.isEnabled(); }
    double getImageChargePlaneZ() const { return imageScatterController_m.getZPlane(); }

    /// @brief Configure the shifted Green's function Dirichlet correction (alternative to
    /// image charges). Mutually exclusive with @c setImageChargeConfiguration(true, ...).
    /// Requires the OPEN field solver; the solver-type check happens at runtime in
    /// @c FieldSolver::runShiftedOpenSolver when the correction pass fires.
    void setShiftedGreensConfiguration(bool enabled, double zPlane);
    bool isShiftedGreensEnabled() const { return shiftedGreensEnabled_m; }
    double getShiftedGreensPlaneZ() const { return shiftedGreensPlaneZ_m; }

    /// @brief Set the maximum number of timesteps for which image charges are active (0 = unlimited).
    void setZerofaceMaxSteps(int maxSteps);
    int getZerofaceMaxSteps() const { return zerofaceMaxSteps_m; }

    /// @brief Check whether the explicit image-charge pass should run for a given timestep.
    bool isImageChargeActiveForStep(size_t step) const;

    /// @brief Check whether the shifted Green's function correction should run for a given timestep.
    /// Reuses the same step budget (@c zerofaceMaxSteps_m) as the image-charge path.
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
    int tablePrintFrequency_m = 0;
    int zeroFacePlaneDumpFrequency_m = 0;
    int zerofaceMaxSteps_m = 0;
    ImageChargeScatterController<T, Dim> imageScatterController_m;
    bool warnedPlaneDumpParallelUnsupported_m = false;

    // Shifted Green's function Dirichlet correction (alternative to image charges).
    // Mutually exclusive with the image-charge path (enforced at config time).
    bool shiftedGreensEnabled_m    = false;
    double shiftedGreensPlaneZ_m   = 0.0;

    // Scratch view holding the z-axis-flipped local slab of E' for the shifted-GF
    // correction pass. Populated by buildFlippedZSlab before accumulateFieldToTemp's
    // flipped branch. Shape matches the local view of *(this->getE()) (incl. ghosts).
    // Allocated lazily on first use; reused across bins and timesteps.
    using FlippedView_t = Kokkos::View<Vector_t<T, Dim>***>;
    FlippedView_t flippedZSlab_m;

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
    void printBinStatsTable(const std::string& binningCmdName,
                             const std::vector<BinStatsRow>& rows);

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
            std::shared_ptr<PartBunch_t> bunch, const std::string& solveTag);

    /**
     * @brief Compute self-fields using the binned algorithm.
     *
     * Requires that the bunch has a valid bin structure and a temporary electric field
     * buffer (`bunch->getTempEField()`).
     *
     * @param bunch Particle bunch for which to compute self-fields.
     */
    void computeBinnedSelfFields(std::shared_ptr<PartBunch_t> bunch);

    /**
     * @brief Compute self-fields using the legacy monolithic algorithm.
     *
     * This is a direct adaptation of the legacy implementation:
     * scatter (all particles) -> solve once -> gather electric field directly to particles.
     *
     * @param bunch Particle bunch for which to compute self-fields.
     */
    void computeLegacySelfFields(std::shared_ptr<PartBunch_t> bunch);

    /**
     * @brief Build and prepare adaptive bins for the current step.
     *
     * This performs a full rebin to the maximum bin count, sorts the particle container by bin,
     * generates the adaptive histogram (merging), and restores bin configuration.
     *
     * @param bunch Bunch whose bins are updated.
     * @param bins  Adaptive bin structure owned/managed by the bunch.
     */
    void rebinAndPrepare(std::shared_ptr<PartBunch_t> bunch,
                          std::shared_ptr<AdaptBins_t> bins);

public:

    /**
     * @brief Compute per-bin global mean momentum and gamma.
     *
     * The function computes the global mean momentum vector `pmean` across all
     * particles in the merged bin and derives:
     * `gammaBin = sqrt(1 + dot(pmean, pmean))`.
     *
     * @param bunch        Bunch providing particle data.
     * @param bins         Bins providing the merged-bin iteration policy and indexing.
     * @param binIndex     Bin index in the merged-bin space.
     * @param nPartGlobal  Global particle count for this merged bin.
     *
     * @return Per-bin kinematics bundle (`pmean`, `gammaBin`).
     */
    BinKinematics computeGammaBinGlobal(std::shared_ptr<PartBunch_t> bunch,
                                        std::shared_ptr<AdaptBins_t> bins,
                                        const bin_index_type binIndex,
                                        const size_type nPartGlobal) const;

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
     * @param bins         Adaptive bins providing bin iteration and hash indexing.
     * @param binIndex     Merged bin index.
     * @param nPartGlobal  Global number of particles in that merged bin.
     * @param gammaBin     Global average gamma for that merged bin.
     */
    void prepareRhoForBin(std::shared_ptr<PartBunch_t> bunch,
                           std::shared_ptr<AdaptBins_t> bins,
                           const bin_index_type binIndex,
                           const size_type nPartGlobal,
                           const double gammaBin,
                           ImageScatterMode mode = ImageScatterMode::PrimaryAndImage);

    /**
     * @brief Accumulate Lorentz-transformed electric/magnetic fields into temporaries.
     *
     * The solver output field is interpreted as the bin-frame electric field `E'`
     * with `B' = 0`, and transformed to the lab frame with:
     * - `E_lab = gammaBin * E' + (gammaBin - 1) * (E' · w) * w`
     * - `B_lab = (gammaBin / c^2) * (v x E')`
     * where `v = c * pmean / gammaBin` and `w = v / |v|` (or `w = 0` if `|v| = 0`).
     *
     * When @p flipAxis is a valid axis index (0..Dim-1), the solver output is
     * read at an axis-flipped source index @c N-1-k (in that axis only) and the
     * transverse E components are negated before the Lorentz transform, yielding
     * the image-charge contribution produced by a shifted-Green's-function solve:
     * - `E_d -> -E_d` for `d != flipAxis`
     * - `E_d -> +E_d` for `d == flipAxis` (the component parallel to the flip).
     * This is derived from `phi_image(x) = -phi_shifted(R(x))` and `E = -grad(phi)`.
     * The zero-copy read from the flipped source index is the reason this is
     * baked into @c accumulateFieldToTemp instead of an out-of-place kernel.
     *
     * Currently single-rank only when @p flipAxis >= 0 (an axis-flip across MPI
     * ranks would need extra communication). The caller is expected to guard
     * against multi-rank use upstream.
     *
     * @param gammaBin  Global average gamma for the merged bin.
     * @param pmean     Global average normalized momentum for the merged bin.
     * @param EtmpSP    Temporary electric field buffer for accumulation.
     * @param BtmpSP    Temporary magnetic field buffer for accumulation.
     * @param bFieldSign +1 for forward-moving charges, -1 for image charges.
     * @param flipAxis  Axis in which to flip the read index (use -1 for no flip).
     */
    void accumulateFieldToTemp(const double gammaBin,
                               const Vector_t<double, Dim>& pmean,
                               std::shared_ptr<VField_t<T, Dim>> EtmpSP,
                               std::shared_ptr<VField_t<T, Dim>> BtmpSP,
                               double bFieldSign = 1.0,
                               int flipAxis = -1);

private:
    /// @brief Populate @c flippedZSlab_m with the z-axis globally-flipped version of @p src.
    ///
    /// Under `PARFFTZ=true` the global flip `k -> N_z_global-1-k` generally crosses MPI ranks.
    /// This helper does one pairwise-exchange pass over `ippl::Comm`: each rank packs the z-slabs
    /// of @p src that peers need, posts `MPI_Isend`/`MPI_Irecv`, and unpacks the received slabs
    /// into the correct local destination indices of @c flippedZSlab_m. After the call the
    /// lambda in `accumulateFieldToTemp` reads @c flippedZSlab_m(i, j, k) directly without any
    /// cross-rank access.
    ///
    /// @param src  Source vector field (typically `*(this->getE())` after the shifted-GF solve).
    ///             Only the z axis is flipped; x and y stay local to the rank.
    void buildFlippedZSlab(const VField_t<T, Dim>& src);

    /**
     * @brief Gather the accumulated lab-frame fields from temporaries back to particles.
     *
     * @param bunch   Particle bunch to gather into.
     * @param EtmpSP  Temporary electric field buffer holding the accumulated lab-frame field.
     * @param BtmpSP  Temporary magnetic field buffer holding the accumulated lab-frame field.
     */
    void gatherFromTempToParticles(std::shared_ptr<PartBunch_t> bunch,
                                   std::shared_ptr<VField_t<T, Dim>> EtmpSP,
                                   std::shared_ptr<VField_t<T, Dim>> BtmpSP);
};

// Reduce compile-time churn: instantiate the only supported concrete solver in one TU.
extern template class BinnedFieldSolver<double, 3>;

#endif  // OPALX_BINNED_FIELD_SOLVER_H

