// NOTE: This header contains the template declaration and doxygen documentation.
// Implementation lives in `BinnedFieldSolver.tpp`.
#ifndef OPALX_BINNED_FIELD_SOLVER_H
#define OPALX_BINNED_FIELD_SOLVER_H

#include <cmath>
#include <functional>
#include <iomanip>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "PartBunch/FieldSolver.hpp"
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
 *   gathering.
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

    struct BinKinematics {
        Vector_t<double, Dim> pmean = Vector_t<double, Dim>(0.0);
        double gammaBin             = 1.0;
    };

private:
    ScatterAttribute scatterAttribute_m;
    GatherAttribute gatherAttribute_m;
    int tablePrintFrequency_m = 0;

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
                           const double gammaBin);

    /**
     * @brief Accumulate Lorentz-transformed electric/magnetic fields into temporaries.
     *
     * The solver output field is interpreted as the bin-frame electric field `E'`
     * with `B' = 0`, and transformed to the lab frame with:
     * - `E_lab = gammaBin * E' + (gammaBin - 1) * (E' · w) * w`
     * - `B_lab = (gammaBin / c^2) * (v x E')`
     * where `v = c * pmean / gammaBin` and `w = v / |v|` (or `w = 0` if `|v| = 0`).
     *
     * @param gammaBin  Global average gamma for the merged bin.
     * @param pmean     Global average normalized momentum for the merged bin.
     * @param EtmpSP    Temporary electric field buffer for accumulation.
     * @param BtmpSP    Temporary magnetic field buffer for accumulation.
     */
    void accumulateFieldToTemp(const double gammaBin,
                               const Vector_t<double, Dim>& pmean,
                               std::shared_ptr<VField_t<T, Dim>> EtmpSP,
                               std::shared_ptr<VField_t<T, Dim>> BtmpSP);

private:

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

