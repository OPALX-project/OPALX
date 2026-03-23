#ifndef OPALX_BINNED_FIELD_SOLVER_H
#define OPALX_BINNED_FIELD_SOLVER_H

#include <cmath>
#include <functional>
#include <memory>
#include <numeric>
#include <string>
#include <iomanip>

#include "PartBunch/FieldSolver.hpp"
#include "PartBunch/PartBunch.h"
#include "Utilities/OpalException.h"

/**
 * @brief Field solver wrapper that implements the full binned self-field algorithm.
 *
 * Design choices:
 * - The solver owns no particle data. It consumes a `PartBunch` via `std::shared_ptr`
 *   so that the orchestration lives outside (e.g. `ParallelTracker`).
 * - Runtime choice: if `bunch->hasBinning()` / `bunch->getBins()` is available, the solver
 *   executes the GPU-enabled per-bin algorithm. Otherwise, it runs the legacy monolithic
 *   scatter/solve/gather path so non-binned runs remain available.
 * - Attribute modularity: scatter and gather are controlled via small enums. For now only
 *   `Q -> E` is implemented (matching current OPALX behavior).
 */
template <typename T, unsigned Dim>
class BinnedFieldSolver : public FieldSolver<T, Dim> {
    static_assert(Dim == 3, "BinnedFieldSolver currently supports Dim == 3 only.");

public:
    using PartBunch_t = PartBunch<T, Dim>;
    using ParticleCtr_t = typename PartBunch_t::ParticleContainer_t;
    using AdaptBins_t = typename PartBunch_t::AdaptBins_t;
    using BCHandler_t = BCHandler<Dim>;
    using bin_index_type = typename AdaptBins_t::bin_index_type;
    using size_type = typename AdaptBins_t::size_type;

    /**
     * @brief Which particle attribute to scatter from to build `rho`.
     *
     * Only `ChargeQ` is implemented for now.
     */
    enum class ScatterAttribute { ChargeQ };

    /**
     * @brief Which particle attribute to gather the accumulated field into.
     *
     * Only `ElectricFieldE` is implemented for now.
     */
    enum class GatherAttribute { ElectricFieldE };

    /**
     * @brief Construct a binned/legacy-compatible solver.
     *
     * @param solver    The concrete solver type name (e.g. `FFT`, `OPEN`, `CG`, `NONE`).
     * @param rho       Pointer to the mesh charge-density field storage.
     * @param E         Pointer to the mesh electric-field storage (solver output).
     * @param phi       Pointer to the potential-field storage (solver internal use).
     * @param bcHandler Shared pointer to the boundary-condition handler.
     */
    BinnedFieldSolver(std::string solver,
                      Field_t<Dim>* rho,
                      VField_t<T, Dim>* E,
                      Field_t<Dim>* phi,
                      std::shared_ptr<BCHandler_t> bcHandler)
        : FieldSolver<T, Dim>(solver, rho, E, phi, bcHandler) {
        scatterAttribute_m = ScatterAttribute::ChargeQ;
        gatherAttribute_m = GatherAttribute::ElectricFieldE;
    }

    /**
     * @brief Compute space-charge self-fields for the given bunch.
     *
     * If `bunch->getBins()` exists, runs the per-bin binned algorithm:
     * scatter->rho corrections->solve->Lorentz scaling->accumulate->gather.
     *
     * Otherwise runs the legacy monolithic algorithm (scatter all particles, solve once,
     * gather directly).
     */
    void computeSelfFields(std::shared_ptr<PartBunch_t> bunch) {
        // Major step: validate inputs and decide between binned vs legacy solver.
        if (!bunch) {
            throw OpalException("BinnedFieldSolver::computeSelfFields",
                                "Passed nullptr bunch.");
        }

        std::shared_ptr<ParticleCtr_t> pc = bunch->getParticleContainer();
        if (!pc) {
            throw OpalException("BinnedFieldSolver::computeSelfFields",
                                "Bunch particle container is not available.");
        }

        if (ippl::Comm->size() == 1 && pc->getLocalNum() <= 1) {
            pc->E = 0.0;
            return;
        }

        const bool hasBins =
            bunch->hasBinning() && (bunch->getBins() != nullptr);

        // Level-5 debug: entry logging.
        Inform dbg("BinnedFieldSolver::computeSelfFields");
        dbg << level4 << "Entry: rank=" << ippl::Comm->rank()
            << ", localParticles=" << pc->getLocalNum()
            << ", totalParticles=" << pc->getTotalNum()
            << ", hasBins=" << (hasBins ? 1 : 0)
            << ", stype=" << this->getStype() << endl;

        if (hasBins) {
            computeBinnedSelfFields(bunch);
        } else {
            computeLegacySelfFields(bunch);
        }
    }

    /// Set particle scatter attribute (extensible; default `Q`).
    void setScatterAttribute(const ScatterAttribute attr) { scatterAttribute_m = attr; }

    /// Set particle gather attribute (extensible; default `E`).
    void setGatherAttribute(const GatherAttribute attr) { gatherAttribute_m = attr; }

private:
    ScatterAttribute scatterAttribute_m;
    GatherAttribute gatherAttribute_m;

private:
    void computeBinnedSelfFields(std::shared_ptr<PartBunch_t> bunch) {
        // Major step: execute full binned self-field algorithm.
        std::shared_ptr<AdaptBins_t> bins = bunch->getBins();
        if (!bins) {
            // Defensive: runtime selection above should prevent this.
            computeLegacySelfFields(bunch);
            return;
        }

        // Major step: build and merge adaptive bins for this step.
        rebinAndPrepare(bunch, bins);

        std::shared_ptr<VField_t<T, Dim>> EtmpSP = bunch->getTempEField();
        if (!EtmpSP) {
            throw OpalException("BinnedFieldSolver::computeBinnedSelfFields",
                                "Temporary E field (Etmp) is not initialized.");
        }

        VField_t<T, Dim>& Etmp = *EtmpSP;
        // Major step: reset the lab-frame accumulation buffer.
        Etmp = 0.0;

        const bin_index_type nBins = bins->getCurrentBinCount();
        // Level-5 debug: per-step overview before entering the bin loop.
        Inform dbg("BinnedFieldSolver::computeBinnedSelfFields");
        dbg << level4 << "Binned mode: nBins=" << static_cast<int>(nBins)
            << ", stype=" << this->getStype() << endl;

        // Level-3 per-call table: gammaBin / nParticles / binNumber.
        Inform table("BinnedFieldSolver::BinStats");
        table << level3 << "BinStats | nBins=" << static_cast<int>(nBins)
              << " | stype=" << this->getStype() << endl;
        table << level3 << std::setw(9) << "bin"
              << " | " << std::setw(13) << "nParticles"
              << " | " << std::setw(15) << "gammaBin" << endl;
        for (bin_index_type binIndex = 0; binIndex < nBins; ++binIndex) {
            // Major step: process a single merged bin (gamma->rho->solve->accumulate).
            const size_type nPartGlobal = bins->getNPartInBin(binIndex, true);
            if (nPartGlobal == 0) {
                continue;
            }

            dbg << level4 << "binIndex=" << static_cast<int>(binIndex)
                << " nPartGlobal=" << static_cast<unsigned long long>(nPartGlobal) << endl;

            // Major step: compute global average gamma for this bin.
            const double gammaBin = computeGammaBinGlobal(bunch, bins, binIndex, nPartGlobal);
            if (gammaBin <= 0.0) {
                throw OpalException("BinnedFieldSolver::computeBinnedSelfFields",
                                    "Computed non-positive gamma for bin.");
            }

            dbg << level4 << "binIndex=" << static_cast<int>(binIndex)
                << " gammaBin=" << std::setprecision(10) << gammaBin << endl;

            table << level3 << std::setw(9) << static_cast<int>(binIndex)
                  << " | " << std::setw(13) << static_cast<unsigned long long>(nPartGlobal)
                  << " | " << std::setw(15) << std::setprecision(10) << gammaBin << endl;

            // Major step: build rho for this bin and apply lab->solver corrections.
            prepareRhoForBin(bunch, bins, binIndex, nPartGlobal, gammaBin);
            // Ensure deterministic accumulation even for solver types that do not update `E`.
            *(this->getE()) = 0.0;
            dbg << level4 << "binIndex=" << static_cast<int>(binIndex)
                << " runSolver(true) start" << endl;
            this->runSolver(true);
            dbg << level4 << "binIndex=" << static_cast<int>(binIndex)
                << " runSolver(true) done; accumulate->Etmp" << endl;

            accumulateFieldToTemp(gammaBin, EtmpSP);
        }

        // Major step: gather accumulated lab-frame E onto particle attributes.
        gatherFromTempToParticles(bunch, EtmpSP);
    }

    void computeLegacySelfFields(std::shared_ptr<PartBunch_t> bunch) {
        // This code is a direct move of the legacy implementation from
        // PartBunch::computeSelfFields (scatter/solve/gather for all particles).

        std::shared_ptr<ParticleCtr_t> pc = bunch->getParticleContainer();
        if (!pc) {
            throw OpalException("BinnedFieldSolver::computeLegacySelfFields",
                                "Bunch particle container is not available.");
        }

        // Level-5 debug: legacy mode entry.
        Inform dbg("BinnedFieldSolver::computeLegacySelfFields");
        dbg << level4 << "Legacy mode entry: localP=" << pc->getLocalNum()
            << ", totalP=" << pc->getTotalNum() << ", stype=" << this->getStype() << endl;

        if (scatterAttribute_m != ScatterAttribute::ChargeQ) {
            throw OpalException("BinnedFieldSolver::computeLegacySelfFields",
                                "Unsupported scatter attribute in legacy solver.");
        }

        // Level-3 legacy-mode table (binNumber=-1).
        Inform table("BinnedFieldSolver::LegacyStats");
        table << level3 << "BinStats | nBins=0 | stype=" << this->getStype() << endl;
        table << level3 << std::setw(9) << "bin"
              << " | " << std::setw(13) << "nParticles"
              << " | " << std::setw(15) << "gammaBin" << endl;
        table << level3 << std::setw(9) << -1
              << " | " << std::setw(13) << static_cast<unsigned long long>(pc->getTotalNum())
              << " | " << std::setw(15) << std::setprecision(10) << bunch->get_gamma()
              << endl;

        ippl::ParticleAttrib<T>* Qattrib = &pc->Q;
        typename PartBunch_t::Base::particle_position_type* R = &pc->R;

        Field_t<Dim>* rho = this->getRho();
        *rho = 0.0;

        // Scatter only depends on charges; scale by per-particle dt like legacy code.
        *Qattrib = (*Qattrib) * pc->dt;
        scatter(*Qattrib, *rho, *R);
        *Qattrib = (*Qattrib) / pc->dt;

        // Rho normalization for fractional time steps.
        (*rho) = (*rho) / bunch->getdT();

        const std::string stype = this->getStype();
        if (stype != "FEM" && stype != "FEM_PRECON") {
            const double cellVolume =
                std::reduce(bunch->hr_m.begin(), bunch->hr_m.end(), 1.0, std::multiplies<double>());
            (*rho) = (*rho) / cellVolume;
        }

        // Alpine uses net-0 charge for non-OPEN solvers (periodic BCs).
        if (stype != "OPEN") {
            double size = 1.0;
            for (size_t d = 0; d < Dim; ++d) {
                size *= bunch->rmax_m[d] - bunch->rmin_m[d];
            }

            const double totalQ = bunch->getCharge();
            (*rho) = (*rho) - (totalQ / size);
        }

        (*rho) = (*rho) * this->getCouplingConstant();
        // Ensure deterministic output even for solver types that do not update `E`.
        *(this->getE()) = 0.0;
        dbg << level4 << "Legacy mode: runSolver() start" << endl;
        this->runSolver();
        dbg << level4 << "Legacy mode: gather E->particles" << endl;

        // Gather solver output directly (legacy path does not use Etmp).
        if (gatherAttribute_m == GatherAttribute::ElectricFieldE) {
            gather(pc->E, *this->getE(), *R);
        } else {
            throw OpalException("BinnedFieldSolver::computeLegacySelfFields",
                                "Unsupported gather attribute in legacy solver.");
        }
    }

    void rebinAndPrepare(std::shared_ptr<PartBunch_t> bunch,
                          std::shared_ptr<AdaptBins_t> bins) {
        // Major step: adaptive histogram configuration.
        Inform dbg("BinnedFieldSolver::rebinAndPrepare");
        dbg << level4 << "Rebin start: maxBins=" << static_cast<int>(bins->getMaxBinCount()) << endl;
        bins->doFullRebin(bins->getMaxBinCount());
        bunch->dumpBinConfig(true);
        bins->sortContainerByBin();
        bins->genAdaptiveHistogram();
        bunch->dumpBinConfig(false);
        dbg << level4 << "Rebin done: currentBins=" << static_cast<int>(bins->getCurrentBinCount()) << endl;
    }

    /**
     * @brief Compute the per-bin global average gamma factor.
     *
     * We compute `gamma_i = sqrt(1 + Pz_i^2)` for particles in the (merged) bin,
     * then average across particles and MPI ranks by using the global particle count.
     */
    double computeGammaBinGlobal(std::shared_ptr<PartBunch_t> bunch,
                                  std::shared_ptr<AdaptBins_t> bins,
                                  const bin_index_type binIndex,
                                  const size_type nPartGlobal) const {
        // Major step: compute global mean gamma for the merged bin.
        Inform dbg("BinnedFieldSolver::computeGammaBinGlobal");
        dbg << level4 << "gammaBinGlobal: binIndex=" << static_cast<int>(binIndex)
            << ", nPartGlobal=" << static_cast<unsigned long long>(nPartGlobal) << endl;
        typename AdaptBins_t::position_view_type pView = bunch->getParticleContainer()->P.getView();
        typename AdaptBins_t::hash_type indices = bins->getHashArray();

        double localGammaSum = 0.0;
        Kokkos::parallel_reduce(
            "BinnedFieldSolver::gammaPerBin",
            bins->getBinIterationPolicy(binIndex),
            KOKKOS_LAMBDA(const size_type i, double& sum) {
                const double pz = pView(indices(i))[2];
                sum += Kokkos::sqrt(1.0 + pz * pz);
            },
            localGammaSum);

        double globalGammaSum = 0.0;
        ippl::Comm->allreduce(&localGammaSum, &globalGammaSum, 1, std::plus<double>());

        if (nPartGlobal == 0) {
            return 1.0;
        }
        return globalGammaSum / static_cast<double>(nPartGlobal);
    }

    /**
     * @brief Build rho on the mesh for a specific bin and apply all lab->solver corrections.
     *
     * Corrections follow the legacy ordering from `PartBunch::computeSelfFields`, plus:
     * - Lorentz rest-frame scaling: rho /= gamma_bin
     */
    void prepareRhoForBin(std::shared_ptr<PartBunch_t> bunch,
                           std::shared_ptr<AdaptBins_t> bins,
                           const bin_index_type binIndex,
                           const size_type nPartGlobal,
                           const double gammaBin) {
        // Major step: scatter bin charges and apply rho corrections (dt, cellVolume, background, Lorentz scaling).
        Inform dbg("BinnedFieldSolver::prepareRhoForBin");
        dbg << level4 << "prepareRho: binIndex=" << static_cast<int>(binIndex)
            << ", nPartGlobal=" << static_cast<unsigned long long>(nPartGlobal)
            << ", gammaBin=" << std::setprecision(10) << gammaBin << endl;
        Field_t<Dim>* rho = this->getRho();
        *rho = 0.0;

        std::shared_ptr<ParticleCtr_t> pc = bunch->getParticleContainer();
        ippl::ParticleAttrib<T>* Qattrib = &pc->Q;
        typename PartBunch_t::Base::particle_position_type* R = &pc->R;

        if (scatterAttribute_m != ScatterAttribute::ChargeQ) {
            throw OpalException("BinnedFieldSolver::prepareRhoForBin",
                                "Unsupported scatter attribute in binned solver.");
        }

        // Charge 'unit' is charge per macroparticle; scale by per-particle dt before scatter.
        *Qattrib = (*Qattrib) * pc->dt;
        scatter(*Qattrib, *rho, *R, bins->getBinIterationPolicy(binIndex), bins->getHashArray());
        *Qattrib = (*Qattrib) / pc->dt;

        (*rho) = (*rho) / bunch->getdT();

        const std::string stype = this->getStype();
        if (stype != "FEM" && stype != "FEM_PRECON") {
            const double cellVolume =
                std::reduce(bunch->hr_m.begin(), bunch->hr_m.end(), 1.0, std::multiplies<double>());
            (*rho) = (*rho) / cellVolume;
        }

        // Background subtraction for non-OPEN solvers. Here we subtract only the bin's charge.
        if (stype != "OPEN") {
            double size = 1.0;
            for (size_t d = 0; d < Dim; ++d) {
                size *= bunch->rmax_m[d] - bunch->rmin_m[d];
            }

            const double totalQBin =
                bunch->getChargePerParticle() * static_cast<double>(nPartGlobal);
            (*rho) = (*rho) - (totalQBin / size);
        }

        // Lorentz transform of charge density to the bin rest frame (thesis Eq. step 7).
        (*rho) = (*rho) / gammaBin;
        (*rho) = (*rho) * this->getCouplingConstant();
    }

    /**
     * @brief Accumulate the transformed E-field from the solver output into `Etmp`.
     *
     * Lorentz boost for a z-directed boost with `B' = 0` (thesis Eq. (12)):
     * - Ex_lab = gamma_bin * Ex'
     * - Ey_lab = gamma_bin * Ey'
     * - Ez_lab = Ez'
     */
    void accumulateFieldToTemp(const double gammaBin, std::shared_ptr<VField_t<T, Dim>> EtmpSP) {
        // Major step: transform rest-frame E to lab-frame E (z-boost) and accumulate.
        Inform dbg("BinnedFieldSolver::accumulateFieldToTemp");
        dbg << level4 << "accumulate: gammaBin=" << std::setprecision(10) << gammaBin << endl;
        const VField_t<T, Dim>& Eprime = *(this->getE());
        VField_t<T, Dim>& Etmp = *EtmpSP;

        ippl::parallel_for(
            "BinnedFieldSolver::accumulateFieldToTemp",
            Eprime.getFieldRangePolicy(),
            KOKKOS_LAMBDA(const ippl::RangePolicy<Dim>::index_array_type& idx) {
                Vector_t<T, Dim> ePrime = apply(Eprime, idx);
                Vector_t<T, Dim> eTotal = apply(Etmp, idx);

                eTotal[0] += gammaBin * ePrime[0];
                eTotal[1] += gammaBin * ePrime[1];
                eTotal[2] += ePrime[2];
                apply(Etmp, idx) = eTotal;
            });
    }

    void gatherFromTempToParticles(std::shared_ptr<PartBunch_t> bunch,
                                    std::shared_ptr<VField_t<T, Dim>> EtmpSP) {
        // Major step: gather accumulated lab-frame E from mesh back to particles.
        Inform dbg("BinnedFieldSolver::gatherFromTempToParticles");
        dbg << level4 << "gather Etmp->particles" << endl;
        VField_t<T, Dim>& Etmp = *EtmpSP;
        std::shared_ptr<ParticleCtr_t> pc = bunch->getParticleContainer();
        if (gatherAttribute_m == GatherAttribute::ElectricFieldE) {
            gather(pc->E, Etmp, pc->R);
        } else {
            throw OpalException("BinnedFieldSolver::gatherFromTempToParticles",
                                "Unsupported gather attribute in binned solver.");
        }
    }
};

#endif  // OPALX_BINNED_FIELD_SOLVER_H
