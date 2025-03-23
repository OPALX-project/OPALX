/**
 * @file AdaptBins.h
 * @brief Defines a structure to hold particles in energy bins and their associated data.
 * 
 * "AdaptBins" is a improved version of "PartBins" from the old OPAL. AdaptBins uses
 * data structures from IPPL and Kokkos to make improve performance when using MPI and Cuda.
 * In contrast to the old PartBin, this class allows for the use of rebinning during runtime. 
 */

/**
 * @copyright Copyright (c) 2007-2020, Paul Scherrer Institut, Villigen PSI, Switzerland
 * All rights reserved
 *
 * This file is part of OPALX.
 *
 * OPAL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU General Public License
 * along with OPAL. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef ADAPT_BINS_H
#define ADAPT_BINS_H

//#include <memory>
//#include <iostream>
#include "Ippl.h"

#include <Kokkos_DualView.hpp>   // Looks like no one ever used it before, so import it here...
#include "ParallelReduceTools.h" // Has custom reducer objects (--> needed in AdaptBins.h and BinningTools.h)
#include "BinningTools.h"        // Has custom particle selection
#include "BinHisto.h"           // Has custom histogram class for DualView management

namespace ParticleBinning {

    /**
     * @class AdaptBins
     * @brief A class that bins particles in energy bins and allows for runtime rebinning.
     * 
     * This class provides functionality to group particles into energy bins, initialize and reinitialize
     * bin structures, and update bin contents using data structures from IPPL and Kokkos. It is optimized
     * for usage with MPI and CUDA backends.
     *
     * @tparam BunchType The type of particle bunch (container) used in the binning process.
     */
    template <typename BunchType, typename BinningSelector>
    class AdaptBins {
    public:
        using value_type             = typename BinningSelector::value_type;
        using particle_position_type = typename BunchType::particle_position_type;
        using position_view_type     = typename particle_position_type::view_type;
        using size_type              = typename BunchType::size_type;
        using bin_index_type         = typename BunchType::bin_index_type;
        using bin_type               = typename ippl::ParticleAttrib<bin_index_type>;
        using bin_view_type          = typename bin_type::view_type;
        // using bin_histo_type         = Kokkos::View<size_type*>;
        //using bin_host_histo_type    = Kokkos::View<size_type*, Kokkos::HostSpace>;
        // using bin_histo_dual_type    = typename Kokkos::DualView<size_type*>;
        // using binning_var_selector_type = typename BinningVariableSelector<size_type>;
        // using buffer_view_type       = Kokkos::View<int*>;
        using hash_type              = ippl::detail::hash_type<Kokkos::DefaultExecutionSpace::memory_space>;

        //using histo_type      = Histogram<size_type, bin_index_type, value_type>;
        using d_histo_type          = Histogram<size_type, bin_index_type, value_type, true>;
        using dview_type            = typename d_histo_type::dview_type;
        using hview_type            = typename d_histo_type::hview_type;
        using dwidth_view_type      = typename d_histo_type::dwidth_view_type;
        using hindex_transform_type = typename d_histo_type::hindex_transform_type;
        using dindex_transform_type = typename d_histo_type::dindex_transform_type;

        using h_histo_type_g         = Histogram<size_type, bin_index_type, value_type, false, Kokkos::HostSpace>; 
        using hview_type_g           = typename h_histo_type_g::hview_type;
        using hindex_transform_type_g= typename h_histo_type_g::hindex_transform_type;

        /**
         * @brief Constructs an AdaptBins object with a specified maximum number of bins.
         * 
         * @param bunch A shared pointer to the particle container.
         * @param maxBins The maximum number of bins to initialize with (default is 10).
         */
        AdaptBins(std::shared_ptr<BunchType> bunch, BinningSelector var_selector, bin_index_type maxBins,
                  value_type binningAlpha, value_type binningBeta, value_type desiredWidth)
            : bunch_m(bunch)
            , var_selector_m(var_selector)
            , maxBins_m(maxBins)
            , binningAlpha_m(binningAlpha)
            , binningBeta_m(binningBeta)
            , desiredWidth_m(desiredWidth) {

            currentBins_m   = maxBins; // TODO for now...
            // sortingBuffer_m = Kokkos::View<int*>("particlePermutationBuffer", bunch->getLocalNum());

            // Initialize all the timers
            bInitLimitsT            = IpplTimings::getTimer("bInitLimits");
            bAllReduceLimitsT       = IpplTimings::getTimer("bAllReduceLimits");
            bAllReduceGlobalHistoT  = IpplTimings::getTimer("bAllReduceGlobalHisto");
            bAssignUniformBinsT     = IpplTimings::getTimer("bAssignUniformBins");
            bExecuteHistoReductionT = IpplTimings::getTimer("bExecuteHistoReduction");
            bSortContainerByBinT    = IpplTimings::getTimer("bSortContainerByBin");
            bVerifySortingT         = IpplTimings::getTimer("bVerifySorting");

            Inform msg("AdaptBins");
            msg << "AdaptBins initialized with maxBins = " << maxBins_m 
                << ", alpha = " << binningAlpha_m
                << ", beta = " << binningBeta_m
                << ", desiredWidth = " << desiredWidth_m << endl;
        }

        /**
         * @brief Returns a view to the particle bin array.
         * 
         * @note How the bin attribute it named might change, so this function exists!
         */
        bin_view_type getBinView() { return bunch_m->Bin.getView(); }

        /**
         * @brief Gets the current number of bins.
         * 
         * @return The current bin count.
         */
        bin_index_type getCurrentBinCount() const { return currentBins_m; }

        /**
         * @brief Gets the maximum number of bins.
         * 
         * @return The maximum allowed number of bins.
         */
        bin_index_type getMaxBinCount() const { return maxBins_m; }

        /**
         * @brief Gets the average binwidth.
         * 
         * @return Corresponds to (xmax_m - xmin_m)/n_bins.
         */
        value_type getBinWidth() const { return binWidth_m; }

        /**
         * @brief Sets the current number of bins and adjusts the bin width.
         * 
         * @param nBins The new number of bins.
         */
        void setCurrentBinCount(bin_index_type nBins) {
            currentBins_m = (nBins > maxBins_m) ? maxBins_m : nBins; 
            binWidth_m    = (xMax_m - xMin_m) / currentBins_m; // assuming particles did not change!
        }

        hash_type getHashArray() { return sortedIndexArr_m; }

        /**
         * @brief Calculates the bin index for a given position value.
         * 
         * This static method calculates which bin a position value falls into based on the bin
         * boundaries and bin width.
         * 
         * @param x The "position" value.
         * @param xMin Minimum boundary for the bins.
         * @param xMax Maximum boundary for the bins.
         * @param binWidthInv Inverse of the bin width for efficiency.
         * @param numBins The total number of bins.
         * @return The calculated bin index for the position value.
         */
        KOKKOS_INLINE_FUNCTION
        static bin_index_type getBin(value_type x, value_type xMin, value_type xMax, value_type binWidthInv, bin_index_type numBins);

        /**
         * @brief Initializes the limits for binning based on the particle data.
         * 
         * This function calculates the minimum and maximum limits (xMin and xMax) from the
         * particle positions, which are then used to define bin boundaries.
         * 
         * @note Needs to be called _before_ bins and histograms are initialized.
         */
        void initLimits();

        /**
         * @brief Initializes the histogram view for binning and optionally sets it to zero.
         * 
         * @param setToZero If true, initializes the histogram view to zero. Default is false. The 0 initialization is not needed if it is overwritten anyways.
         */
        void instantiateHistogram(bool setToZero = false);

        /**
         * @brief Assigns each particle in the bunch to a bin based on its position.
         * 
         * This function iterates over all particles in the bunch, calculates their bin
         * index, and updates the bin structure accordingly.
         */
        void assignBinsToParticles();

        /**
         * @brief Initializes a local histogram view for particle binning.
         * 
         * This function prepares a local histogram, enabling binning of particles within
         * local data, which can then be reduced into a global histogram.
         * Calls executeInitLocalHistoReduction to perform the reduction (has more information).
         */
        void initLocalHisto(HistoReductionMode modePreference = HistoReductionMode::Standard);

        /**
        * @brief Initializes and performs a team-based histogram reduction for particle bins.
        *
        * This function allocates scratch memory on each team, initializes a local
        * histogram for each team in shared memory, updates it based on bin indices of particles,
        * and finally reduces the team-local histograms into a global histogram in device memory
        * using pure atomics.
        * 
        * @details The process consists of the following steps:
        * - Allocating scratch memory for each team's local histogram.
        * - Initializing the local histogram to zero.
        * - Assigning particles to bins in parallel within each team.
        * - Reducing each team's local histogram into a global histogram (atomics).
        *
        * ### Parameters
        * - **binIndex**: A view of bin indices for each particle.
        * - **localBinHisto**: A global histogram where the final reduction result is stored.
        * - **binCount**: The total number of histogram bins.
        * - **localNumParticles**: The number of particles in the local process.
        *
        * ### Memory and Execution
        * - **Scratch Memory**: Scratch memory is allocated per team for a local histogram, with size `binCount`.
        * - **Concurrency**: `team_size` specifies the number of threads per team, and each team processes a `block_size`.
        *
        * @note This function is optimized for GPU execution using team-based parallelism,
        *       it does not work on Host (since team_size is hardcoded and to big). 
        *       If you want to run this on Host, change team_size=1 and increase block_size.
        *
        * @pre `localBinHisto` and `binIndex` must be initialized with appropriate sizes before calling this function.
        * 
        * @post `localBinHisto` contains the reduced histogram for the local data. 
        *       Next step is to reduce across all MPI ranks. @see getGlobalHistogram
        */
        void executeInitLocalHistoReductionTeamFor();

        /**
        * @brief Executes a parallel reduction to initialize the local histogram for particle bins.
        *
        * This function performs a Kokkos parallel reduction over the particles in the bunch, incrementing 
        * counts in the reduction array `to_reduce` based on the bin index for each particle. 
        * After the reduction, the results are copied to the final histogram `localBinHisto`.
        *
        * @tparam ReducerType The type of the reduction object, which should support `the_array` for bin counts.
        * @param to_reduce A reduction object that accumulates bin counts for the histogram.
        *
        * The function performs the following steps:
        * - Executes a Kokkos parallel reduction loop where each particle increments the bin count 
        *   corresponding to its bin index.
        * - Executes a parallel loop to copy the reduced bin counts from `to_reduce` to `localBinHisto`.
        *
        * @note This function uses the Kokkos parallel programming model and assumes that `to_reduce` 
        * has a `the_array` member which stores the histogram counts. So far, only ParticleBinning::ArrayReduction
        * is implemented in that way (to work together with Kokkos::Sum reducer). `the_array` needs to have a known
        * size at compile time.
        */
        template<typename ReducerType>
        void executeInitLocalHistoReduction(ReducerType& to_reduce);

        /**
         * @brief Retrieves the global histogram across all processes.
         * 
         * This function reduces the local histograms across all MPI processes into
         * a single global histogram view.
         * 
         * @return A view of the global histogram in host space (used for debugging output).
         */
        void initGlobalHistogram();

        /**
         * @brief Performs a full rebinning of particles with a specified number of bins.
         * 
         * @param nBins The new number of bins to use for rebinning.
         */
        void doFullRebin(bin_index_type nBins, bool recalculateLimits = true, HistoReductionMode modePreference = HistoReductionMode::Standard) {
            if (recalculateLimits) initLimits();
            setCurrentBinCount(nBins);
            assignBinsToParticles();
            initHistogram(modePreference);
        }

        void initHistogram(HistoReductionMode modePreference = HistoReductionMode::Standard) {
            instantiateHistogram(true); // Init histogram (no need to set to 0, since executeInitLocalHistoReduction overwrites values from reduction...) --> true, since it is necessary for atomics option...
            initLocalHisto(modePreference);
            //initLocalPostSum();
            initGlobalHistogram();

            // Init both histograms --> buils postSums and widths arrays
            localBinHisto_m.init();
            globalBinHisto_m.init();
        }

        /**
         * @brief Initializes the prefix sum for the local histogram (used e.g. for sorting and indexing in `scatter(...)`).
         */
        /*void initLocalPostSum() { 
            localBinHistoPostSum_m = bin_histo_dual_type("localBinHistoPostSum_m", getCurrentBinCount() + 1);

            computeFixSum<size_type>(localBinHisto_m.view_device(), localBinHistoPostSum_m.view_device());
            localBinHistoPostSum_m.modify_device(); 
            localBinHistoPostSum_m.sync_host(); 

            // Print the prefix sum for debugging
            Inform msg("AdaptBins");
            auto localPostSumHost = localBinHistoPostSum_m.view_host();
            for (bin_index_type i = 0; i < getCurrentBinCount() + 1; ++i) {
                msg << "PrefixSum[" << i << "] = " << localPostSumHost(i) << endl;
            }
        }*/

        size_type getNPartInBin(bin_index_type binIndex, bool global = false) {
            /**
             * Assume DualView was properly synchronized.
             * Might create some overhead from .view_host() call if called often (solved: not anymore with DualView).
             * However, it is only called on host (max nBins times per iteration), so should be fine. You can make it
             * more efficient by avoiding the Kokkos:View "copying-action" with e.g. dualView.h_view(binIndex)
             */
            if (binIndex < 0 || binIndex >= getCurrentBinCount()) { return bunch_m->getTotalNum(); } // shouldn't happen..., "binIndex < 0" unnecessary, since binIndex is usually unsigned
            if (global) {
                return globalBinHisto_m.getNPartInBin(binIndex);
            } else {
                return localBinHisto_m.getNPartInBin(binIndex);
            }
        }

        void sortContainerByBin();

        Kokkos::RangePolicy<> getBinIterationPolicy(const bin_index_type& binIndex) {
            return localBinHisto_m.getBinIterationPolicy(binIndex);
            //auto localPostSumHost = localBinHistoPostSum_m.view_host();
            //return Kokkos::RangePolicy<>(localPostSumHost(binIndex), localPostSumHost(binIndex + 1));
        }

        void genAdaptiveHistogram() {
            // 1. Run merging algorithm on globalHisto --> generates global binWidths array and postSum.
            // Note: Assumes that the histograms are properly initialized.
            //double tmp_ratio                        = (xMax_m - xMin_m) * sqrt(bunch_m->getTotalNum()) / 20; // (xMax_m - xMin_m) / 10 * (xMax_m - xMin_m)/bunch_m->getTotalNum(); // should be ~10 bins // bunch_m->getTotalNum()  /  
            // double tmp_ratio                        = bunch_m->getTotalNum() / 10; // (xMax_m - xMin_m) / 10;

            var_selector_m.updateDataArr(bunch_m); // Probably not necessary, since it is called before updating particles; but just in case!
            hindex_transform_type adaptLookup       = globalBinHisto_m.mergeBins(/*sortedIndexArr_m, var_selector_m*/);
            dindex_transform_type adaptLookupDevice = dindex_transform_type("adaptLookupDevice", currentBins_m);
            Kokkos::deep_copy(adaptLookupDevice, adaptLookup);
            
            bin_view_type binIndex                  = getBinView();

            setCurrentBinCount(globalBinHisto_m.getCurrentBinCount());

            // 2. Map old indices to the new histogram ("Rebin")
            Kokkos::parallel_for("RebinParticles", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const size_type& i) {
                bin_index_type oldBin = binIndex(i);
                binIndex(i) = adaptLookupDevice(oldBin);
            });

            // 3. Update local histogram with new indices
            instantiateHistogram(true);
            initLocalHisto(); // Runs reducer on new bin indices (also does the sync)

            localBinHisto_m.initPostSum(); // only init postsum, since the widths are not constant anymore
            localBinHisto_m.copyBinWidths(globalBinHisto_m);
        }

        /**
         * @brief Prints the current global histogram to the Inform output stream.
         * 
         * This function outputs the global histogram data (bin counts) to the standard output.
         * Note: Only works correctly for rank 0 in an MPI environment.
         */
        void print() {
            /*Inform os("AdaptBins");
            // Only works correct for rank 0
            os << "-----------------------------------------" << endl;
            os << "     Output Global Binning Structure     " << endl;

            bin_index_type numBins = getCurrentBinCount();
            os << "Bins = " << numBins << " hBin = " << binWidth_m << endl;
            os << "Bin #;Val" << endl;

            hview_type globalHostHisto = globalBinHisto_m.template getHostView<hview_type>(globalBinHisto_m.getHistogram());

            // Only rank 0 prints the global histogram
            size_type total = 0;
            for (bin_index_type i = 0; i < numBins; ++i) {
                size_type val = globalHostHisto(i); // Can do it like this, since DualView knows it is on host
                os << i << ";" << val << endl;
                total += val;
            }   
            os << "Total = " << total << endl;
            os << "-----------------------------------------" << endl;*/

            // TODO
            
            // globalBinHisto_m.printHistogram();

            globalBinHisto_m.printPythonArrays();
        }

        /**
         * @brief Outputs debug information related to Kokkos and MPI configurations.
         * 
         * This function prints information about the number of threads (in OpenMP) or GPUs
         * (in CUDA) available on the current MPI rank, along with other debug information.
         */
        void debug() {
            Inform msg("KOKKOS DEBUG"); // , INFORM_ALL_NODES

            int rank = ippl::Comm->rank();
            msg << "=====================================" << endl;
            msg << " Kokkos Debug Information (Rank " << rank << ")" << endl;
            msg << "=====================================" << endl;

            // Check number of CPU threads (OpenMP or other CPU execution spaces)
            #ifdef KOKKOS_ENABLE_OPENMP
            int num_threads = Kokkos::OpenMP::concurrency();
            msg << "CPU Threads (OpenMP): " << num_threads << endl;
            #elif defined(KOKKOS_ENABLE_THREADS)
            int num_threads = Kokkos::Threads::concurrency();
            msg << "CPU Threads (Kokkos::Threads): " << num_threads << endl;
            #else
            msg << "CPU Threads: No multi-threaded CPU execution space enabled." << endl;
            #endif

            // Check number of GPUs (CUDA devices)
            #ifdef KOKKOS_ENABLE_CUDA
            int num_gpus = Kokkos::Cuda::detect_device_count();
            msg << "CUDA Enabled: Rank " << rank << " sees " << num_gpus << " GPU(s) available." << endl;
            Kokkos::Cuda cuda_instance;  
            std::stringstream ss;
            cuda_instance.print_configuration(ss);
            msg << ss.str();
            #else
            msg << "CUDA: GPU support disabled.\n";
            #endif

            // Additional information on concurrency in the default execution space
            int default_concurrency = Kokkos::DefaultExecutionSpace::concurrency();
            msg << "Default Execution Space Concurrency: " << default_concurrency << endl;
            msg << "Binning cost function parameters: alpha = " << binningAlpha_m << ", beta = " << binningBeta_m << ", desiredWidth = " << desiredWidth_m << endl;

            msg << "=====================================" << endl;
        }

        /**
        * @brief Applies a Lorentz transformation to a given vector field based on particle velocities.
        * 
        * @tparam T The data type of the field components (e.g., `double`, `float`).
        * @tparam Dim The dimensionality of the field.
        * @param field A reference to the vector field to be transformed.
        * @return A reference to the transformed field.
        * 
        * @details
        * - **Gamma Factor Calculation:** The gamma factor is derived from the velocity components of particles.
        *   For a given particle velocity \(v\), the gamma factor is computed as:
        *   \[
        *   \gamma = \frac{1}{\sqrt{1 - \vec{v} \cdot \vec{v}}}
        *   \]
        *   where \(\vec{v}\) is the velocity vector of a particle.
        * 
        * - **Field Transformation:** After computing the gamma factor, each component of the field is
        *   multiplied by the corresponding gamma factor.
        * 
        * ### Example Usage:
        * ```cpp
        * VField_t<double, 3> field = ...; // Initialize the field
        * this->bins_m->LTrans(field);     // Apply Lorentz transformation
        * ```
        */
        template <typename T, unsigned Dim>
        VField_t<T, Dim>& LTrans(VField_t<T, Dim>& field, const bin_index_type& currentBin); // TODO: may want to add usage of c constant when it exists...

    private:
        std::shared_ptr<BunchType> bunch_m;    ///< Shared pointer to the particle container.
        BinningSelector var_selector_m;        ///< Variable selector for binning.
        const bin_index_type maxBins_m;              ///< Maximum number of bins.
        bin_index_type currentBins_m;          ///< Current number of bins in use.
        value_type xMin_m;                     ///< Minimum boundary for bins.
        value_type xMax_m;                     ///< Maximum boundary for bins.
        value_type binWidth_m;                 ///< Width of each bin.

        value_type binningAlpha_m;               ///< Alpha parameter for binning cost function. 
        value_type binningBeta_m;                ///< Beta parameter for binning cost function.
        value_type desiredWidth_m;               ///< Desired bin width for binning cost function.

        // Histograms 
        d_histo_type localBinHisto_m;          ///< Local histogram view for bin counts.
        h_histo_type_g globalBinHisto_m;         ///< Global histogram view (over ranks reduced local histograms).


        //bin_histo_dual_type localBinHisto_m;          ///< Local histogram view for bin counts.
        //bin_histo_dual_type localBinHistoPostSum_m; ///< Local prefix sum view for bin counts.
        //bin_host_histo_type localBinHistoHost_m; // TODO: Use DualView instead!!!!
        //bin_histo_dual_type globalBinHisto_m;         ///< Global histogram view (over ranks reduced local histograms).
        hash_type sortedIndexArr_m;                  ///< Hash table for sorting particles by bin index.
        // buffer_view_type sortingBuffer_m;      ///< Buffer for permutating particles after sorting by bin index.

        // Here are all the timer
        IpplTimings::TimerRef bInitLimitsT;
        IpplTimings::TimerRef bAllReduceLimitsT;
        IpplTimings::TimerRef bAllReduceGlobalHistoT;
        IpplTimings::TimerRef bAssignUniformBinsT;
        IpplTimings::TimerRef bExecuteHistoReductionT;
        IpplTimings::TimerRef bSortContainerByBinT;
        IpplTimings::TimerRef bVerifySortingT;
    };

}

#include "AdaptBins.tpp"

#endif  // ADAPT_BINS_H


