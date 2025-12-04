/**
 * @file AdaptBins.h
 * @brief Defines a structure to hold particles in energy bins and their associated data.
 * 
 * "AdaptBins" is developed for OPAL-X and replaces the "PartBins" clas from the old OPAL.
 * In contrast to the old PartBin, this class re-bins adaptively on the fly without
 * copying data. 
 */

#ifndef ADAPT_BINS_H
#define ADAPT_BINS_H

#include "Ippl.h"

#include <Kokkos_DualView.hpp>
#include "ParallelReduceTools.h" 
#include "BinningTools.h"        
#include "BinHisto.h"           

namespace ParticleBinning {

    /**
     * @class AdaptBins
     * @brief A class that bins particles in energy bins and allows for adaptive runtime rebinning.
     * 
     * This class provides functionality to group particles into energy bins, initialize and reinitialize
     * bin structures, and update bin contents using data structures from IPPL and Kokkos. It is optimized
     * for usage on CPU and GPU. It is explained in detail in Alexander Liemen's master thesis and was
     * developed for OPAL-X to improve the performance when solving the electric field in a relativistic
     * particle bunch.  
     *
     * @tparam BunchType The type of particle bunch (container) used in the binning process.
     * @tparam BinningSelector A function or functor that selects the variable for binning (based on attributes).
     */
    template <typename BunchType, typename BinningSelector>
    class AdaptBins {
    public:
        // Variable type definitions
        using value_type             = typename BinningSelector::value_type;
        using particle_position_type = typename BunchType::particle_position_type;
        using position_view_type     = typename particle_position_type::view_type;
        using size_type              = typename BunchType::size_type;
        using bin_index_type         = typename BunchType::bin_index_type;
        using bin_type               = typename ippl::ParticleAttrib<bin_index_type>;
        using bin_view_type          = typename bin_type::view_type;
        using hash_type              = ippl::detail::hash_type<Kokkos::DefaultExecutionSpace::memory_space>;

        // Types from BinHisto.h to manage histogram structures without "auto"
        using d_histo_type          = Histogram<size_type, bin_index_type, value_type, true>;
        using dview_type            = typename d_histo_type::dview_type;
        using hview_type            = typename d_histo_type::hview_type;
        using dwidth_view_type      = typename d_histo_type::dwidth_view_type;
        using hindex_transform_type = typename d_histo_type::hindex_transform_type;
        using dindex_transform_type = typename d_histo_type::dindex_transform_type;

        using h_histo_type_g          = Histogram<size_type, bin_index_type, value_type, false, Kokkos::HostSpace>; 
        using hview_type_g            = typename h_histo_type_g::hview_type;
        using hindex_transform_type_g = typename h_histo_type_g::hindex_transform_type;

        /**
         * @brief Constructs an AdaptBins object with a specified maximum number of bins and a selector.
         * 
         * @param bunch A shared pointer to the underlying particle container.
         * @param var_selector A function or functor that selects the variable for binning (based on attributes).
         * @param maxBins The maximum number of bins to initialize with (usually 128).
         * @param binningAlpha The alpha parameter for the adaptive binning (merging) algorithm.
         * @param binningBeta The beta parameter for the adaptive binning (merging) algorithm.
         * @param desiredWidth The desiredWidth for the adaptive binning (merging) algorithm.
         * 
         * @note The cost function parameters for the adaptive binning are outlined in 
         *       the AdaptBins::partialMergedCDFIntegralCost function and explained in detail
         *       in Alexander Liemen's master thesis "Adaptive Energy Binning in OPAL-X".
         */
        AdaptBins(std::shared_ptr<BunchType> bunch, BinningSelector var_selector, bin_index_type maxBins,
                  value_type binningAlpha, value_type binningBeta, value_type desiredWidth)
            : bunch_m(bunch)
            , var_selector_m(var_selector)
            , maxBins_m(maxBins)
            , binningAlpha_m(binningAlpha)
            , binningBeta_m(binningBeta)
            , desiredWidth_m(desiredWidth) {

            // Will be set during usage/rebinning
            currentBins_m = maxBins;

            initTimers();

            Inform msg("AdaptBins");
            msg << "AdaptBins initialized with maxBins = " << maxBins_m 
                << ", alpha = " << binningAlpha_m
                << ", beta = " << binningBeta_m
                << ", desiredWidth = " << desiredWidth_m << endl;
        }

        
        /**
         * @brief Initializes timers for various operations in the binning process.
         */
        void initTimers() {
            bInitLimitsT            = IpplTimings::getTimer("bInitLimits");
            bAllReduceLimitsT       = IpplTimings::getTimer("bAllReduceLimits");
            bAllReduceGlobalHistoT  = IpplTimings::getTimer("bAllReduceGlobalHisto");
            bAssignUniformBinsT     = IpplTimings::getTimer("bAssignUniformBins");
            bExecuteHistoReductionT = IpplTimings::getTimer("bExecuteHistoReduction");
            bSortContainerByBinT    = IpplTimings::getTimer("bSortContainerByBin");
            bVerifySortingT         = IpplTimings::getTimer("bVerifySorting");
        }

        /**
         * @brief Returns a view to the particle bin array.
         * 
         * @note: Change this function if the name of the Bin attribute in the container is changed.
         */
        bin_view_type getBinView() { return bunch_m->Bin.getView(); }

        /**
         * @brief Returns the current number of bins.
         * 
         * @return The current bin count.
         */
        bin_index_type getCurrentBinCount() const { return currentBins_m; }

        /**
         * @brief Gets the maximum number of bins. Will be used for the fine uniform histogram before merging.
         * 
         * @return The maximum allowed number of bins.
         */
        bin_index_type getMaxBinCount() const { return maxBins_m; }

        /**
         * @brief Returns the average binwidth.
         * 
         * @return Corresponds to (xmax_m - xmin_m)/n_bins. Used in the fine histogram.
         */
        value_type getBinWidth() const { return binWidth_m; }

        /**
         * @brief Sets the current number of bins and adjusts the bin width.
         * 
         * @param nBins The new number of currently used bins (locally).
         * 
         * @note The other parameters (limits) are set before this function is called in doFullRebin().
         */
        void setCurrentBinCount(bin_index_type nBins) {
            currentBins_m = (nBins > maxBins_m) ? maxBins_m : nBins; 
            binWidth_m    = (xMax_m - xMin_m) / currentBins_m; 
        }

        /**
         * @brief Returns the index map that sorts the particle container by bin number.
         * 
         * @return hash_type sorting the container.
         */
        hash_type getHashArray() { return sortedIndexArr_m; }

        /**
         * @brief Calculates the bin index for a given position value in a uniform histogram.
         * 
         * This static method calculates which bin a position value falls into based on the bin
         * boundaries and bin width. Uses integer casting to get bin index and clamps it to the 
         * biggest bin if out of range (only happens if limits and the particle bunch are out
         * of sync -- should not happen since particle push is not done between limit calculation
         * and bin assignment).
         * 
         * @param x The binning parameter value (e.g. z-velocity).
         * @param xMin Minimum boundary for the bins.
         * @param xMax Maximum boundary for the bins.
         * @param binWidthInv Inverse of the bin width for efficiency.
         * @param numBins The total number of bins.
         * @return The calculated bin index for the given x value.
         */
        KOKKOS_INLINE_FUNCTION
        static bin_index_type getBin(value_type x, value_type xMin, value_type xMax, value_type binWidthInv, bin_index_type numBins);

        /**
         * @brief Initializes the limits for binning based on the particle data.
         * 
         * This function calculates the minimum and maximum limits (xMin and xMax) from the
         * binning attribute (e.g. velocity), which are then used to define bin boundaries.
         * 
         * @note Called _before_ bins and histograms are initialized.
         */
        void initLimits();

        /**
         * @brief Initializes the histogram view for binning and optionally sets it to zero.
         * 
         * @param setToZero If true, initializes the histogram view to zero. Default is false. 
         *                  The 0 initialization is not needed if data is overwritten anyways.
         */
        void instantiateHistogram(bool setToZero = false);

        /**
         * @brief Assigns each particle in the bunch to a bin based on its position.
         * 
         * This function iterates over all particles in the bunch, calculates their bin
         * index using getBin(...), and updates the "Bin" attribute in the container
         * structure accordingly.
         * 
         * @note It only calculates the fine uniform histogram using maxBins bins.
         */
        void assignBinsToParticles();

        /**
         * @brief Initializes a local histogram for particle binning.
         * 
         * This function initializes the local histogram (instance of BinHisto) after bin indices
         * were assigned using assignBinsToParticles(). It reduces the indices to a histogram using
         * on the the reduction functions in this class. It chooses which function to use based on the
         * HistoReductionMode parameter (defined in ParallelReduceTools.h) or the optimal method based
         * on the current architecture/bin number.  
         */
        void initLocalHisto(HistoReductionMode modePreference = HistoReductionMode::Standard);

        /**
        * @brief Initializes and performs a team-based histogram reduction for particle bins.
        *
        * This function allocates scratch memory for each team on device, initializes a local
        * histogram for each team in shared memory, updates it based on bin indices of particles
        * assigned to that team and finally reduces the team-local histograms into a global 
        * histogram in device memory using pure atomics.
        * 
        * @details The process consists of the following steps:
        * - Allocating scratch memory for each team's local histogram.
        * - Initializing the local histogram to zero.
        * - Assigning particles to bins in parallel within each team.
        * - Reducing each team's local histogram into a global histogram (atomics).
        *
        * ### Memory and Execution
        * For GPU architecture, this function does not need to be changed. Should problems with unusal execution
        * time occur (e.g. if this function is supposed to be called on CPU -- which it shouldn't), it might make
        * sense to look at the following parameters calculated in the function:
        * - **Scratch Memory**: Scratch memory is allocated per team for a local histogram, with size `binCount`.
        * - **Concurrency**: `team_size` specifies the number of threads per team, and each team processes a `block_size`.
        *
        * @note This function is optimized for GPU execution using team-based parallelism,
        *       it does not work on Host (since `team_size` is hardcoded and to a big number). 
        *       If you want to run this on Host, change `team_size=1` and increase `block_size`.
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
        * This function performs a `Kokkos::parallel_reduce` over the particles in the bunch, incrementing 
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
        * has a `the_array` member which stores the histogram counts. So far, only `ParticleBinning::ArrayReduction`
        * is implemented in that way (to work together with `Kokkos::Sum` reducer). `the_array` needs to have a known
        * size at compile time. To use it on GPU anyways, ParallelReduceTools.h contains a function `createReductionObject()`
        * to select the reducer object from a `std::variant` of pre-compiled reducer sizes. 
        */
        template<typename ReducerType>
        void executeInitLocalHistoReduction(ReducerType& to_reduce);

        /**
         * @brief Retrieves the global histogram across all processes.
         * 
         * This function reduces the local histograms across all MPI processes into
         * a single global histogram view.
         * 
         * @return A view of the global histogram in host space (used for merging/the adaptive histogram).
         */
        void initGlobalHistogram();

        /**
         * @brief Performs a full rebinning operation on the bunch.
         * 
         * It does the following steps:
         * - Initializes the limits of the binning attribute.
         * - Sets the current number of bins to nBins (usually 128 for the fine histogram).
         * - Assigns uniform bins to particles based on the new limits and bin count.
         * - (Re-)initializes the local **and** global histogram based on the new bins.
         * 
         * @param nBins The new number of bins to use for rebinning.
         * @param recalculateLimits If true, the limits are recalculated based on the current particle data.
         * @param modePreference The preferred mode for histogram reduction (default is Standard and choses 
         *                       the best method itself).
         */
        void doFullRebin(bin_index_type nBins, bool recalculateLimits = true, HistoReductionMode modePreference = HistoReductionMode::Standard) {
            if (recalculateLimits) initLimits();
            setCurrentBinCount(nBins);
            assignBinsToParticles();
            initHistogram(modePreference);
        }

        /**
         * @brief Initializes the local/global histogram based on the Bin attribute.
         * 
         * @param modePreference The preferred mode for histogram reduction (default is Standard).
         */
        void initHistogram(HistoReductionMode modePreference = HistoReductionMode::Standard) {
            instantiateHistogram(true); 
            initLocalHisto(modePreference);
            initGlobalHistogram();

            // Init both histograms --> buils postSums and widths arrays
            localBinHisto_m.init();
            globalBinHisto_m.init();
        }

        /**
         * @brief Retrieves the number of particles in a specified bin.
         * 
         * This function returns the number of particles in the bin specified by the given index.
         * It assumes that the DualView (if used) has been properly synchronized and initialized. If the function is called
         * frequently, it might create some overhead due to the .view_host() call. However, since
         * it is only called on the host (a maximum of nBins times per iteration), the overhead
         * should be minimal. For better efficiency, one can avoid the Kokkos::View "copying-action"
         * by using dualView.h_view(binIndex). Just be careful, since this might change the underlying
         * data structure of the BinHisto class.
         * 
         * @param binIndex The index of the bin for which the number of particles is to be retrieved.
         * @param global If true, retrieves from global histogram, otherwise from local histogram.
         * @return The number of particles in the specified bin.
         */
        size_type getNPartInBin(bin_index_type binIndex, bool global = false) {
            // shouldn't happen..., "binIndex < 0" unnecessary, since binIndex is usually unsigned; but just in case the type is changed
            if (binIndex < 0 || binIndex >= getCurrentBinCount()) { return bunch_m->getTotalNum(); } 

            if (global) {
                return globalBinHisto_m.getNPartInBin(binIndex);
            } else {
                return localBinHisto_m.getNPartInBin(binIndex);
            }
        }

        /**
         * @brief Sorts the container of particles by their bin indices.
         *
         * This function assumes that the prefix sum (post-sum) has already been initialized.
         * It rearranges a particle indices array into the correct order based on their Bin attribute
         * and verifies the sorting. The function uses a parallelized bin sort.
         *
         * Steps:
         * 1. Retrieves the bin view and initializes the bin offsets using the prefix sum
         *    that is saved inside the localBinHisto_m instance (should be of type BinHisto
         *    in DualView mode). 
         * 2. Allocates a `hash_array` to store the sorted indices of particles.
         * 3. Uses a Kokkos loop to place each particle into its target position based 
         *    on its bin index.
         * 4. (In Debug mode) verifies that the sorting was successful.
         * 
         * @note This function is most efficient with a large number of bins, since this minimizes
         *       overhead from the atomic operations. Therefore, it should be called before merging
         *       (even though calling it afterwards also works). 
         */
        void sortContainerByBin();

        /**
         * @brief Returns the bin iteration policy for a given bin index.
         * 
         * This function generates a range policy for iterating over the elements within a given bin index.
         * If no DualView is used, it might need to copy some values to host, which might cause overhead.
         * However, this only happens if the localBinHisto_m is implemented as a device only view, which
         * is not its intended type.
         * 
         * @param binIndex The index of the bin for which the iteration policy is to be generated.
         * @return Kokkos::RangePolicy<> The range policy for iterating over the elements in the specified bin.
         * 
         * @note It returns an iteration policy that can be used together with sortedIndexArr_m inside
         *       `scatter()` to only iterate and scatter particles in a specific bin. 
         */
        Kokkos::RangePolicy<> getBinIterationPolicy(const bin_index_type& binIndex) {
            return localBinHisto_m.getBinIterationPolicy(binIndex);
        }

        /**
         * @brief Generates an adaptive histogram based on a fine global histogram.
         * 
         * This function is used to create an adaptive histogram from a fine global histogram. The fine global
         * histogram (usually 128 bins) has to be created beforehand, but does not necessarily have to be 
         * uniform (even though this is recommended). The algorithm works perfectly fine as long as the 
         * globalBinHisto_m.binWidths_m field is initialized correctly. If not specified otherwise, it will
         * be initialized either uniformly at initialization, or changed by the `mergeBins()` call or
         * copied from another `BinHisto` instance. 
         * 
         * It performs the following steps:
         * 1. Merges bins in the global histogram to create an adaptive histogram (see `BinHisto::mergeBins()`).
         * 2. Maps old bin indices to new bin indices using the lookup table returned by `mergeBins()`.
         * 3. Updates the local histogram with the new bin indices and widths.
         */
        void genAdaptiveHistogram();

        /**
         * @brief Prints the current global histogram to the Inform output stream.
         * 
         * This function outputs the global histogram data (bin counts) to the standard output.
         * Note: Only for rank 0 in an MPI environment.
         */
        void print() {
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
        * 
        * @note This function should not be used inside OPAL-X, since it makes assumptions regarding the actual
        *       particle bunch. It is just for testing during development. 
        */
        template <typename T, unsigned Dim>
        VField_t<T, Dim>& LTrans(VField_t<T, Dim>& field, const bin_index_type& currentBin); // TODO: may want to add usage of c constant when it exists...

    private:
        std::shared_ptr<BunchType> bunch_m;    ///< Shared pointer to the particle container.
        BinningSelector var_selector_m;        ///< Variable selector for binning.
        const bin_index_type maxBins_m;        ///< Maximum number of bins. 
        bin_index_type currentBins_m;          ///< Current number of bins in use.
        value_type xMin_m;                     ///< Minimum value of bin attribute.
        value_type xMax_m;                     ///< Maximum value of bin attribute.
        value_type binWidth_m;                 ///< Width of each bin (assumes a uniform histogram).

        // Merging cost function parameters
        value_type binningAlpha_m;             ///< Alpha parameter for binning cost function. 
        value_type binningBeta_m;              ///< Beta parameter for binning cost function.
        value_type desiredWidth_m;             ///< Desired bin width for binning cost function.

        // Histograms 
        d_histo_type localBinHisto_m;          ///< Local histogram view for bin counts.
        h_histo_type_g globalBinHisto_m;       ///< Global histogram view (over ranks reduced local histograms).

        hash_type sortedIndexArr_m;            ///< Particle index map that sorts the bunch by bin index.

        // Timers...
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


