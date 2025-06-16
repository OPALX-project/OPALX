#ifndef BINNINGTOOLS_H
#define BINNINGTOOLS_H

#include "ParallelReduceTools.h" // needed for HistoReductionMode and maxArrSize

namespace ParticleBinning {

    /**
    * @brief Determines the appropriate histogram reduction mode based on user preference, 
    *        bin count, and execution environment.
    *
    * This function selects the optimal histogram reduction method. If the code is compiled 
    * with a host execution space (e.g., Serial or OpenMP), it forces the reduction mode to 
    * `HistoReductionMode::HostOnly`, disregarding the user's preference. Otherwise, if the 
    * user preference is `HistoReductionMode::Standard`, it automatically chooses between 
    * `ParallelReduce` or `TeamBased` based on the `binCount`. If a specific preference is 
    * provided (not `Standard`), that preference is respected.
    *
    * @param modePreference The user's preferred reduction mode.
    *                       - `Standard` to select a mode based on bin count.
    *                       - `ParallelReduce`, `TeamBased`, or `HostOnly` to force a specific mode.
    * @return HistoReductionMode The selected histogram reduction mode:
    *         - `HostOnly` if the default execution space is a host space.
    *         - `ParallelReduce` if `binCount` is within `maxArrSize<bin_index_type>`.
    *         - `TeamBased` if `binCount` exceeds `maxArrSize<bin_index_type>`.
    *         - Otherwise, respects the specified `modePreference`.
    *
    * @note If compiled for a host-only execution environment, the returned mode will always be 
    *       `HistoReductionMode::HostOnly`, regardless of `modePreference`.
    * @see HistoReductionMode
    */
    template <typename bin_index_type>
    HistoReductionMode determineHistoReductionMode(HistoReductionMode modePreference, bin_index_type binCount) {
        // Overwrite standard mode if compiled with default host execution space!
        if (std::is_same<Kokkos::DefaultExecutionSpace, Kokkos::DefaultHostExecutionSpace>::value) return HistoReductionMode::HostOnly;

        // Otherwise choose automatically if Standard and respect preference if not on host and not standard!
        if (modePreference == HistoReductionMode::Standard) { //  || modePreference == HistoReductionMode::HostOnly
            //if (modePreference == HistoReductionMode::HostOnly) {
            //    std::cerr << "Warning: HostOnly mode is not supported on CUDA! Switching to Standard mode." << std::endl;
            //}
            return (binCount <= maxArrSize<bin_index_type>) ? HistoReductionMode::ParallelReduce : HistoReductionMode::TeamBased;
        } else {
            return modePreference;
        }
    } 


    /**
    * @brief Example struct used to access the binning variable for each particle.
    *
    * This struct provides a flexible way to select a specific coordinate axis
    * (e.g., position or velocity) for binning. It allows specifying the axis
    * index, making it versatile for different coordinate-based binning operations.
    * Alternatively, custom structs could be defined for other binning variables, such
    * as calculating energy with saved mass and velocity values.
    *
    * Requirement:
    * - `operator()` must be defined with the exact signature shown below.
    * - `value_type` must be defined in the struct.
    * - If `operator()` does not return `value_type`, `AdaptBins` will not compile,
    *   and may cause segmentation faults.
    *
    * Example usage:
    * ```
    * CoordinateSelector<bunch_type> selector(bunch->R.getView(), 2); // Access z-axis of position
    * ```
    *
    * @note Using a shared pointer does not work, since these are managed on host only! Using a simple
    *       passed Kokkos::View also does not work, since create() might be called which calls Kokkos::realloc(),
    *       which is a problem, since the copied view does not know the reallocated data. 
    *       Therefore, the view is updated with `updateDataArr`!
    *
    * @tparam bunch_type Type of the particle bunch (derived from ParticleBase).
    */
    template <typename bunch_type>
    struct CoordinateSelector {
        /// Type representing the value of the binning variable (e.g., position or velocity).
        using value_type = typename bunch_type::Layout_t::value_type;

        /// Type representing the size of the particle bunch.
        using size_type = typename bunch_type::size_type;

        /// Type representing the view of particle positions.
        using position_view_type = typename bunch_type::particle_position_type::view_type;

        position_view_type data_arr; ///< Kokkos view of the particle data array.
        const int axis;              ///< Index of the coordinate axis to use for binning.

        /**
        * @brief Constructs a CoordinateSelector for a specific axis.
        *
        * @param data_ Kokkos view representing the particle data array.
        * @param axis_ Index of the axis to use for binning (e.g., 0 for x, 1 for y, 2 for z).
        */
        CoordinateSelector(int axis_) 
            : axis(axis_) {}

        /**
        * @brief Updates the data array view with the latest particle data.
        *
        * This function updates `data_arr` to reflect the latest particle data in the container
        * by retrieving the view from `bunch->R`. This ensures `data_arr` is synchronized with any
        * recent changes to the particle data (if bunch->create() is called between binnings!).
        * AdaptBins calls this function automatically everytime the binning parameter is used.
        *
        * @note If you have your own class, you need to implement this signature and potentially
        *       update all data arrays that operator() uses (e.g. also velocity or E).
        */
        void updateDataArr(std::shared_ptr<bunch_type> bunch) { data_arr = bunch->P.getView(); }

        /**
        * @brief Returns the value of the binning variable for a given particle index.
        *
        * This function is called by AdaptBins to obtain the binning value for each particle.
        *
        * @param i Index of the particle in the data array.
        * @return value_type Value of the specified coordinate axis for the particle at index `i`.
        */
        KOKKOS_INLINE_FUNCTION
        value_type operator()(const size_type& i) const {
            //std::cout << "CoordinateSelector: " << i << std::endl; // TODO: debug, remove later
            const value_type value = fabs(data_arr(i)[axis]);
            return value / sqrt(1 + value * value); // Normalize to v/c, so v in [0, 1]
        }
    };


    /**
     * @brief Computes the post- or prefix-sum of the input view and stores the result in the ...-sum view.
     *
     * @tparam SizeType The type of the elements in the input and ...-sum views.
     * @param input_view The input view containing the elements to be summed.
     * @param post_sum_view The output view where the ...-sum results will be stored. It must have a size of input_view.extent(0) + 1.
     * @param postSum If true, the post-sum is computed; otherwise, the prefix-sum is computed.
     *
     * @throws `ippl::Comm->abort();` if the size of the post_sum_view is not equal to input_view.extent(0) + 1.
     */
    template <typename ViewType>
    void computeFixSum(const ViewType& input_view, const ViewType& post_sum_view) {
        using execution_space = typename ViewType::execution_space;
        using size_type       = typename ViewType::size_type;
        using value_type      = typename ViewType::value_type;

        // Ensure the output view has the correct size
        if (post_sum_view.extent(0) != input_view.extent(0) + 1) {
            Inform m("computePostSum");
            m << "Output view must have size input_view.extent(0) + 1" << endl;
            ippl::Comm->abort();
        }

        // Initialize the first element to 0
        Kokkos::parallel_for("InitPostSum", Kokkos::RangePolicy<execution_space>(0, 1),
            KOKKOS_LAMBDA(const size_type) {
                post_sum_view(0) = 0;
            });

        // Compute the fix sum
        Kokkos::parallel_scan("ComputePostSum", Kokkos::RangePolicy<execution_space>(0, input_view.extent(0)),
            KOKKOS_LAMBDA(const size_type& i, value_type& partial_sum, bool final) {
                partial_sum += input_view(i);
                if (final) {
                    post_sum_view(i + 1) = partial_sum;
                }
            });
    }


    /**
     * @brief Checks if the elements in a Kokkos::View are sorted in non-decreasing order.
     *
     * @tparam ValueType The type of the elements in the Kokkos::View so be checked.
     * @tparam SizeType The type used for indexing and size.
     * @param view The Kokkos::View containing the elements to be checked.
     * @param indices Argsort of view.
     * @param npart The number of elements in the view to be checked.
     * @return true if the elements are sorted in non-decreasing order, false otherwise.
     * 
     * @note This function does not leverage short circuiting.
     */
    template <typename ValueType, typename SizeType, typename HashType>
    bool viewIsSorted (const Kokkos::View<ValueType*> view, HashType indices, SizeType npart) {
        bool sorted = true;
        Kokkos::parallel_reduce("CheckSorted", npart - 1, KOKKOS_LAMBDA(const SizeType& i, bool& update) {
            if (view(indices(i)) > view(indices(i + 1))) update = false;
        }, Kokkos::LAnd<bool>(sorted));
        return sorted;
    }


} // namespace ParticleBinning

#endif // BINNINGTOOLS_H