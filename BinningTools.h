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
        if (modePreference == HistoReductionMode::Standard) {
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

        const int axis; ///< Index of the coordinate axis to use for binning.
        position_view_type data_arr; ///< Kokkos view of the particle data array.

        /**
        * @brief Constructs a CoordinateSelector for a specific axis.
        *
        * @param data_ Kokkos view representing the particle data array.
        * @param axis_ Index of the axis to use for binning (e.g., 0 for x, 1 for y, 2 for z).
        */
        CoordinateSelector(position_view_type data_, int axis_) 
            : data_arr(data_)
            , axis(axis_) {}

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
            return data_arr(i)[axis];
        }
    };



} // namespace ParticleBinning

#endif // BINNINGTOOLS_H