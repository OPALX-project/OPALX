#ifndef BIN_HISTO_H
#define BIN_HISTO_H

#include "Ippl.h"
#include "BinningTools.h" // For postSum computation 

#include <iomanip>  // for std::setw, std::setprecision, etc. (debug output)

namespace ParticleBinning {

    // Define a traits class for obtaining the device view type.
    template <bool UseDualView, typename ViewType>
    struct DeviceViewTraits;

    // Specialization when UseDualView is true
    template <typename ViewType>
    struct DeviceViewTraits<true, ViewType> {
        using h_type = typename ViewType::t_host;
        using d_type = typename ViewType::t_dev;
    };

    // Specialization when UseDualView is false
    template <typename ViewType>
    struct DeviceViewTraits<false, ViewType> {
        using h_type = ViewType;
        using d_type = ViewType;
    };


    template <typename size_type, typename bin_index_type, typename value_type, 
              bool UseDualView = false, class... Properties>
    class Histogram {
    public:

        // Histogram counts view type(s)
        using view_type  = std::conditional_t<UseDualView, 
                                             Kokkos::DualView<size_type*, Properties...>, 
                                             Kokkos::View<size_type*, Properties...>>;
        using dview_type = typename DeviceViewTraits<UseDualView, view_type>::d_type;
        using hview_type = typename DeviceViewTraits<UseDualView, view_type>::h_type;        

        // Histogram widths view type(s)
        using width_view_type = std::conditional_t<UseDualView,
                                                   Kokkos::DualView<value_type*, Properties...>,
                                                   Kokkos::View<value_type*, Properties...>>;
        using hwidth_view_type = typename DeviceViewTraits<UseDualView, width_view_type>::h_type;
        using dwidth_view_type = typename DeviceViewTraits<UseDualView, width_view_type>::d_type;

        // View types to map bin indices
        template <class... Args>
        using index_transform_type  = Kokkos::View<bin_index_type*, Args...>;
        using dindex_transform_type = index_transform_type<Kokkos::DefaultExecutionSpace>; // Do not remove: needed inside AdaptBins::genAdaptiveHistogram()
        using hindex_transform_type = index_transform_type<Kokkos::HostSpace>; 
        // using hash_type             = ippl::detail::hash_type<Kokkos::DefaultExecutionSpace::memory_space>;

        /**
         * @brief Default constructor for the Histogram class.
         */
        Histogram() = default;

        /**
         * @brief Constructor for the Histogram class with a given name, number of bins, and total bin width.
         *
         * @param debug_name The name of the histogram for debugging purposes. Is passed as a name for the Kokkos::...View.
         * @param numBins The number of bins in the histogram. Might change once the merging algorithm is used.
         * @param totalBinWidth The total width of the value range covered by the particles, so $x_\mathrm{max} - x_\mathrm{min}$.
         * @param binningAlpha The alpha parameter for the adaptive binning (merging) cost function.
         * @param binningBeta The beta parameter for the adaptive binning (merging) cost function.
         * @param desiredWidth The desired width for the adaptive binning (merging) cost function.
         */
        Histogram(std::string debug_name, bin_index_type numBins, value_type totalBinWidth,
                  value_type binningAlpha, value_type binningBeta, value_type desiredWidth)
            : debug_name_m(debug_name)
            , numBins_m(numBins)
            , totalBinWidth_m(totalBinWidth)
            , binningAlpha_m(binningAlpha)
            , binningBeta_m(binningBeta)
            , desiredWidth_m(desiredWidth) {
            // Initialize the histogram, bin widths and post sum view
            instantiateHistograms();
            initTimers();
        }

        /**
         * @brief Default destructor for the Histogram class.
         */
        /*~Histogram() {
            //std::cout << "Histogram " << debug_name_m << " destroyed." << std::endl;
        } */

        /**
         * @brief Copy constructor for copying the fields from another Histogram object.
         * 
         * @see copyFields() for more information on how the fields are copied (uses Kokkos shallow copy).
         */
        Histogram(const Histogram& other) {
            copyFields(other);
        }

        /**
         * @brief Assignment operator for copying the fields from another Histogram object.
         * 
         * @see copyFields() for more information on how the fields are copied (uses Kokkos shallow copy).
         */
        Histogram& operator=(const Histogram& other) {
            if (this == &other) return *this;
            copyFields(other);
            return *this;
        }

        /**
         * @brief Initializes timers for various operations in the binning process.
         *
         * Timers initialized:
         * - bDeviceSyncronizationT: Measures the time taken for device synchronization.
         * - bHistogramInitT: Measures the time taken to initialize histograms.
         * - bMergeBinsT: Measures the time taken to merge bins.
         */
        void initTimers() {
            bDeviceSyncronizationT = IpplTimings::getTimer("bDeviceSyncronization");
            bHistogramInitT        = IpplTimings::getTimer("bHistogramInit");
            bMergeBinsT            = IpplTimings::getTimer("bMergeBins");
        }

        /**
         * @brief Retrieves the number of particles in a specified bin.
         *
         * This function returns the number of particles in the bin specified by the given index.
         * It assumes that the DualView has been properly synchronized and initialized. If the function is called
         * frequently, it might create some overhead due to the .view_host() call. However, since
         * it is only called on the host (a maximum of nBins times per iteration), the overhead
         * should be minimal. For better efficiency, one can avoid the Kokkos::View "copying-action"
         * by using dualView.h_view(binIndex).
         *
         * @tparam UseDualView A boolean template parameter indicating whether DualView is used.
         * @param binIndex The index of the bin for which the number of particles is to be retrieved.
         * @return The number of particles in the specified bin.
         */
        size_type getNPartInBin(bin_index_type binIndex) {
            if constexpr (UseDualView) {
                return histogram_m.h_view(binIndex);
            } else if (std::is_same<typename hview_type::memory_space, Kokkos::HostSpace>::value) {
                return histogram_m(binIndex);
            } else {
                std::cerr << "Warning: Accessing BinHisto.getNPartInBin without DualView might be inefficient!" << std::endl;
                Kokkos::View<size_type, Kokkos::HostSpace> host_scalar("host_scalar");
                Kokkos::deep_copy(host_scalar, Kokkos::subview(histogram_m, binIndex));
                return host_scalar();
            }
        }

        size_type getCurrentBinCount() const { return numBins_m; }

        view_type getHistogram() { return histogram_m; }

        view_type getPostSum() { return postSum_m; }      
        
        width_view_type getBinWidths() const { return binWidths_m; }

        /**
         * @brief Sets the bin widths by copying them from a different Histogram instance (usually neccessary after
         *        merging the global histogram to copy the new non-uniform widths to the local histogram).
         * 
         * @note This function needs to be templated to allow copying from BinHisto instances that are from a 
         *       different type than "this". It is neccessary in this context, since localBinHisto_m in AdaptBins is
         *       usually a DualView setup, while globalBinHisto_m is modified to use a host-only view.
         */
        template <typename Histogram_t>
        void copyBinWidths(const Histogram_t& other) {
            using other_dwidth_view_type = typename Histogram_t::dwidth_view_type;
            Kokkos::deep_copy(getDeviceView<dwidth_view_type>(binWidths_m), 
                              other.template getDeviceView<other_dwidth_view_type>(other.getBinWidths()));
            if constexpr (UseDualView) {
                binWidths_m.modify_device();

                IpplTimings::startTimer(bDeviceSyncronizationT);
                binWidths_m.sync_host();
                IpplTimings::stopTimer(bDeviceSyncronizationT);
            }
        }

        
        /**
         * @brief Instantiates the histogram, bin widths, and post-sum views (Possibly DualView).
         */
        void instantiateHistograms() {
            histogram_m = view_type("histogram", numBins_m);
            binWidths_m = width_view_type("binWidths", numBins_m);
            postSum_m   = view_type("postSum", numBins_m + 1);
        }


        /**
         * @brief Synchronizes the histogram view and initializes the bin widths and post-sum.
         *
         * @note The bin widths are assumed to be constant. Should only be called the first time the histogram is created.
         */
        void init() { // const value_type constBinWidth
            // Assumes you have initialized histogram_m from the outside!
            sync();
            initConstBinWidths(totalBinWidth_m);
            initPostSum();
        }


        /**
         * @brief Initializes the bin widths with a constant value.
         *
         * @param constBinWidth The constant value to set for all bin widths.
         * 
         * @note Should not be called again after merging bins, since the bin widths won't be uniform.
         */
        void initConstBinWidths(const value_type constBinWidth) {
            dwidth_view_type dWidthView = getDeviceView<dwidth_view_type>(binWidths_m);
            const value_type binWidth   = constBinWidth / numBins_m;
            using execution_space       = typename dwidth_view_type::execution_space;

            IpplTimings::startTimer(bHistogramInitT);
            // Note: "Kokkos::deep_copy(getDeviceView<dwidth_view_type>(binWidths_m), constBinWidth / numBins_m);" would work too!
            Kokkos::parallel_for("InitConstBinWidths", 
                Kokkos::RangePolicy<execution_space>(0, numBins_m), KOKKOS_LAMBDA(const size_t i) {
                    dWidthView(i) = binWidth;
                }
            );
            IpplTimings::stopTimer(bHistogramInitT);

            if constexpr (UseDualView) {
                binWidths_m.modify_device();
                IpplTimings::startTimer(bDeviceSyncronizationT);
                binWidths_m.sync_host();
                IpplTimings::stopTimer(bDeviceSyncronizationT);
            }
        }


        /**
         * @brief Initializes and computes the post-sum for the histogram.
         *
         * This function initializes the post-sum by computing the fixed sum on the device view of the post-sum member.
         * If the UseDualView constant is true, it modifies the device view and synchronizes it with the host.
         */
        void initPostSum() {
            IpplTimings::startTimer(bHistogramInitT);
            computeFixSum<dview_type>(getDeviceView<dview_type>(histogram_m), getDeviceView<dview_type>(postSum_m));
            IpplTimings::stopTimer(bHistogramInitT);
            
            if constexpr (UseDualView) {
                postSum_m.modify_device();
                IpplTimings::startTimer(bDeviceSyncronizationT);
                postSum_m.sync_host();
                IpplTimings::stopTimer(bDeviceSyncronizationT);
            }
        }


        /**
         * @brief Returns a `Kokkos::RangePolicy` for iterating over the elements in a specified bin.
         *
         * This function generates a range policy for iterating over the elements within a given bin 
         * index range. 
         *
         * @tparam bin_index_type The type of the bin index.
         * @param binIndex1 The index of the bin for which the iteration policy is to be generated.
         * @param numBins The number of bins to iterate over (default is 1) starting at `binIndex1`.
         * @return Kokkos::RangePolicy<> The range policy for iterating over the elements in the specified bin.
         * 
         * @note The bounds are needed on host. So if no `DualView` is used or the values are only available on 
         *       device, then there will be some values copied from device to host.
         */
        Kokkos::RangePolicy<> getBinIterationPolicy(const bin_index_type& binIndex1, const bin_index_type numBins = 1) {
            if constexpr (UseDualView) {
                return Kokkos::RangePolicy<>(postSum_m.h_view(binIndex1), postSum_m.h_view(binIndex1 + numBins));
            } else {
                std::cerr << "Warning: Accessing BinHisto.getBinIterationPolicy without DualView might be inefficient!" << std::endl;
                Kokkos::View<bin_index_type[2], Kokkos::HostSpace> host_ranges("host_scalar");
                Kokkos::deep_copy(host_ranges, Kokkos::subview(postSum_m, std::make_pair(binIndex1, binIndex1 + numBins)));
                return Kokkos::RangePolicy<>(host_ranges(0), host_ranges(1));
            }
        }


        /*
        Below are methods used for syncing the histogram view between host and device.
        If a normal View is used, they have no effect.
        Only necessary for the histogram view, since the binWidths and postSum views 
        are only modified inside this class. 
        */

        /**
         * @brief Synchronizes the histogram data between host and device.
         *
         * This function checks if the histogram data needs to be synchronized between
         * the host and the device. If both the host and device have modifications, it
         * issues a warning and overwrites the changes on the host. It then performs
         * the necessary synchronization based on where the modifications occurred.
         *
         * @note This function only performs synchronization if the `UseDualView` 
         *       template parameter is true. Otherwise it does nothing.
         */
        void sync() {
            IpplTimings::startTimer(bDeviceSyncronizationT);
            if constexpr (UseDualView) {
                if (histogram_m.need_sync_host() && histogram_m.need_sync_device()) {
                    std::cerr << "Warning: Histogram was modified on host AND device -- overwriting changes on host." << std::endl;
                } 
                if (histogram_m.need_sync_host()) {
                    histogram_m.sync_host();
                } else if (histogram_m.need_sync_device()) {
                    histogram_m.sync_device();
                } // else do nothing
            }
            IpplTimings::stopTimer(bDeviceSyncronizationT);
        }


        /**
         * @brief If a DualView is used, it sets the flag on the view that the device has been modified.
         * 
         * @see sync() After this function you might call sync at some point.
         */
        void modify_device() { if constexpr (UseDualView) histogram_m.modify_device(); }


        /**
         * @brief If a DualView is used, it sets the flag on the view that the host has been modified.
         * 
         * @see sync() After this function you might call sync at some point.
         */
        void modify_host() { if constexpr (UseDualView) histogram_m.modify_host(); }


        /**
         * @brief Retrieves the device view of the histogram.
         *
         * This function returns the device view of the histogram if the `UseDualView`
         * flag is set to true. Otherwise, it returns the histogram itself (which might not be on device!).
         *
         * @tparam HistogramType The type of the histogram.
         * @param histo Reference to the histogram.
         * @return The device view of the histogram if `UseDualView` is true, otherwise the histogram itself.
         */
        template <typename return_type, typename HistogramType>
        static constexpr return_type getDeviceView(HistogramType histo) {
            if constexpr (UseDualView) {
                return histo.view_device();
            } else {
                return histo;
            }
        }


        /**
         * @brief Retrieves a host view of the given histogram.
         *
         * This function returns a host view of the provided histogram object. If a 
         * DualView is used, it calls the `view_host()`, otherwise it returns the normal view.
         *
         * @tparam HistogramType The type of the histogram object.
         * @param histo The histogram object from which to retrieve the host view.
         * @return A host view of the histogram object.
         */
        template <typename return_type, typename HistogramType>
        static constexpr return_type getHostView(HistogramType histo) {
            if constexpr (UseDualView) {
                return histo.view_host();
            } else {
                return histo;
            }
        }

        
        /**
         * @brief Merges bins in a histogram to reduce the number of bins while minimizing a cost function.
         *
         * This function performs adaptive binning by merging adjacent bins in a histogram. The merging process
         * is guided by a cost function that evaluates the deviation introduced by merging bins. The goal is to
         * minimize the total cost while reducing the number of bins.
         *
         * @return A Kokkos View mapping old bin indices to new bin indices after merging.
         *
         * @details
         * The function performs the following steps:
         * 1. Ensures the histogram is available on the host (if using DualView).
         * 2. Computes prefix sums for bin counts and widths on the host.
         * 3. Uses dynamic programming to find the optimal bin merging strategy that minimizes the total cost.
         * 4. Reconstructs the boundaries of the merged bins based on the dynamic programming results.
         * 5. Builds new arrays for the merged bins, including counts and widths.
         * 6. Generates a lookup table mapping old bin indices to new bin indices.
         * 7. Overwrites the old histogram arrays with the new merged ones.
         * 8. If using DualView, marks the host data as modified and synchronizes with the device.
         * 9. Recomputes the post-sum for the new merged histogram.
         *
         * The cost function used for merging is defined by `adaptiveBinningCostFunction()`, which evaluates
         * the deviation introduced by merging bins based on their counts and widths.
         *
         * @note This function assumes that the histogram has been generated with a sufficiently large number
         * of bins (e.g., 128 bins). If the number of bins is too small (<= 1), merging is skipped.
         *
         * @warning If no feasible merges are found, the function falls back to a no-merge scenario and logs
         * a warning.
         */
        hindex_transform_type mergeBins();

        /**
         * @brief Computes the cost function for adaptive binning in a histogram.
         *
         * This function calculates a cost value based on several factors, including 
         * Shannon entropy, bin width penalties, and sparsity penalties. It is used 
         * to optimize the binning process by balancing the distribution of particles 
         * across bins while considering desired bin width and sparsity constraints.
         *
         * @param sumCount The total count of particles in the current to-be-merged bin.
         * @param sumWidth The total width of the current to-be-merged bin.
         * @param totalNumParticles The total number of particles across all bins.
         * 
         * @return The computed cost value for the current to-be-merged bin configuration.
         *
         * @note The function assumes `sumCount > 0` as a precondition. This is only checked
         *       in debug mode. However, the cost should not even be computed in that case.
         *
         * @details
         * The cost function is composed of the following terms:
         * - **Shannon Entropy**: Encourages a balanced distribution of particles.
         * - **Wide Bin Penalty**: Penalizes bins that are too wide or too narrow, 
         *   based on the `binningAlpha` parameter.
         * - **Bin Size Bias**: Biases the bin size towards `desiredWidth`, controlled 
         *   by the `binningBeta` parameter.
         * - **Sparse Penalty**: Penalizes bins with too few particles, particularly 
         *   in the distribution tails, normalized by the desired width.
         * 
         * The cost function is defined as:
         * \f[
         * c[i, j] \overset{!}{=} 
         * \underbrace{N_\mathrm{tot}\, \log\left(N_\mathrm{tot}\right)\, \Delta v_\mathrm{tot}}_{\text{Shannon Entropy}}
         * + \underbrace{\alpha \, \Delta v_\mathrm{tot}^2}_{\text{Width Bias}}
         * + \underbrace{\beta \, \left(\Delta v_\mathrm{tot} - \Delta v_\mathrm{bias}\right)^2}_{\text{Deviation Penalty}}
         * + \underbrace{\frac{\Delta v_\mathrm{bias}}{N_\mathrm{tot}} \Theta\left(\Delta v_\mathrm{bias}-N_\mathrm{tot}\right)}_{\text{Sparse Penalty}}
         * \f]
         * 
         * Where:
         * - \( N_\mathrm{tot} \): Total number of particles in the bin.
         * - \( \Delta v_\mathrm{tot} \): Total width of the bin.
         * - \( \Delta v_\mathrm{bias} \): Desired bin width.
         * - \( \alpha \): Weight for the width bias term.
         * - \( \beta \): Weight for the deviation penalty term.
         * - \( \Theta(x) \): Heaviside step function.
         * 
         * The function is designed to judge how exactly an "optimally merged histogram"
         * looks like. For a more detailed explanation, see Alexander Liemen's master thesis
         * ("Adaptive Energy Binning in OPAL-X").
         */
        value_type adaptiveBinningCostFunction(
            const size_type& sumCount,
            const value_type& sumWidth,
            const size_type& totalNumParticles
        );

    private:
        std::string debug_name_m;    /// \brief Debug name for identifying the histogram instance.
        bin_index_type numBins_m;    /// \brief Number of bins in the histogram.
        value_type totalBinWidth_m;  /// \brief Total width of all bins combined.

        value_type binningAlpha_m;   /// \brief Alpha parameter for the adaptive binning (merging) cost function.
        value_type binningBeta_m;    /// \brief Beta parameter for the adaptive binning (merging) cost function.
        value_type desiredWidth_m;   /// \brief Desired width for the adaptive binning (merging) cost function.

        view_type       histogram_m; /// \brief View storing the particle counts in each bin.
        width_view_type binWidths_m; /// \brief View storing the widths of the bins.
        view_type       postSum_m;   /// \brief View storing the cumulative sum of bin counts (used in sorting, generating range policies).

        IpplTimings::TimerRef bDeviceSyncronizationT;
        IpplTimings::TimerRef bHistogramInitT;
        IpplTimings::TimerRef bMergeBinsT;

        /**
         * @brief Copies the fields from another Histogram object.
         *
         * This function copies the internal fields from the provided Histogram object
         * to the current object. The fields are copied using Kokkos' shallow copy.
         *
         * @param other The Histogram object from which to copy the fields.
         */
        void copyFields(const Histogram& other);


    /*
    Here are just some debug functions, like a nice histogram output formatted as python numpy arrays.
    */
    public:
        /**
         * @brief Prints a nicely formatted table of bin indices, counts, and widths.
         *
         * @param os The output stream to write to (defaults to std::cout).
         */
        void printHistogram(std::ostream &os = std::cout) {
            if (ippl::Comm->rank() != 0) return;
            hview_type countsHost       = getHostView<hview_type>(histogram_m); 
            hwidth_view_type widthsHost = getHostView<hwidth_view_type>(binWidths_m);

            // 3) Print header
            os << "Histogram \"" << debug_name_m << "\" with " << numBins_m << " bins. BinWidth = " << totalBinWidth_m << ".\n\n";

            // Format columns: BinIndex, Count, Width
            // Adjust widths as needed
            os << std::left 
            << std::setw(10) << "Bin" 
            << std::right 
            << std::setw(12) << "Count"
            << std::setw(16) << "Width\n";

            os << std::string(38, '-') << "\n"; 
            // (38 dashes or however many you prefer to underline)

            // 4) Print each bin
            for (bin_index_type i = 0; i < numBins_m; ++i) {
                os << std::left << std::setw(10) << i    // bin index left-aligned
                << std::right << std::setw(12) << countsHost(i)
                << std::fixed << std::setw(16) << std::setprecision(6) 
                << static_cast<double>(widthsHost(i))  // in case 'value_type' is double/float
                << "\n";
            }

            //os << "-----------------------------------------" << endl;
            os << std::endl; // extra newline at the end
        }

        void printPythonArrays() const {
            if (ippl::Comm->rank() != 0) return;
            hview_type hostCounts = getHostView<hview_type>(histogram_m);
            hwidth_view_type hostWidths = getHostView<hwidth_view_type>(binWidths_m);
            // TODO: if I leave this here, it may need a deep_copy to make it save for every execution space

            // Output counts as a Python NumPy array
            std::cout << "bin_counts = np.array([";
            for (bin_index_type i = 0; i < numBins_m; ++i) {
                std::cout << hostCounts(i);
                if (i < numBins_m - 1) std::cout << ", ";
            }
            std::cout << "])" << std::endl;

            // Output widths as a Python NumPy array
            std::cout << "bin_widths = np.array([";
            for (bin_index_type i = 0; i < numBins_m; ++i) {
                std::cout << std::fixed << std::setprecision(6) << hostWidths(i);
                if (i < numBins_m - 1) std::cout << ", ";
            }
            std::cout << "])" << std::endl;
        }

    };
    
}

#include "BinHisto.tpp"

#endif // BIN_HISTO_H
