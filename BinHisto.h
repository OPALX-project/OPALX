#ifndef BIN_HISTO_H
#define BIN_HISTO_H

#include "Ippl.h"
#include "BinningTools.h" // For postSum computation 

#include <iomanip>  // for std::setw, std::setprecision, etc. (debug output)

namespace ParticleBinning {

    /**
     * \enum HistoTypeIdentifier
     * \brief Enum to identify the type of histogram.
     * \details Provides identifiers for different histogram-related views.
     *
     * \var HistoTypeIdentifier::Histogram
     * Indicates the main histogram view, storing the particle counts in bins.
     *
     * \var HistoTypeIdentifier::BinWidth
     * Indicates the view storing the widths of each bin.
     *
     * \var HistoTypeIdentifier::PostSum
     * Indicates the view storing the cumulative sums for post-processing.
     */
    /*enum class HistoTypeIdentifier {
        Histogram, ///< Main histogram view.
        BinWidth,  ///< Bin width view.
        PostSum    ///< Post-sum view.
    };*/

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

        using view_type  = std::conditional_t<UseDualView, 
                                             Kokkos::DualView<size_type*, Properties...>, 
                                             Kokkos::View<size_type*, Properties...>>;
        /*using dview_type = std::conditional_t<UseDualView, 
                                              typename view_type::t_dev, 
                                              view_type>;
        using hview_type = std::conditional_t<UseDualView, 
                                              typename view_type::t_host, 
                                              view_type>;*/
        using dview_type = typename DeviceViewTraits<UseDualView, view_type>::d_type;
        using hview_type = typename DeviceViewTraits<UseDualView, view_type>::h_type;        

        using width_view_type = std::conditional_t<UseDualView,
                                                   Kokkos::DualView<value_type*, Properties...>,
                                                   Kokkos::View<value_type*, Properties...>>;
        /*using hwidth_view_type = std::conditional_t<UseDualView, 
                                                   typename width_view_type::t_host, 
                                                   width_view_type>;
        using dwidth_view_type = std::conditional_t<UseDualView, 
                                                    typename width_view_type::t_dev, 
                                                    width_view_type>;*/
        using hwidth_view_type = typename DeviceViewTraits<UseDualView, width_view_type>::h_type;
        using dwidth_view_type = typename DeviceViewTraits<UseDualView, width_view_type>::d_type;

        template <class... Args>
        using index_transform_type = Kokkos::View<bin_index_type*, Args...>;
        
        using dindex_transform_type = index_transform_type<Kokkos::DefaultExecutionSpace>;
        using hindex_transform_type = index_transform_type<Kokkos::HostSpace>;

        using hash_type             = ippl::detail::hash_type<Kokkos::DefaultExecutionSpace::memory_space>;

        /**
         * @brief Default constructor for the Histogram class.
         */
        Histogram() = default;

        /**
         * @brief Constructor for the Histogram class with a given name, number of bins, and total bin width.
         *
         * @param debug_name The name of the histogram for debugging purposes. Is passed as a name for the Kokkos::...View.
         * @param numBins The number of bins in the histogram. Might change once the adaptive algorithm is used.
         * @param totalBinWidth The total width of the value range covered by the particles, so $x_\mathrm{max} - x_\mathrm{min}$.
         */
        Histogram(std::string debug_name, bin_index_type numBins, value_type totalBinWidth,
                  value_type binningAlpha, value_type binningBeta, value_type desiredWidth)
            : debug_name_m(debug_name)
            , numBins_m(numBins)
            , totalBinWidth_m(totalBinWidth)
            , binningAlpha_m(binningAlpha)
            , binningBeta_m(binningBeta)
            , desiredWidth_m(desiredWidth) {
            
            // binningAlpha_m = Options::binningAlpha;
            // binningBeta_m  = Options::binningBeta;
            // desiredWidth_m = Options::desiredWidth;
            instantiateHistograms();

            // Initialize timers
            initTimers();
        }

        /**
         * @brief Default destructor for the Histogram class.
         */
        ~Histogram() {
            //std::cout << "Histogram " << debug_name_m << " destroyed." << std::endl;
        } 

        /**
         * @brief Copy constructor for copying the fields from another Histogram object.
         * 
         * @see copyFields() for more information on how the fields are copied (Kokkos shallow copy).
         */
        Histogram(const Histogram& other) {
            copyFields(other);
        }

        /**
         * @brief Assignment operator for copying the fields from another Histogram object.
         * 
         * @see copyFields() for more information on how the fields are copied (Kokkos shallow copy).
         */
        Histogram& operator=(const Histogram& other) {
            if (this == &other) return *this;
            copyFields(other);
            return *this;
        }

        void initTimers() {
            bDeviceSyncronizationT = IpplTimings::getTimer("bDeviceSyncronization");
            bHistogramInitT        = IpplTimings::getTimer("bHistogramInit");
            bMergeBinsT            = IpplTimings::getTimer("bMergeBins");
        }
        
        /*
        Some functions to access only single elements of the histogram.
        */

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
            } else if (std::is_same<typename hview_type::memory_space, Kokkos::HostSpace>::value) { // No DualView, but on host anyways
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
         * @brief Sets the bin widths by copying them from a different Histogram instance (usually global -> local)
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

        /*
        Some function for initialization.
        */
        
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
            //static IpplTimings::TimerRef histoInitTimer = IpplTimings::getTimer("syncInitHistoTools");

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
         * @note Should not be called again after merging bins, since the bin widths will all be different.
         */
        void initConstBinWidths(const value_type constBinWidth) {
            dwidth_view_type dWidthView = getDeviceView<dwidth_view_type>(binWidths_m);
            const value_type binWidth   = constBinWidth / numBins_m;
            using execution_space       = typename dwidth_view_type::execution_space;
            //Kokkos::deep_copy(dWidthView, binWidth);
            IpplTimings::startTimer(bHistogramInitT);
            Kokkos::parallel_for("InitConstBinWidths", 
                Kokkos::RangePolicy<execution_space>(0, numBins_m), KOKKOS_LAMBDA(const size_t i) {
                    dWidthView(i) = binWidth;
                }
            );
            //std::cout << "Hey!" << std::endl;
            IpplTimings::stopTimer(bHistogramInitT);
            /*
            Note: DON'T use "Kokkos::deep_copy(getDeviceView<dwidth_view_type>(binWidths_m), constBinWidth / numBins_m);"!
            For some reason, this resulted in a huge overhead (always 0.3s just for this function)
            */
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
         *
         * @tparam size_type The type used for the size of the elements.
         */
        void initPostSum() {
            //auto postSumView = constexpr UseDualView ? postSum_m.view_device() : postSum_m;
            // dview_type postSumView = getDeviceView(postSum_m);
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
         * @brief Returns a Kokkos::RangePolicy for iterating over the elements in a specified bin.
         *
         * This function generates a range policy for iterating over the elements within a given bin index.
         * If no DualView is used, it needs to copy some values to host, which might cause overhead.
         *
         * @tparam bin_index_type The type of the bin index.
         * @param binIndex1 The index of the bin for which the iteration policy is to be generated.
         * @param numBins The number of bins to iterate over (default is 1) starting at `binIndex1`.
         * @return Kokkos::RangePolicy<> The range policy for iterating over the elements in the specified bin.
         */
        Kokkos::RangePolicy<> getBinIterationPolicy(const bin_index_type& binIndex1, const bin_index_type numBins = 1) {
            if constexpr (UseDualView) {
                // localPostSumHost = postSum_m.view_host();
                //hview_type localPostSumHost = getHostView<hview_type>(postSum_m);
                //return Kokkos::RangePolicy<>(localPostSumHost(binIndex), localPostSumHost(binIndex + 1));
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
        Will only do this for the histogram view, since the binWidths and postSum views 
        are modified inside this class only. 
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
            //static IpplTimings::TimerRef histoSyncOperation = IpplTimings::getTimer("histoSyncOperation");
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
            //IpplTimings::stopTimer(histoSyncOperation);
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
         * flag is set to true. Otherwise, it returns the histogram itself.
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

        /*
        The following contain functions that are used to make the histogram adaptive.
        */

        //template <typename BinningSelector_t>
        hindex_transform_type mergeBins(/*const hash_type sortedIndexArr, const BinningSelector_t var_selector*/);

        /*KOKKOS_INLINE_FUNCTION // in case it is needed...
        static value_type computeDeviationCost(const size_type& sumCount,
                                              const value_type& sumWidth,
                                              const value_type& maxBinRatio,
                                              const value_type& alpha, 
                                              const value_type& mergedStd);*/

        /*value_type mergedBinStd(
            const bin_index_type& i, const bin_index_type& k,
            const size_type& sumCount, const value_type& varPerBin, 
            const hwidth_view_type& prefixWidth, 
            const hview_type& fineCounts, const hwidth_view_type& fineWidths
        );*/

        //template <typename BinningSelector_t>
        value_type partialMergedCDFIntegralCost(
            //const bin_index_type& i, const bin_index_type& k,
            const size_type& sumCount,
            const value_type& sumWidth,
            const size_type& totalNumParticles
            //const hash_type sortedIndexArr, 
            //const BinningSelector_t var_selector
        );

    private:
        std::string debug_name_m;   /// \brief Debug name for identifying the histogram instance.
        bin_index_type numBins_m;   /// \brief Number of bins in the histogram.
        value_type totalBinWidth_m; /// \brief Total width of all bins combined.

        value_type binningAlpha_m; 
        value_type binningBeta_m;
        value_type desiredWidth_m;

        view_type       histogram_m;      /// \brief View storing the particle counts in each bin.
        width_view_type binWidths_m;      /// \brief View storing the widths of the bins.
        view_type       postSum_m;        /// \brief View storing the cumulative sum of bin counts (used in sorting, generating range policies).

        // Some timers
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
    Here are just some debug functions, like a nice output.
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
