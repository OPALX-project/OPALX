#ifndef BIN_HISTO_H
#define BIN_HISTO_H

#include "Ippl.h"
#include "BinningTools.h" // For postSum computation 

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
    enum class HistoTypeIdentifier {
        Histogram, ///< Main histogram view.
        BinWidth,  ///< Bin width view.
        PostSum    ///< Post-sum view.
    };


    template <typename size_type, typename bin_index_type, typename value_type, 
              bool UseDualView = false, class... Properties>
    class Histogram {
    public:

        using view_type  = std::conditional_t<UseDualView, 
                                             Kokkos::DualView<size_type*, Properties...>, 
                                             Kokkos::View<size_type*, Properties...>>;
        using dview_type = std::conditional_t<UseDualView, 
                                              typename view_type::t_dev, 
                                              view_type>;
        using hview_type = std::conditional_t<UseDualView, 
                                              typename view_type::t_host, 
                                              view_type>;

        using width_view_type = std::conditional_t<UseDualView,
                                                   Kokkos::DualView<value_type*, Properties...>,
                                                   Kokkos::View<value_type*, Properties...>>;
                                                   
        using hwidth_view_type = std::conditional_t<UseDualView, 
                                                   typename width_view_type::t_host, 
                                                   width_view_type>;

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
        Histogram(std::string debug_name, bin_index_type numBins, value_type totalBinWidth)
            : debug_name_m(debug_name)
            , numBins_m(numBins)
            , totalBinWidth_m(totalBinWidth) {
            instantiateHistograms();
        }

        /**
         * @brief Default destructor for the Histogram class.
         */
        ~Histogram() = default; 

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
         * should be manageable. For better efficiency, one can avoid the Kokkos::View "copying-action"
         * by using dualView.h_view(binIndex).
         *
         * @tparam UseDualView A boolean template parameter indicating whether DualView is used.
         * @param binIndex The index of the bin for which the number of particles is to be retrieved.
         * @return The number of particles in the specified bin.
         */
        size_type getNPartInBin(bin_index_type binIndex) {
            if constexpr (UseDualView) {
                histogram_m.hview(binIndex);
            } else {
                std::cerr << "Warning: Accessing BinHisto.getNPartInBin without DualView might be inefficient!" << std::endl;
                Kokkos::View<size_type, Kokkos::HostSpace> host_scalar("host_scalar");
                Kokkos::deep_copy(host_scalar, Kokkos::subview(histogram_m, binIndex));
                return host_scalar();
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
            static IpplTimings::TimerRef histoInitTimer = IpplTimings::getTimer("syncInitHistoTools");

            // Assumes you have initialized histogram_m from the outside!
            IpplTimings::startTimer(histoInitTimer);
            sync();
            initConstBinWidths(totalBinWidth_m);
            initPostSum();
            IpplTimings::stopTimer(histoInitTimer);
        }

        /**
         * @brief Initializes the bin widths with a constant value.
         *
         * @param constBinWidth The constant value to set for all bin widths.
         * 
         * @note Should not be called again after merging bins, since the bin widths will all be different.
         */
        void initConstBinWidths(const value_type constBinWidth) {
            Kokkos::deep_copy(getDeviceView(binWidths_m), constBinWidth);
            if constexpr (UseDualView) {
                binWidths_m.modify_device();
                binWidths_m.sync_host();
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
            computeFixSum<size_type>(view_device(), getDeviceView(postSum_m));
            if constexpr (UseDualView) {
                postSum_m.modify_device();
                postSum_m.sync_host();
            }
        }

        /**
         * @brief Returns a Kokkos::RangePolicy for iterating over the elements in a specified bin.
         *
         * This function generates a range policy for iterating over the elements within a given bin index.
         * If no DualView is used, it needs to copy some values to host, which might cause overhead.
         *
         * @tparam bin_index_type The type of the bin index.
         * @param binIndex The index of the bin for which the iteration policy is to be generated.
         * @return Kokkos::RangePolicy<> The range policy for iterating over the elements in the specified bin.
         */
        Kokkos::RangePolicy<> getBinIterationPolicy(const bin_index_type& binIndex) {
            if constexpr (UseDualView) {
                // localPostSumHost = postSum_m.view_host();
                hview_type localPostSumHost = getHostView(postSum_m);
                return Kokkos::RangePolicy<>(localPostSumHost(binIndex), localPostSumHost(binIndex + 1));
            } else {
                std::cerr << "Warning: Accessing BinHisto.getBinIterationPolicy without DualView might be inefficient!" << std::endl;
                Kokkos::View<bin_index_type[2], Kokkos::HostSpace> host_ranges("host_scalar");
                Kokkos::deep_copy(host_ranges, Kokkos::subview(postSum_m, std::make_pair(binIndex, binIndex + 1)));
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
         * @brief This function provides access to different types of histogram data stored on the host. Depending on the provided histogram type identifier, it returns the corresponding host view.
         * 
         * @param histo_type The type of histogram data to retrieve. It can be one of the following:
         *                   - HistoTypeIdentifier::Histogram: Returns the host view of the main histogram data.
         *                   - HistoTypeIdentifier::BinWidth: Returns the host view of the bin widths.
         *                   - HistoTypeIdentifier::PostSum: Returns the host view of the post-sum data.
         * 
         * @return hview_type The host view corresponding to the specified histogram type.
         * 
         * @note If an unknown histogram type identifier is provided, an error message is printed,
         *       the program is aborted.
         */
        template <typename return_type = hview_type>
        return_type view_host(HistoTypeIdentifier histo_type = HistoTypeIdentifier::Histogram) { 
            switch (histo_type) {
                case HistoTypeIdentifier::Histogram:
                    return getHostView<return_type>(histogram_m);
                case HistoTypeIdentifier::BinWidth:
                    return getHostView<return_type>(binWidths_m);
                case HistoTypeIdentifier::PostSum:
                    return getHostView<return_type>(postSum_m);
                default:
                    std::cerr << "Error: Unknown histogram type identifier!" << std::endl;
                    ippl::Comm->abort();
                    return getHostView<return_type>(histogram_m); // just so it compiles...
            }
        }

        /**
         * @brief Returns a device view of the histogram data based on the specified histogram type.
         *
         * @param histo_type The type of histogram data to view. Default is `HistoTypeIdentifier::Histogram`.
         * @return dview_type A device view of the requested histogram data.
         *
         * @see view_host() for similar functionality on the host side and more explanation.
         */
        template <typename return_type = dview_type>
        return_type view_device(HistoTypeIdentifier histo_type = HistoTypeIdentifier::Histogram) { 
            switch (histo_type) {
                case HistoTypeIdentifier::Histogram:
                    return getDeviceView<return_type>(histogram_m);
                case HistoTypeIdentifier::BinWidth:
                    return getDeviceView<return_type>(binWidths_m);
                case HistoTypeIdentifier::PostSum:
                    return getDeviceView<return_type>(postSum_m);
                default:
                    std::cerr << "Error: Unknown histogram type identifier!" << std::endl;
                    ippl::Comm->abort();
                    return getDeviceView<return_type>(histogram_m); // just so it compiles...
            }
        }

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
        static constexpr return_type getDeviceView(HistogramType& histo) {
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
        static constexpr return_type getHostView(HistogramType& histo) {
            if constexpr (UseDualView) {
                return histo.view_host();
            } else {
                return histo;
            }
        }

        /*
        The following contain functions that are used to make the histogram adaptive.
        */

        void mergeBins(const value_type maxBinRatio);

    private:
        std::string debug_name_m;   /// \brief Debug name for identifying the histogram instance.
        bin_index_type numBins_m;   /// \brief Number of bins in the histogram.
        value_type totalBinWidth_m; /// \brief Total width of all bins combined.

        view_type       histogram_m;      /// \brief View storing the particle counts in each bin.
        width_view_type binWidths_m;      /// \brief View storing the widths of the bins.
        view_type       postSum_m;        /// \brief View storing the cumulative sum of bin counts (used in sorting, generating range policies).

        /**
         * @brief Copies the fields from another Histogram object.
         *
         * This function copies the internal fields from the provided Histogram object
         * to the current object. The fields are copied using Kokkos' shallow copy.
         *
         * @param other The Histogram object from which to copy the fields.
         */
        void copyFields(const Histogram& other);
    };
    
}

#include "BinHisto.hpp"

#endif // BIN_HISTO_H
