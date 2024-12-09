#ifndef BIN_HISTO_H
#define BIN_HISTO_H

#include "Ippl.h"
#include "BinningTools.h" // For postSum computation 

namespace ParticleBinning {

    enum class HistoTypeIdentifier {
        Histogram,
        BinWidth,
        PostSum
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

        Histogram() = default;

        Histogram(std::string debug_name, bin_index_type numBins, value_type totalBinWidth)
            : debug_name_m(debug_name)
            , numBins_m(numBins)
            , totalBinWidth_m(totalBinWidth) {
            instantiateHistograms();
        }
        
        
        /*
        Some functions to access only single elements of the histogram.
        */

        size_type getNPartInBin(bin_index_type binIndex) {
            /**
             * Assume DualView was properly synchronized.
             * Might create some overhead from .view_host() call if called often.
             * However, it is only called on host (max nBins times per iteration), so should be fine. You can make it
             * more efficient by avoiding the Kokkos:View "copying-action" with e.g. dualView.h_view(binIndex).
             */
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
        
        // Called after initializing with constant bins width from inside (e.g.) AdaptBins
        void instantiateHistograms() {
            histogram_m = view_type("histogram", numBins_m);
            binWidths_m = view_type("binWidths", numBins_m - 1);
            postSum_m   = view_type("postSum", numBins_m + 1);
        }

        // Possibly updates the other histograms --> meant for when the histogram is built during the first call!
        void init() { // const value_type constBinWidth
            static IpplTimings::TimerRef histoInitTimer = IpplTimings::getTimer("syncInitHistoTools");

            // Assumes you have initialized histogram_m from the outside!
            IpplTimings::startTimer(histoInitTimer);
            sync();
            initConstBinWidths(totalBinWidth_m);
            initPostSum();
            IpplTimings::stopTimer(histoInitTimer);
        }

        void initConstBinWidths(const value_type constBinWidth) {
            dview_type binWidthsView = getDeviceView(binWidths_m);
            //constexpr auto binWidthsView = UseDualView ? binWidths_m.view_device() : binWidths_m;

            Kokkos::deep_copy(binWidthsView, constBinWidth);
            if constexpr (UseDualView) {
                binWidths_m.modify_device();
                binWidths_m.sync_host();
            }
        }

        void initPostSum() {
            //auto postSumView = constexpr UseDualView ? postSum_m.view_device() : postSum_m;
            dview_type postSumView = getDeviceView(postSum_m);
            computeFixSum<size_type>(view_device(), postSumView);
            if constexpr (UseDualView) {
                postSum_m.modify_device();
                postSum_m.sync_host();
            }
        }

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

        void modify_device() {
            if constexpr (UseDualView) histogram_m.modify_device();
        }

        void modify_host() {
            if constexpr (UseDualView) histogram_m.modify_host();
        }

        hview_type view_host(HistoTypeIdentifier histo_type = HistoTypeIdentifier::Histogram) { 
            switch (histo_type) {
                case HistoTypeIdentifier::Histogram:
                    return getHostView(histogram_m);
                case HistoTypeIdentifier::BinWidth:
                    return getHostView(binWidths_m);
                case HistoTypeIdentifier::PostSum:
                    return getHostView(postSum_m);
                default:
                    std::cerr << "Error: Unknown histogram type identifier!" << std::endl;
                    ippl::Comm->abort();
                    return getHostView(histogram_m); // just so it compiles...
            }
        }

        dview_type view_device(HistoTypeIdentifier histo_type = HistoTypeIdentifier::Histogram) { 
            switch (histo_type) {
                case HistoTypeIdentifier::Histogram:
                    return getDeviceView(histogram_m);
                case HistoTypeIdentifier::BinWidth:
                    return getDeviceView(binWidths_m);
                case HistoTypeIdentifier::PostSum:
                    return getDeviceView(postSum_m);
                default:
                    std::cerr << "Error: Unknown histogram type identifier!" << std::endl;
                    ippl::Comm->abort();
                    return getDeviceView(histogram_m); // just so it compiles...
            }
        }

        template <typename HistogramType>
        static constexpr dview_type getDeviceView(HistogramType& histo) {
            if constexpr (UseDualView) {
                return histo.view_device();
            } else {
                return histo;
            }
        }

        template <typename HistogramType>
        static constexpr hview_type getHostView(HistogramType& histo) {
            if constexpr (UseDualView) {
                return histo.view_host();
            } else {
                return histo;
            }
        }

    private:
        std::string debug_name_m;
        bin_index_type numBins_m;
        value_type totalBinWidth_m;

        view_type histogram_m;
        view_type binWidths_m;
        view_type postSum_m;
    };
    
}

#include "BinHisto.hpp"

#endif // BIN_HISTO_H
