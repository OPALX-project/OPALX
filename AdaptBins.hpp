#ifndef ADAPT_BINS_HPP
#define ADAPT_BINS_HPP

// #include "AdaptBins.h"

namespace ParticleBinning {

    template <typename BunchType>
    void AdaptBins<BunchType>::initLimits() {
        Inform msg("AdaptBins");  // INFORM_ALL_NODES

        Kokkos::MinMaxScalar<value_type> localMinMax;
        position_view_type localData = bunch_m->R.getView();
        
        Kokkos::parallel_reduce("localBinLimitReduction", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const int i, Kokkos::MinMaxScalar<value_type>& update) {
            value_type val = localData(i)[2]; // use z axis for binning!
            update.min_val = Kokkos::min(update.min_val, val);
            update.max_val = Kokkos::max(update.max_val, val);
        }, Kokkos::MinMax<value_type>(localMinMax));

        xMin_m = localMinMax.min_val;
        xMax_m = localMinMax.max_val;

        // Putting the same to-reduce variable as an argument ensures that every node gets the correct min/max and not just the root node!
        // Note: boradcast does not exist, use allreduce for reduce+broadcast together!
        ippl::Comm->allreduce(xMax_m, 1, std::greater<value_type>());
        ippl::Comm->allreduce(xMin_m, 1, std::less<value_type>());

        binWidth_m = (xMax_m - xMin_m) / currentBins_m;
        
        msg << "Initialized limits. Min: " << xMin_m << ", max: " << xMax_m << ", binWidth: " << binWidth_m << endl;
    }

    template <typename BunchType>
    void AdaptBins<BunchType>::initializeHistogram(bool setToZero) {
        // Reinitialize the histogram view with the new size (numBins)
        const bin_index_type numBins = getCurrentBinCount();
        localBinHisto_m = bin_histo_type("binHisto_m", numBins);
        
        // Optionally, initialize the histogram to zero
        if (setToZero) {
            bin_histo_type localBinHisto = localBinHisto_m; // avoid implicit "this" capture
            Kokkos::parallel_for("initHistogram", numBins, KOKKOS_LAMBDA(const bin_index_type i) {
                localBinHisto(i) = 0;
            });
        }
    }

    template <typename BunchType>
    KOKKOS_INLINE_FUNCTION typename AdaptBins<BunchType>::bin_index_type 
    AdaptBins<BunchType>::getBin(value_type x, value_type xMin, value_type xMax, value_type binWidthInv, bin_index_type numBins) {
        // Explanation: Don't access xMin, binWidth, ... through the members to avoid implicit
        // variable capture by Kokkos and potential copying overhead. Instead, pass them as an 
        // argument, s.t. Kokkos can capture them explicitly!
        // Make it static to avoid implicit capture of this inside Kokkos lambda! 

        // Ensure x is within bounds (clamp it between xMin and xMax --> this is only for bin assignment)
        // x = (x < xMin) ? xMin : ((x > xMax) ? xMax : x);
        x += (x < xMin) * (xMin - x) + (x > xMax) * (xMax - x); // puts x in the bin or nearest bin if out of bounds

        bin_index_type bin = (x - xMin) * binWidthInv; // multiply with inverse of binwidth
        return (bin >= numBins) ? (numBins - 1) : bin;  // Clamp to the maximum bin
    }

    template <typename BunchType>    
    void AdaptBins<BunchType>::assignBinsToParticles() {
        // Set the bin attribute for the given particle
        Inform msg("AdaptBins");

        position_view_type localData = bunch_m->R.getView();
        bin_view_type binIndex       = bunch_m->bin.getView();  

        // Declare the variables locally before the Kokkos::parallel_for (to avoid implicit this capture in Kokkos lambda)
        value_type xMin = xMin_m, xMax = xMax_m, binWidthInv = 1.0/binWidth_m;
        bin_index_type numBins = currentBins_m;
        // Alternatively explicit capture: [xMin = xMin_m, xMax = xMax_m, binWidth = binWidth_m, numBins = currentBins_m, localData = localData, binIndex = binIndex]

        static IpplTimings::TimerRef assignParticleBins = IpplTimings::getTimer("assignParticleBins");
        IpplTimings::startTimer(assignParticleBins);

        Kokkos::parallel_for("assignParticleBins", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const size_type i) {
                // Access the z-axis position of the i-th particle
                value_type v = localData(i)[2];  
                
                // Assign the bin index to the particle (directly on device)
                bin_index_type bin = getBin(v, xMin, xMax, binWidthInv, numBins);
                binIndex(i)        = bin;
        });
        IpplTimings::stopTimer(assignParticleBins);
        msg << "All bins assigned." << endl; 

        initLocalHisto();
        msg << "Local Histogram initialized." << endl;
    }

    template <typename BunchType>
    void AdaptBins<BunchType>::initLocalHisto() {
        Inform msg("AdaptBins");
        initializeHistogram(); // Init histogram (no need to set to 0, since contribute_into overwrites values...)
        msg << "Histogram initialized to 0" << endl;

        bin_view_type binIndex       = bunch_m->bin.getView();
        bin_histo_type localBinHisto = localBinHisto_m;
        //bin_index_type binCount      = getCurrentBinCount();
        
        msg << "Starting reducer." << endl;

        static IpplTimings::TimerRef initLocalHisto = IpplTimings::getTimer("initLocalHisto");
        IpplTimings::startTimer(initLocalHisto);

        //Kokkos::Experimental::ScatterView<size_type*> scatter_view("scatter_view", binCount);
        Kokkos::Experimental::ScatterView<size_type*> scatter_view(localBinHisto); // does not contribute if not done like this...
        Kokkos::parallel_for("initLocalHist", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const size_type i) {
                auto scatter = scatter_view.access();
                bin_index_type ndx = binIndex(i);
                scatter(ndx)++; // scatter(ndx)++;
            });

        msg << "Reduced, now contributing." << endl;
        scatter_view.contribute_into(localBinHisto); // does nothing if "scatter_view.subview().data() == localBinHisto_m.data()", so can as well remain just in case
        
        IpplTimings::stopTimer(initLocalHisto);

        msg << "Reducer ran without error." << endl;
    }

    template <typename BunchType>
    AdaptBins<BunchType>::bin_host_histo_type AdaptBins<BunchType>::getGlobalHistogram() {
        Inform msg("GetGlobalHistogram");

        // Get the current number of bins
        bin_index_type numBins = getCurrentBinCount(); // number of local bins = number of global bins!
        
        // Create a view to hold the global histogram on all ranks
        bin_host_histo_type globalBinHisto("globalBinHisto", numBins);

        // Need host mirror, otherwise the data is not available when the histogram is created using CUDA
        auto localBinHistoHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), localBinHisto_m);
        
        // Loop through each bin and reduce its value across all nodes
        for (bin_index_type i = 0; i < numBins; ++i) {
            size_type localValue = localBinHistoHost(i);
            size_type globalValue = 0;
            
            ippl::Comm->allreduce(localValue, globalValue, 1, std::plus<size_type>());  // Reduce to all nodes
            globalBinHisto(i) = globalValue; // Saves global histogram on ALL ranks...
        }

        // Return the global histogram (note: only meaningful on rank 0)
        return globalBinHisto;
    }

}

#endif // ADAPT_BINS_HPP


