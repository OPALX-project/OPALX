#ifndef BINHISTO_TPP
#define BINHISTO_TPP

#include "BinHisto.h"

// #include <random>

namespace ParticleBinning {

    template <typename size_type, typename bin_index_type, typename value_type, bool UseDualView, class... Properties>
    void Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::copyFields(const Histogram& other) {
        debug_name_m    = other.debug_name_m;
        numBins_m       = other.numBins_m;
        totalBinWidth_m = other.totalBinWidth_m;

        binningAlpha_m = other.binningAlpha_m;
        binningBeta_m = other.binningBeta_m;
        desiredWidth_m = other.desiredWidth_m;

        histogram_m = other.histogram_m;
        binWidths_m = other.binWidths_m;
        postSum_m   = other.postSum_m;

        initTimers();
    }


    template <typename size_type, typename bin_index_type, typename value_type, bool UseDualView, class... Properties>
    value_type 
    Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::adaptiveBinningCostFunction(
        const size_type& sumCount,
        const value_type& sumWidth,
        const size_type& totalNumParticles
    ) {
        # ifdef DEBUG
        if (sumCount == 0) {
            Inform err("mergeBins");
            err << "Error in adaptiveBinningCostFunction: " 
                << "sumCount = " << sumCount
                << ", sumWidth = " << sumWidth
                << ", totalNumParticles = " << totalNumParticles << endl;
            ippl::Comm->abort();
        }
        # endif

        value_type totalSum     = static_cast<value_type>(totalNumParticles);
        value_type sumCountNorm = sumCount / totalSum; 

        value_type penalty        = sumWidth - desiredWidth_m;
        value_type wideBinPenalty = binningAlpha_m;
        value_type binSizeBias    = binningBeta_m; 

        value_type sparse_penalty = (sumCountNorm < desiredWidth_m) 
                                     ? desiredWidth_m / sumCountNorm 
                                     : 0.0;

        return sumCountNorm*log(sumCountNorm)*sumWidth // minimize shannon entropy as a basis
                + wideBinPenalty * pow(sumWidth, 2)    // >0 wants smallest possible bin
                                                       // <0 wants largest possible bin
                                                       // Use ^3 to make it reasonably sensitive
                + binSizeBias * pow(penalty, 2)        // bias towards desiredWidth
                + sparse_penalty;                      // penalty for too few particles (specifically "distribution tails")
    }


    template <typename size_type, typename bin_index_type, typename value_type, bool UseDualView, class... Properties>
    Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::hindex_transform_type
    Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::mergeBins() {
        Inform m("Histogram");

        /*
        The following if makes sure that the mergeBins function is only called if the histogram is
        actually available on host!.
        */
        if constexpr (!std::is_same<typename hview_type::memory_space, Kokkos::HostSpace>::value) {
            m << "This does not work if the histogram is not saved in a DualView, since it needs host access to the data." << endl;
            ippl::Comm->abort();
            return hindex_transform_type("error", 0);
        }

        // Get host views
        hview_type oldHistHost       = getHostView<hview_type>(histogram_m);
        hwidth_view_type oldBinWHost = getHostView<hwidth_view_type>(binWidths_m);

        const bin_index_type n = numBins_m;
        if (n < 2) {
            // Should not happen, since this function is to be called after generating a very fine histogram, e.g. 128 bins
            m << "Not merging, since n_bins = " << n << " is too small!" << endl;
            hindex_transform_type oldToNewBinsView("oldToNewBinsView", n);
            Kokkos::deep_copy(oldToNewBinsView, 0);
            return oldToNewBinsView;
        }


        IpplTimings::startTimer(bMergeBinsT);
        // ----------------------------------------------------------------
        // 1) Build prefix sums on the host
        //    prefixCount[k] = sum of counts in bins [0..k-1]
        //    prefixWidth[k] = sum of widths in bins [0..k-1]
        // ----------------------------------------------------------------
        hview_type       prefixCount("prefixCount", n+1);
        hwidth_view_type prefixWidth("prefixWidth", n+1);
        
        prefixCount(0)  = 0;
        prefixWidth(0)  = 0;
        
        for (bin_index_type i = 0; i < n; ++i) {
            prefixCount(i+1) = prefixCount(i) + oldHistHost(i);
            prefixWidth(i+1) = prefixWidth(i) + oldBinWHost(i);
        }
        const size_type totalNumParticles = prefixCount(n); // Last value in prefixCount is the total number of particles


        // ----------------------------------------------------------------
        // 2) Dynamic Programming arrays:
        //    dp(k)      = minimal total cost covering [0..k-1]
        //    prevIdx(k) = the index i that yields that minimal cost
        // ----------------------------------------------------------------
        Kokkos::View<value_type*, Kokkos::HostSpace> dp("dp", n+1);
        // Kokkos::View<value_type*, Kokkos::HostSpace> dpMoment("dpMoment", n+1); // Store cumulative moments --> allow first order moment error estimation
        Kokkos::View<int*,        Kokkos::HostSpace> prevIdx("prevIdx", n+1);

        // Initialize dp with something large
        value_type largeVal = std::numeric_limits<value_type>::max() / value_type(2);
        for (bin_index_type k = 0; k <= n; ++k) {
            dp(k)       = largeVal;
            // dpMoment(k) = 0; 
            prevIdx(k)  = -1;
        }
        dp(0) = value_type(0);  // 0 cost to cover an empty set 
        //m << "DP arrays initialized." << endl;


        // ----------------------------------------------------------------
        // 3) Fill dp with an O(n^2) algorithm to find the minimal total cost
        // ----------------------------------------------------------------
        for (bin_index_type k = 1; k <= n; ++k) {
            // Try all possible start indices i for the last merged bin
            for (bin_index_type i = 0; i < k; ++i) {
                size_type  sumCount      = prefixCount(k) - prefixCount(i);
                value_type sumWidth      = prefixWidth(k) - prefixWidth(i);
                value_type segCost       = largeVal;
                if (sumCount > 0) {
                    // Division by 0 if called with sumCount=0. Merge anyways --> set cost to a "very large" value instead
                    segCost = adaptiveBinningCostFunction(sumCount, sumWidth, totalNumParticles);
                }
                value_type candidate = dp(i) + segCost;
                if (candidate < dp(k)) {
                    dp(k)      = candidate;
                    prevIdx(k) = i;
                }
            }
        }
        //m << "DP arrays filled." << endl;

        // dp(n) is the minimal total cost for covering [0..n-1].
        value_type totalCost = dp(n);
        if (totalCost >= largeVal) {
            // Means everything was effectively "impossible" => fallback (should not happen!)
            std::cerr << "Warning: no feasible merges found. Setting cost=0, no merges." << std::endl;
            totalCost = value_type(0);
        }


        // ----------------------------------------------------------------
        // 4) Reconstruct boundaries from prevIdx
        //    We start from k=n and step backwards until k=0
        // ----------------------------------------------------------------
        std::vector<int> boundaries;
        boundaries.reserve(20); // should be sufficient for most use cases (usually aim for <10)
        int cur = n;
        // We'll just push them in reverse
        while (cur > 0) {
            int start = prevIdx(cur);
            if (start < 0) {
                std::cerr << "Error: prevIdx(" << cur << ") < 0. "
                            << "Merging not successful, aborted loop." << std::endl;
                // fallback, break out
                break;
            }
            boundaries.push_back(start);
            cur = start;
        }
        // boundaries are reversed (e.g. [startK, i2, i1, 0])
        std::reverse(boundaries.begin(), boundaries.end());
        // final boundary is n
        boundaries.push_back(n);

        // Now the number of merged bins is boundaries.size() - 1
        size_type mergedBinsCount = static_cast<size_type>(boundaries.size()) - 1;
        //m << "Merged bins (based on minimal cost partition): " << mergedBinsCount << ". Minimal total cost = " << totalCost << endl;


        // ----------------------------------------------------------------
        // 5) Build new arrays for the merged bins
        // ----------------------------------------------------------------
        Kokkos::View<size_type*,  Kokkos::HostSpace> newCounts("newCounts", mergedBinsCount);
        Kokkos::View<value_type*, Kokkos::HostSpace> newWidths("newWidths", mergedBinsCount);

        for (size_type j = 0; j < mergedBinsCount; ++j) {
            bin_index_type start = boundaries[j];
            bin_index_type end   = boundaries[j+1] - 1;  // inclusive
            size_type  sumCount  = prefixCount(end+1) - prefixCount(start);
            value_type sumWidth  = prefixWidth(end+1) - prefixWidth(start);
            newCounts(j) = sumCount;
            newWidths(j) = sumWidth;
        }
        //m << "New bins computed." << endl;

        // Also generate a lookup table that maps the old bin index to the new bin index (for "rebin")
        hindex_transform_type oldToNewBinsView("oldToNewBinsView", n);
        for (size_type j = 0; j < mergedBinsCount; ++j) {
            bin_index_type startIdx = boundaries[j];
            bin_index_type endIdx   = boundaries[j+1]; // exclusive
            for (bin_index_type i = startIdx; i < endIdx; ++i) {
                oldToNewBinsView(i) = j;
            }
        }
        //m << "Lookup table generated." << endl;


        // ----------------------------------------------------------------
        // 6) Overwrite the old histogram arrays with the new merged ones
        // ----------------------------------------------------------------
        numBins_m = static_cast<bin_index_type>(mergedBinsCount);
        instantiateHistograms();
        //m << "New histograms instantiated." << endl;

        // Copy the data into the new Kokkos Views (on host)
        hview_type newHistHost        = getHostView<hview_type>(histogram_m);
        hwidth_view_type newWidthHost = getHostView<hwidth_view_type>(binWidths_m);
        Kokkos::deep_copy(newHistHost, newCounts);
        Kokkos::deep_copy(newWidthHost, newWidths);
        //m << "New histograms filled." << endl;


        // ----------------------------------------------------------------
        // 7) If using DualView, mark host as modified & sync
        // ----------------------------------------------------------------
        if constexpr (UseDualView) {
            modify_host(); 
            sync();

            binWidths_m.modify_host();

            IpplTimings::startTimer(bDeviceSyncronizationT);
            binWidths_m.sync_device();
            IpplTimings::stopTimer(bDeviceSyncronizationT);
            //m << "Host views modified/synced." << endl;
        }


        // ----------------------------------------------------------------
        // 8) Recompute postSum for the new merged histogram
        // ----------------------------------------------------------------
        initPostSum();
        IpplTimings::stopTimer(bMergeBinsT);

        m << "Re-binned from " << n << " bins down to "
          << numBins_m << " bins. Total deviation cost = "
          << totalCost << endl;

        // Return the old->new index transform
        return oldToNewBinsView;
    }


} // namespace ParticleBinning

#endif // BINHISTO_TPP