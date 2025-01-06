#ifndef BINHISTO_HPP
#define BINHISTO_HPP

namespace ParticleBinning {


    template <typename size_type, typename bin_index_type, typename value_type, bool UseDualView, class... Properties>
    void Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::copyFields(const Histogram& other) {
        debug_name_m    = other.debug_name_m;
        numBins_m       = other.numBins_m;
        totalBinWidth_m = other.totalBinWidth_m;

        histogram_m = other.histogram_m;
        binWidths_m = other.binWidths_m;
        postSum_m   = other.postSum_m;
    }

    template <typename size_type, typename bin_index_type, typename value_type, bool UseDualView, class... Properties>
    Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::index_transform_type
    Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::mergeBins(value_type maxBinRatio) {
        // TODO 
        // Should merge neighbouring bins such that the width/N_part ratio is roughly maxBinRatio.
        // TODO: Find algorithm for that
        // std::cout << "Warning: mergeBins not implemented yet!" << std::endl;
        Inform m("Histogram");
        m << "Merging bins with maxBinRatio = " << maxBinRatio << endl;

        if constexpr (!UseDualView) {
            m << "This does not work if the histogram is not saved in a DualView, since it needs host access to the data." << endl;
            ippl::Comm->abort();
            return index_transform_type("error", 0);
        }

        // Get host views
        hview_type oldHistHost       = getHostView<hview_type>(histogram_m);
        hwidth_view_type oldBinWHost = getHostView<hwidth_view_type>(binWidths_m);

        const bin_index_type n = numBins_m;
        if (n < 2) {
            m << "Not merging, since n_bins = " << n << " is too small!" << endl;
            index_transform_type oldToNewBinsView("oldToNewBinsView", n);
            Kokkos::deep_copy(oldToNewBinsView, 0);
            return oldToNewBinsView;
        }


        // ----------------------------------------------------------------
        // 2) Build prefix sums on the host
        //
        //    prefixCount[k] = sum of counts in bins [0..k-1]
        //    prefixWidth[k] = sum of widths in bins [0..k-1]
        //
        //    i.e. prefixCount[0] = 0, prefixCount[1] = oldHistHost(0), etc.
        //
        //    Note: everything on host, since algorithm is not parallelized.
        // ----------------------------------------------------------------
        Kokkos::View<size_type*,  Kokkos::HostSpace> prefixCount("prefixCount", n+1);
        Kokkos::View<value_type*, Kokkos::HostSpace> prefixWidth("prefixWidth", n+1);
        prefixCount(0) = 0;
        prefixWidth(0) = 0;
        for (bin_index_type i = 0; i < n; ++i) {
            prefixCount(i+1) = prefixCount(i) + oldHistHost(i);
            prefixWidth(i+1) = prefixWidth(i) + oldBinWHost(i);
        }


        // ----------------------------------------------------------------
        // 3) Dynamic Programming arrays:
        //    dp(k)      = minimum number of merged bins covering [0..k-1]
        //    prevIdx(k) = the index i that yields that minimal partition
        // ----------------------------------------------------------------
        Kokkos::View<int*, Kokkos::HostSpace> dp("dp", n+1);
        Kokkos::View<int*, Kokkos::HostSpace> prevIdx("prevIdx", n+1);

        // Initialize dp with something large
        for (bin_index_type k = 0; k <= n; ++k) {
            dp(k)      = std::numeric_limits<int>::max();
            prevIdx(k) = -1;
        }
        dp(0) = 0;  // 0 bins needed to cover an empty set (base case)


        // ----------------------------------------------------------------
        // 4) Fill dp with an O(n^2) algorithm to find the optimal partition
        // ----------------------------------------------------------------
        for (bin_index_type k = 1; k <= n; ++k) {
            // Try all possible start indices i for the last merged bin
            for (bin_index_type i = 0; i < k; ++i) {
                size_type  sumCount  = prefixCount(k) - prefixCount(i);
                value_type sumWidth  = prefixWidth(k) - prefixWidth(i);

                // If sumCount==0 but sumWidth!=0, ratio = infinite => not valid
                // If sumCount==0 and sumWidth==0, ratio is 0/0 => interpret as 0 or skip
                if (sumCount == 0) {
                    if (sumWidth != value_type(0)) {
                        continue; // can't form a valid bin with ratio <= maxBinRatio
                    }
                    // else sumWidth=0 => ratio=0 => valid
                } else {
                    value_type ratio = sumWidth / static_cast<value_type>(sumCount);
                    if (ratio > maxBinRatio) {
                        continue;  // doesn't meet ratio constraint
                    }
                }

                // If we can form one valid merged bin from [i..k-1], 
                // check if that yields a better (fewer bins) partition
                int cost = dp(i) + 1; 
                if (cost < dp(k)) {
                    dp(k)      = cost;
                    prevIdx(k) = i;
                }
            }
        }

        int mergedBinsCount = dp(n);
        if (mergedBinsCount < 1) {
            // Shouldn't happen unless everything is zero?
            mergedBinsCount = 1;
        }


        // ----------------------------------------------------------------
        // 5) Reconstruct boundaries from prevIdx
        //
        //    We start from k=n and step backwards until k=0
        //    Each boundary = i => means the last bin covered [i..k-1]
        // ----------------------------------------------------------------
        std::vector<int> boundaries; 
        boundaries.reserve(mergedBinsCount+1);
        {
            int cur = n;
            while (cur > 0) {
                int start = prevIdx(cur); 
                if (start < 0) {
                    // fallback: something went wrong, break to avoid infinite loop
                    // this should not happen if dp was computed correctly
                    std::cerr << "Error: prevIdx(" << cur << ") < 0. Merging not successful, aborted loop." << std::endl;
                    break;
                }
                boundaries.push_back(start);
                cur = start;
            }
            // boundaries is in reverse order of merges, e.g. [startK, ..., 0]
            std::reverse(boundaries.begin(), boundaries.end());
            // now boundaries = [0, i1, i2, ..., iLast]
            // The final bin covers [iLast..n-1].
            boundaries.push_back(n);
        }

        // The number of bins we actually have is boundaries.size() - 1
        if (static_cast<int>(boundaries.size()) - 1 != mergedBinsCount) {
            // Just a consistency check
            mergedBinsCount = static_cast<int>(boundaries.size()) - 1;
        }


        // ----------------------------------------------------------------
        // 6) Build new arrays for the merged bins
        //
        //    For each merged bin j from 0..(mergedBinsCount-1):
        //       start = boundaries[j], end = boundaries[j+1]-1
        // ----------------------------------------------------------------
        std::vector<size_type>  newCounts(mergedBinsCount, 0);
        std::vector<value_type> newWidths(mergedBinsCount, value_type(0));

        for (int j = 0; j < mergedBinsCount; ++j) {
            int start = boundaries[j];
            int end   = boundaries[j+1] - 1;  // inclusive
            size_type  sumCount  = prefixCount(end+1) - prefixCount(start);
            value_type sumWidth  = prefixWidth(end+1) - prefixWidth(start);
            newCounts[j] = sumCount;
            newWidths[j] = sumWidth;
        }

        // Also generate a lookup table that maps the old bin index to the new bin index
        index_transform_type oldToNewBinsView("oldToNewBinsView", n);
        // For j in [0 .. mergedBinsCount-1], the range of old bins is [boundaries[j], boundaries[j+1]-1].
        // So fill oldToNewBinsView(i) = j for i in that range.
        for (bin_index_type j = 0; j < mergedBinsCount; ++j) {
            bin_index_type startIdx = boundaries[j];
            bin_index_type endIdx   = boundaries[j+1]; // exclusive
            for (bin_index_type i = startIdx; i < endIdx; ++i) {
                oldToNewBinsView(i) = j;
            }
        }


        // ----------------------------------------------------------------
        // 7) Overwrite the old histogram arrays with the new merged ones
        // ----------------------------------------------------------------
        numBins_m = static_cast<bin_index_type>(mergedBinsCount);

        instantiateHistograms();
        //histogram_m = view_type("histogram", newNumBins);
        //binWidths_m = width_view_type("binWidths", newNumBins); 

        // Copy data into the new Kokkos Views (on host)
        {
            // These return the views initialized view lines above
            hview_type newHistHost        = getHostView<hview_type>(histogram_m);
            hwidth_view_type newWidthHost = getHostView<hwidth_view_type>(binWidths_m);
            for (bin_index_type i = 0; i < numBins_m; ++i) {
                newHistHost(i)   = newCounts[i];
                newWidthHost(i)  = newWidths[i];
            }
        }


        // ----------------------------------------------------------------
        // 8) If using DualView, mark host as modified so that a future sync 
        //    will push to device if needed (should always happen)
        // ----------------------------------------------------------------
        if constexpr (UseDualView) {
            modify_host(); // histogram
            sync();

            binWidths_m.modify_host();
            binWidths_m.sync_device();
        }

        // ----------------------------------------------------------------
        // 9) Recompute postSum for the new histogram
        // ----------------------------------------------------------------
        initPostSum();

        m << "Re-binned from " << n << " bins down to " 
          << numBins_m << " bins (optimal partition for, ratio <= " 
          << maxBinRatio << ")." << endl;

        return oldToNewBinsView;
    }


} // namespace ParticleBinning

#endif // BINHISTO_HPP