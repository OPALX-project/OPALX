#ifndef BINHISTO_HPP
#define BINHISTO_HPP

namespace ParticleBinning {

    // We'll reduce to find the minimal candidate and the index i that achieves it.
    /*template <typename value_type, typename bin_index_type>
    struct MinPair {
        value_type val;
        bin_index_type idx;
    };

    template <typename value_type, typename bin_index_type>
    static KOKKOS_INLINE_FUNCTION
    MinPair<value_type, bin_index_type> combineMinPair(const MinPair<value_type, bin_index_type>& a, 
                                                       const MinPair<value_type, bin_index_type>& b) {
        return (a.val < b.val) ? a : b;
    }*/


    template <typename size_type, typename bin_index_type, typename value_type, bool UseDualView, class... Properties>
    void Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::copyFields(const Histogram& other) {
        debug_name_m    = other.debug_name_m;
        numBins_m       = other.numBins_m;
        totalBinWidth_m = other.totalBinWidth_m;

        histogram_m = other.histogram_m;
        binWidths_m = other.binWidths_m;
        postSum_m   = other.postSum_m;
    }


    // Implementation of the cost function
    /*template <typename size_type, typename bin_index_type, typename value_type, bool UseDualView, class... Properties>
    KOKKOS_INLINE_FUNCTION value_type
    Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::computeDeviationCost(
        const size_type& sumCount, const value_type& sumWidth, const value_type& maxBinRatio,
        const value_type& alpha, const value_type& mergedStd
    ) {
        value_type h = alpha * mergedStd * pow(sumCount, 3); // from "/ sumCount^1/3"
        return pow(maxBinRatio - sumWidth * h, 2);

        // Written such that there if no warp divergence
        // Compute the ratio
        //value_type ratio = sumWidth * sqrt(static_cast<value_type>(sumCount));
        
        // Select the cost based on the conditions
        // Explanation: when the bin is empty, the cost is "inf" to guarantee merging this bin
        //value_type costWhenSumCountZero = (sumWidth == value_type(0)) ? std::fabs(0.0 - maxBinRatio) : largeVal;
        //value_type costWhenSumCountNonZero = std::fabs(ratio - maxBinRatio);
        
        // Use a conditional select to avoid branching
        //return (sumCount == 0) ? costWhenSumCountZero : costWhenSumCountNonZero;
    }*/

    template <typename size_type, typename bin_index_type, typename value_type, bool UseDualView, class... Properties>
    template <typename BinningSelector_t>
    value_type 
    Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::partialMergedCDFIntegralCost(
        const bin_index_type& i, const bin_index_type& k,
        const size_type& sumCount,
        const value_type& sumWidth,
        const value_type& wideBinPenalty,
        const hash_type sortedIndexArr,
        const BinningSelector_t var_selector
        // const value_type& alpha, // TODO: maybe later
        //const value_type& mergedStd,
        //const value_type& segFineMoment
    ) {
        //std::cout << "i = " << i << ", k = " << k << ", sumCount = " << sumCount << ", sumWidth = " << sumWidth << std::endl;
        // 1. Calculate mean of segment to be merged
        value_type valueMean = 0.0;
        Kokkos::parallel_reduce(getBinIterationPolicy(i, k-i), KOKKOS_LAMBDA(const bin_index_type& j, value_type& valueMean) {
            const bin_index_type idx = sortedIndexArr(j);
            valueMean += var_selector(idx);
        }, valueMean);
        valueMean /= sumCount;
        //std::cout << "Hey there!" << std::endl;

        // 2. Calculate variance of segment to be merged
        value_type squaredDiff = 0.0;
        Kokkos::parallel_reduce(getBinIterationPolicy(i, k-i), KOKKOS_LAMBDA(const bin_index_type& j, value_type& squaredDiff) {
            const bin_index_type idx = sortedIndexArr(j);
            squaredDiff += pow(var_selector(idx) - valueMean, 2);
        }, squaredDiff);

        //std::cout << "Hey there!" << std::endl;

        value_type varEstimate = squaredDiff / sumCount;

        // Basically assuming Gaussian distribution per bin and calculating the log-likelihood. 
        // Aims to give a penalty for bins too wide and bins with not enough particles.
        return sumCount / 2 * log(2*3.1415926) 
                + sumCount * log(sqrt(varEstimate))
                + 1 / (2 * varEstimate) * squaredDiff
                + wideBinPenalty * sumWidth;
        /*
        // Compute the representative value as the weighted average. 
        const value_type representative      = segFineMoment / static_cast<value_type>(sumCount);
        const value_type mergedSegmentMoment = representative * sumCount;

        // Without any offset, this is just 0... (TODO: later)
        const value_type costMoment = pow(segFineMoment - mergedSegmentMoment, 2);

        // higher internal spread (mergedStd) and a larger (sumWidth)
        // lead to a larger integrated difference between the fine and merged CDFs 
        const value_type costCDF = mergedStd * mergedStd; // * sumWidth;

        // Add a penalty to discourage segments that are too narrow --> TODO: might want this to balance modeling error!
        const value_type dmin = 2.34e-4; 
        const value_type narrowPenalty = (sumWidth < dmin) ? pow(dmin - sumWidth, 2) : 0.0;

        // Note: At the moment, costMoment is zero :( -- TODO: what can I do here???
        // return alpha * costCDF + costMoment + narrowPenalty;
        return alpha * costCDF + mergedStd / pow(sumCount, 1.0/3.0) + narrowPenalty;*/
    }

    template <typename size_type, typename bin_index_type, typename value_type, bool UseDualView, class... Properties>
    value_type
    Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::mergedBinStd(
        const bin_index_type& i, const bin_index_type& k,
        const size_type& sumCount, const value_type& varPerBin, 
        const hwidth_view_type& prefixWidth, 
        const hview_type& fineCounts, const hwidth_view_type& fineWidths
    ) {
        // Calculate the weighted mean
        value_type weightedMean = 0;
        for (bin_index_type j = i; j < k; ++j) {
            value_type binMidpoint = prefixWidth(j) + 0.5 * fineWidths(j);
            weightedMean += fineCounts(j) * binMidpoint;
        }
        weightedMean /= sumCount;

        // Calculate the variance
        value_type variance = 0;
        for (bin_index_type j = i; j < k; ++j) {
            value_type binMidpoint = prefixWidth(j) + 0.5 * fineWidths(j);
            variance += fineCounts(j) * ( pow(binMidpoint - weightedMean, 2)+ varPerBin );
        }
        variance /= sumCount;

        return sqrt(variance);
    }


    template <typename size_type, typename bin_index_type, typename value_type, bool UseDualView, class... Properties>
    template <typename BinningSelector_t>
    Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::hindex_transform_type
    Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::mergeBins(
        const hash_type sortedIndexArr,
        const BinningSelector_t var_selector
    ) {
        static IpplTimings::TimerRef mergeBinsTimer = IpplTimings::getTimer("mergeBins");

        // Maybe set this later as a parameter
        value_type alpha = 1.0; // scotts normal reference rule, see https://en.wikipedia.org/wiki/Histogram#Scott's_normal_reference_rule
        
        // TODO 
        // Should merge neighbouring bins such that the width/N_part ratio is roughly maxBinRatio.
        // TODO: Find algorithm for that
        // std::cout << "Warning: mergeBins not implemented yet!" << std::endl;
        Inform m("Histogram");
        //m << "Merging bins with cost-based approach (minimize deviation from maxBinRatio = "
        //  << maxBinRatio << ")" << endl;

        if constexpr (!UseDualView) {
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

        IpplTimings::startTimer(mergeBinsTimer);
        // ----------------------------------------------------------------
        // 1) Build prefix sums on the host
        //    prefixCount[k] = sum of counts in bins [0..k-1]
        //    prefixWidth[k] = sum of widths in bins [0..k-1]
        // ----------------------------------------------------------------
        hview_type       prefixCount("prefixCount", n+1);
        hwidth_view_type prefixWidth("prefixWidth", n+1);
        Kokkos::View<value_type*, Kokkos::HostSpace> prefixMoment("prefixMoment", n+1); // Needed for first order moment error estimation
        prefixCount(0)  = 0;
        prefixWidth(0)  = 0;
        prefixMoment(0) = 0;
        for (bin_index_type i = 0; i < n; ++i) {
            prefixCount(i+1) = prefixCount(i) + oldHistHost(i);
            prefixWidth(i+1) = prefixWidth(i) + oldBinWHost(i);

            value_type binCenter = prefixWidth(i) + 0.5 * oldBinWHost(i); // Technically not necessary, but more general for non-uniform bins...
            prefixMoment(i+1) = prefixMoment(i) + oldHistHost(i) * binCenter; // Something like the cumulative distribution function for the "actual" histogram (fine bins...)
                                                                              // Basically the "integral" \int_{0}^{x} x f(x) dx
                                                                              // TODO: might want to use different integration rule?
        }
        //computeFixSum<hview_type>(oldHistHost, prefixCount);
        //computeFixSum<hwidth_view_type>(oldBinWHost, prefixWidth);

        m << "Prefix sums computed." << endl;


        // ----------------------------------------------------------------
        // 2) Dynamic Programming arrays:
        //    dp(k)      = minimal total cost covering [0..k-1]
        //    prevIdx(k) = the index i that yields that minimal cost
        // ----------------------------------------------------------------
        // We'll store dp as a floating-point "value_type" array
        Kokkos::View<value_type*, Kokkos::HostSpace> dp("dp", n+1);
        Kokkos::View<value_type*, Kokkos::HostSpace> dpMoment("dpMoment", n+1); // Store cumulative moments --> allow first order moment error estimation
        Kokkos::View<int*,        Kokkos::HostSpace> prevIdx("prevIdx", n+1);

        // Initialize dp with something large
        value_type largeVal = std::numeric_limits<value_type>::max() / value_type(2);
        for (bin_index_type k = 0; k <= n; ++k) {
            dp(k)       = largeVal;
            dpMoment(k) = 0; // Added this!
            prevIdx(k)  = -1;
        }
        dp(0) = value_type(0);  // 0 cost to cover an empty set (dpMoment(0) = 0, for no bins...)
        m << "DP arrays initialized." << endl;


        // ----------------------------------------------------------------
        // 3) Fill dp with an O(n^2) algorithm to find the minimal total cost
        // ----------------------------------------------------------------
        // value_type varPerBin = pow(oldBinWHost(0), 2) / 12; // assume equal width, assume uniform distribution per fine bin
        for (bin_index_type k = 1; k <= n; ++k) {
            // Try all possible start indices i for the last merged bin
            for (bin_index_type i = 0; i < k; ++i) {
                size_type  sumCount      = prefixCount(k) - prefixCount(i);
                value_type sumWidth      = prefixWidth(k) - prefixWidth(i);
                value_type segCost       = largeVal;
                if (sumCount > 0) {
                    //value_type segFineMoment = prefixMoment(k) - prefixMoment(i); // "exact" integral value for first order moment (from fine histo)
                    //value_type mergedStd     = mergedBinStd(i, k, sumCount, varPerBin, prefixWidth, oldHistHost, oldBinWHost);
                    //segCost              = computeDeviationCost(sumCount, sumWidth, maxBinRatio, alpha, mergedStd);
                    //segCost = partialMergedCDFIntegralCost(sumCount, sumWidth, alpha, mergedStd, segFineMoment);
                    segCost = partialMergedCDFIntegralCost(i, k, sumCount, sumWidth, alpha, sortedIndexArr, var_selector);

                    if (k % 10 == 0 && i % 10 == 0) {
                        m << "k = " << k << ", i = " << i << ", sumCount = " << sumCount << ", sumWidth = " << sumWidth
                          << ", segCost = " << segCost << endl; // ", mergedStd = " << mergedStd << endl;
                    }
                }
                //value_type segCost   = computeDeviationCost(sumCount, sumWidth, maxBinRatio, largeVal, alpha);
                value_type candidate = dp(i) + segCost;
                if (candidate < dp(k)) {
                    dp(k)      = candidate;
                    prevIdx(k) = i;
                }
            }
        }

        m << "DP arrays filled." << endl;

        // dp(n) is the minimal total cost for covering [0..n-1].
        value_type totalCost = dp(n);
        if (totalCost >= largeVal) {
            // Means everything was effectively "impossible" => fallback
            std::cerr << "Warning: no feasible merges found. Setting cost=0, no merges." << std::endl;
            totalCost = value_type(0);
        }

        //for (bin_index_type k = 0; k <= n; ++k) {
        //    m << "dp(" << k << ") = " << dp(k) << ", prevIdx(" << k << ") = " << prevIdx(k) << endl;
        //}


        // ----------------------------------------------------------------
        // 4) Reconstruct boundaries from prevIdx
        //    We start from k=n and step backwards until k=0
        // ----------------------------------------------------------------
        std::vector<int> boundaries;
        boundaries.reserve(20); // should be sufficient for most use cases
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
        // boundaries is reversed (e.g. [startK, i2, i1, 0])
        std::reverse(boundaries.begin(), boundaries.end());
        // final boundary is n
        boundaries.push_back(n);

        // Now the number of merged bins is boundaries.size() - 1
        size_type mergedBinsCount = static_cast<size_type>(boundaries.size()) - 1;
        m << "Merged bins (based on minimal cost partition): " << mergedBinsCount << ". Minimal total cost = " << totalCost << endl;



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
        m << "New bins computed." << endl;



        // Also generate a lookup table that maps the old bin index
        // to the new bin index
        hindex_transform_type oldToNewBinsView("oldToNewBinsView", n);
        for (size_type j = 0; j < mergedBinsCount; ++j) {
            bin_index_type startIdx = boundaries[j];
            bin_index_type endIdx   = boundaries[j+1]; // exclusive
            for (bin_index_type i = startIdx; i < endIdx; ++i) {
                oldToNewBinsView(i) = j;
            }
        }
        m << "Lookup table generated." << endl;


        // ----------------------------------------------------------------
        // 6) Overwrite the old histogram arrays with the new merged ones
        // ----------------------------------------------------------------
        numBins_m = static_cast<bin_index_type>(mergedBinsCount);

        instantiateHistograms();
        m << "New histograms instantiated." << endl;

        // Copy the data into the new Kokkos Views (on host)
        hview_type newHistHost       = getHostView<hview_type>(histogram_m);
        hwidth_view_type newWidthHost= getHostView<hwidth_view_type>(binWidths_m);
        Kokkos::deep_copy(newHistHost, newCounts);
        Kokkos::deep_copy(newWidthHost, newWidths);
        m << "New histograms filled." << endl;


        // ----------------------------------------------------------------
        // 7) If using DualView, mark host as modified & sync
        // ----------------------------------------------------------------
        if constexpr (UseDualView) {
            modify_host(); 
            sync();

            binWidths_m.modify_host();
            binWidths_m.sync_device();
            m << "Host views modified/synced." << endl;
        }


        // ----------------------------------------------------------------
        // 8) Recompute postSum for the new histogram
        // ----------------------------------------------------------------
        initPostSum();
        IpplTimings::stopTimer(mergeBinsTimer);

        m << "Re-binned from " << n << " bins down to "
          << numBins_m << " bins. Total deviation cost = "
          << totalCost << endl;

        // Return the old->new index transform
        return oldToNewBinsView;
    }


} // namespace ParticleBinning

#endif // BINHISTO_HPP