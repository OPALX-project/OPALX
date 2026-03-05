#ifndef ADAPT_BINS_TPP
#define ADAPT_BINS_TPP

#include "AdaptBins.h"

namespace ParticleBinning {

    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::initLimits() {
        Inform msg("AdaptBins");  // INFORM_ALL_NODES

        // Update needed if bunch->create() is called between binnings!
        var_selector_m.updateDataArr(bunch_m); 
        BinningSelector var_selector = var_selector_m;  
        size_type nlocal             = bunch_m->getLocalNum();

        // Start limit reduction if number of particles is big enough.
        IpplTimings::startTimer(bInitLimitsT);
        if (nlocal <= 0) {
            msg << "Particles in the bunch = " << nlocal << ". Overwriting limits manually." << endl;
            xMin_m = xMax_m = (nlocal == 0) ? 0 : 0; 
        } else {
            Kokkos::MinMaxScalar<value_type> localMinMax;
            // Sadly this is necessary, since Kokkos seems to have a problem when nlocal == 1, where it does not update localMinMax
            if (nlocal == 1) {
                // Access and copy values explicitly
                Kokkos::View<value_type, Kokkos::HostSpace> host_scalar("host_scalar"); 
                Kokkos::View<value_type> tmp_dvalue("tmp_dvalue");
                Kokkos::parallel_for("tmp_dvalue", 1, KOKKOS_LAMBDA(const size_type) { tmp_dvalue() = var_selector(0); });
                Kokkos::deep_copy(host_scalar, tmp_dvalue);
                localMinMax.max_val = localMinMax.min_val = host_scalar();
            } else {
                Kokkos::parallel_reduce("localBinLimitReduction", nlocal, 
                                        KOKKOS_LAMBDA(const size_type i, Kokkos::MinMaxScalar<value_type>& update) {
                    value_type val = var_selector(i);
                    update.min_val = Kokkos::min(update.min_val, val);
                    update.max_val = Kokkos::max(update.max_val, val);
                }, Kokkos::MinMax<value_type>(localMinMax));
            }
            xMin_m = localMinMax.min_val;
            xMax_m = localMinMax.max_val;
        }
        IpplTimings::stopTimer(bInitLimitsT);

        IpplTimings::startTimer(bAllReduceLimitsT);
        ippl::Comm->allreduce(xMax_m, 1, std::greater<value_type>());
        ippl::Comm->allreduce(xMin_m, 1, std::less<value_type>());
        IpplTimings::stopTimer(bAllReduceLimitsT);

        // Update bin width variable with new limits
        binWidth_m = (xMax_m - xMin_m) / currentBins_m;
        msg << "Initialized limits. Min: " << xMin_m << ", max: " << xMax_m << ", binWidth: " << binWidth_m << endl;
    }


    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::instantiateHistogram(bool setToZero) {
        // Reinitialize the histogram view with the new size (numBins)
        const bin_index_type numBins = getCurrentBinCount();
        localBinHisto_m = d_histo_type("localBinHisto_m", numBins, xMax_m - xMin_m,
                                       binningAlpha_m, binningBeta_m, desiredWidth_m);
        
        // Optionally, initialize the histogram to zero
        if (setToZero) {
            dview_type device_histo = localBinHisto_m.template getDeviceView<dview_type>(localBinHisto_m.getHistogram());
            Kokkos::deep_copy(device_histo, 0);

            // Sync changed (only has an effect if DualView mode is used and the two memory spaces are different!)
            localBinHisto_m.modify_device();
            localBinHisto_m.sync();
        }
    }


    template <typename BunchType, typename BinningSelector>
    KOKKOS_INLINE_FUNCTION typename AdaptBins<BunchType, BinningSelector>::bin_index_type 
    AdaptBins<BunchType, BinningSelector>::getBin(value_type x, value_type xMin, value_type xMax, value_type binWidthInv, bin_index_type numBins) {
        // Explanation: Don't access xMin, binWidth, ... through the members to avoid implicit
        // variable capture by Kokkos and potential copying overhead. Instead, pass them as an 
        // argument, s.t. Kokkos can capture them explicitly!
        // Make it static to avoid implicit capture of this inside Kokkos lambda! 

        // Ensure x is within bounds (clamp it between xMin and xMax --> this is only for bin assignment).
        x += (x < xMin) * (xMin - x) + (x > xMax) * (xMax - x); 

        bin_index_type bin = (x - xMin) * binWidthInv; 
        return (bin >= numBins) ? (numBins - 1) : bin; // Clamp to the maximum bin if "out of bounds".
                                                       // Might happen if binning is called before initLimits()...
    }


    template <typename BunchType, typename BinningSelector>    
    void AdaptBins<BunchType, BinningSelector>::assignBinsToParticles() {
        Inform msg("AdaptBins");

        var_selector_m.updateDataArr(bunch_m);
        BinningSelector var_selector  = var_selector_m; 
        bin_view_type binIndex        = getBinView();

        IpplTimings::startTimer(bAssignUniformBinsT);
        if (bunch_m->getLocalNum() <= 1) {
            msg << "Too few bins, assigning all bins to index 0." << endl;
            Kokkos::deep_copy(binIndex, 0);
            return;
        }

        // Declare the variables locally before the Kokkos::parallel_for (to avoid implicit this capture in Kokkos lambda)
        value_type xMin = xMin_m, xMax = xMax_m, binWidthInv = 1.0/binWidth_m;
        bin_index_type numBins = currentBins_m;

        // Assign the bin index to the particle 
        Kokkos::parallel_for("assignParticleBinsConst", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const size_type& i) {
                value_type v = var_selector(i); 
                
                bin_index_type bin = getBin(v, xMin, xMax, binWidthInv, numBins);
                binIndex(i)        = bin;
        });
        IpplTimings::stopTimer(bAssignUniformBinsT);

        msg << "All bins assigned." << endl; 
    }


    template<typename BunchType, typename BinningSelector>
    template<typename ReducerType>
    void AdaptBins<BunchType, BinningSelector>::executeInitLocalHistoReduction(ReducerType& to_reduce) {
        bin_view_type binIndex        = getBinView(); 
        dview_type device_histo       = localBinHisto_m.template getDeviceView<dview_type>(localBinHisto_m.getHistogram());
        bin_index_type binCount       = getCurrentBinCount();

        // Reduce using an array-reducer object (size needs to be known at compile time for GPU kernels)
        Kokkos::parallel_reduce("initLocalHist", bunch_m->getLocalNum(), 
            KOKKOS_LAMBDA(const size_type& i, ReducerType& update) {
                bin_index_type ndx = binIndex(i);  // Determine the bin index for this particle
                update.the_array[ndx]++;           // Increment the corresponding bin count in the reduction array
            }, Kokkos::Sum<ReducerType>(to_reduce)
        );

        // Copy the reduced results to the final histogram (into the BinHisto class instance)
        Kokkos::parallel_for("finalize_histogram", binCount, 
            KOKKOS_LAMBDA(const bin_index_type& i) {
                device_histo(i) = to_reduce.the_array[i];
            }
        );

        // Note, sync is called inside initLocalHiso, not here
        localBinHisto_m.modify_device();
    }


    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::executeInitLocalHistoReductionTeamFor() {
        bin_view_type binIndex            = getBinView();
        dview_type device_histo           = localBinHisto_m.template getDeviceView<dview_type>(localBinHisto_m.getHistogram());
        const bin_index_type binCount     = getCurrentBinCount();
        const size_type localNumParticles = bunch_m->getLocalNum(); 

        using team_policy = Kokkos::TeamPolicy<>;
        using member_type = team_policy::member_type;

        // Define a scratch space view type
        using scratch_view_type = Kokkos::View<size_type*, Kokkos::DefaultExecutionSpace::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

        // Calculate shared memory size for the histogram (binCount elements)
        const size_t shared_size = scratch_view_type::shmem_size(binCount);

        const size_type team_size   = 128;
        const size_type block_size  = team_size * 8;
        const size_type num_leagues = (localNumParticles + block_size - 1) / block_size; // number of teams!
        
        // Set up team policy with scratch memory allocation for each team
        team_policy policy(num_leagues, team_size); 
        policy = policy.set_scratch_size(0, Kokkos::PerTeam(shared_size));

        // Launch a team parallel_for with the scratch memory setup
        Kokkos::parallel_for("initLocalHist", policy, KOKKOS_LAMBDA(const member_type& teamMember) {
            // Allocate team-local histogram in scratch memory
            scratch_view_type team_local_hist(teamMember.team_scratch(0), binCount);

            // Initialize shared memory histogram to zero
            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, binCount), [&](const bin_index_type b) {
                team_local_hist(b) = 0;
            });
            teamMember.team_barrier();

            const size_type start_i = teamMember.league_rank() * block_size;
            const size_type end_i   = Kokkos::min(start_i + block_size, localNumParticles);

            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, start_i, end_i), [&](const size_type i) {
                bin_index_type ndx = binIndex(i); // Get bin index for the particle
                if (ndx < binCount) Kokkos::atomic_increment(&team_local_hist(ndx)); // Atomic within shared memory
            });
            teamMember.team_barrier();

            // Reduce the team-local histogram into global memory
            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, binCount), [&](const bin_index_type i) {
                Kokkos::atomic_add(&device_histo(i), team_local_hist(i));
            });
        });

        localBinHisto_m.modify_device();
    }
    

    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::initLocalHisto(HistoReductionMode modePreference) {
        Inform msg("AdaptBins");
        
        bin_index_type binCount = getCurrentBinCount();

        // Determine the execution method based on the bin count and the mode...
        HistoReductionMode mode = determineHistoReductionMode(modePreference, binCount); 

        IpplTimings::startTimer(bExecuteHistoReductionT);
        if (mode == HistoReductionMode::HostOnly) {
            msg << "Using host-only parallel_reduce reduction." << endl;
            HostArrayReduction<size_type, bin_index_type>::binCountStatic = binCount;  
            HostArrayReduction<size_type, bin_index_type> reducer_arr;
            executeInitLocalHistoReduction(reducer_arr);
        } else if (mode == HistoReductionMode::TeamBased) {
            msg << "Using team-based + atomic reduction." << endl;
            executeInitLocalHistoReductionTeamFor();
        } else if (mode == HistoReductionMode::ParallelReduce) {
            auto to_reduce = createReductionObject<size_type, bin_index_type>(binCount);
            std::visit([&](auto& reducer_arr) {
                msg << "Starting parallel_reduce, array size = " << sizeof(reducer_arr.the_array) / sizeof(reducer_arr.the_array[0]) << endl;
                executeInitLocalHistoReduction(reducer_arr);
            }, to_reduce);
        } else {
            msg << "No valid execution method defined to initialize local histogram for energy binning." << endl;
            ippl::Comm->abort(); // Exit, since histogram cannot be built (should not happen...)!
        }
        IpplTimings::stopTimer(bExecuteHistoReductionT);

        msg << "Reducer ran without error." << endl;
        
        localBinHisto_m.sync(); // since all reductions happen on device --> sync changes if DualView mode is used
    }


    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::initGlobalHistogram() {
        Inform msg("AdaptBins");

        // Number of local bins == number of global bins!
        // Note that this step only happens after the uniform histogram is computed. 
        bin_index_type numBins = getCurrentBinCount();
        
        // Create a view to hold the global histogram on all ranks
        globalBinHisto_m = h_histo_type_g("globalBinHisto_m", numBins, xMax_m - xMin_m,
                                          binningAlpha_m, binningBeta_m, desiredWidth_m);
        
        // Get host ("mirror" <-- depends on mode of BinHisto class) view of histograms --> reduce on host!
        hview_type localBinHistoHost    = localBinHisto_m.template getHostView<hview_type>(localBinHisto_m.getHistogram()); 
        hview_type_g globalBinHistoHost = globalBinHisto_m.template getHostView<hview_type_g>(globalBinHisto_m.getHistogram()); 

        
        // Note: The allreduce also works when the .data() returns a CUDA space pointer.
        //       However, for some reason, copying manually to host and then allreduce-ing is faster. 
        IpplTimings::startTimer(bAllReduceGlobalHistoT);
        ippl::Comm->allreduce(
            localBinHistoHost.data(),           
            globalBinHistoHost.data(),
            numBins,                  
            std::plus<size_type>()
        );
        IpplTimings::stopTimer(bAllReduceGlobalHistoT);

        // The global histogram is currently on host --> sync to device if used as DualView.
        // Note that the global histogram should be used in a host-only mode, so these calls
        // don't have an effect!
        globalBinHisto_m.modify_host();
        globalBinHisto_m.sync(); 
 
        msg << "Global histogram created." << endl;
    }


    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::genAdaptiveHistogram() {
        var_selector_m.updateDataArr(bunch_m); // Probably not necessary, since it is called before updating particles; but just in case!
        bin_view_type binIndex                  = getBinView();
        hindex_transform_type adaptLookup       = globalBinHisto_m.mergeBins();

        // This step is currently unnecessary, since mergeBins works on host only. However, just in case, we can keep this
        // if should become available on device as well (there is barely a performance hit, since the lookup array is usually 
        // only 128 elements large).
        dindex_transform_type adaptLookupDevice = dindex_transform_type("adaptLookupDevice", currentBins_m);
        Kokkos::deep_copy(adaptLookupDevice, adaptLookup);
        
        setCurrentBinCount(globalBinHisto_m.getCurrentBinCount());

        // Map old indices to the new histogram ("Rebin")
        Kokkos::parallel_for("RebinParticles", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const size_type& i) {
            bin_index_type oldBin = binIndex(i);
            binIndex(i) = adaptLookupDevice(oldBin);
        });

        // Update local histogram with new indices
        instantiateHistogram(true);
        initLocalHisto(); // Runs reducer on new bin indices (also does the sync).

        localBinHisto_m.initPostSum();                   // Only init postsum, since the widths are not constant anymore.
        localBinHisto_m.copyBinWidths(globalBinHisto_m); // Manually copy non-constant widths from global histogram.
    }

    template <typename BunchType, typename BinningSelector>
    template <typename T, unsigned Dim>
    VField_t<T, Dim>& AdaptBins<BunchType, BinningSelector>::LTrans(VField_t<T, Dim>& field, const bin_index_type& currentBin) {
        Inform m("AdaptBins");
        
        position_view_type P = bunch_m->P.getView();
        hash_type indices    = sortedIndexArr_m;

        // Calculate gamma factor for field back transformation 
        // Note that P is saved normalized in OPAL, so technically p/mc
        Vector<T, Dim> gamma_bin2(0.0);
        Kokkos::parallel_reduce("CalculateGammaFactor", getBinIterationPolicy(currentBin),  
            KOKKOS_LAMBDA(const size_type& i, Vector<double, 3>& v) {
                Vector<double, 3> v_comp = P(indices(i)); 
                v                       += v_comp;  
            }, Kokkos::Sum<Vector<T, Dim>>(gamma_bin2));
        bin_index_type npart_bin = getNPartInBin(currentBin);
        
        /**
         * TODO Note: when the load balancer is not called often enough, then the adaptive bin
         * can lead to a phenomenon where the number of particles in a bin is zero, which leads to
         * a division by zero. 
         * So: either check if the number is 0 or make sure ranks always have enough particles!
         */
        gamma_bin2  = (npart_bin == 0) ? Vector<double, 3>(0.0) : gamma_bin2/npart_bin; // Now we have <P> for this bin
        gamma_bin2  = -sqrt(1.0 + gamma_bin2*gamma_bin2); // in these units: gamma=sqrt(1 + <P>^2), assuming <P^2>~0 (since bunch per bin should be "considered constant") // -1.0 / sqrt(1.0 - gamma_bin2 / c2); // negative sign, since we want the inverse transformation
        
        m << "Gamma(binIndex = " << currentBin << ") = -" << gamma_bin2 << endl;

        // Next apply the transformation --> do it manually, since fc->E*gamma does not exist in IPPL...
        ippl::parallel_for("TransformFieldWithVelocity", field.getFieldRangePolicy(), 
                           KOKKOS_LAMBDA(const ippl::RangePolicy<Dim>::index_array_type& idx) {
            apply(field, idx) *= gamma_bin2;
        });

        return field;
    }
    

    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::sortContainerByBin() {
        /*
         * Assume, this function is called after the prefix sum is initialized.
         * Then the particles need to be changed (sorted) in the right order and finally
         * the range_policy can simply be retrieved from the prefix sum for the scatter().
         * This is handled by the BinHisto class.
         */

        bin_view_type bins          = getBinView();
        size_type localNumParticles = bunch_m->getLocalNum();
        size_type numBins           = getCurrentBinCount();
        dview_type bin_counts       = localBinHisto_m.template getDeviceView<dview_type>(localBinHisto_m.getHistogram());

        // Get post sum (already calculated with histogram and saved inside local_bin_histo_post_sum_m), use copy to not modify the original
        Kokkos::View<size_type*> bin_offsets("bin_offsets", numBins + 1);
        Kokkos::deep_copy(bin_offsets, localBinHisto_m.template getDeviceView<dview_type>(localBinHisto_m.getPostSum()));

        // Initialize index array, use local shallow copy of the view
        sortedIndexArr_m  = hash_type("indices", localNumParticles);
        hash_type indices = sortedIndexArr_m;
        
        // Builds sorted index array using bin sort
        IpplTimings::startTimer(bSortContainerByBinT);
        Kokkos::parallel_for("SortIndices", localNumParticles, KOKKOS_LAMBDA(const size_type& i) {
            size_type target_bin = bins(i);
            size_type target_pos = Kokkos::atomic_fetch_add(&bin_offsets(target_bin), 1);

            // Place the current particle directly into its target position
            indices(target_pos) = i;
        });
        IpplTimings::stopTimer(bSortContainerByBinT);

        // Test is sorting was successful (only in debug mode)
        #ifdef DEBUG
        IpplTimings::startTimer(bVerifySortingT);
        if (localNumParticles > 1 && !viewIsSorted<bin_index_type>(getBinView(), indices, localNumParticles)) {
            Inform msg("AdaptBins");
            msg << "Sorting failed." << endl;
            ippl::Comm->abort();
        } 
        IpplTimings::stopTimer(bVerifySortingT);
        #endif
    }

}

#endif // ADAPT_BINS_HPP


