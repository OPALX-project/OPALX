#ifndef ADAPT_BINS_HPP
#define ADAPT_BINS_HPP

#include "AdaptBins.h"

namespace ParticleBinning {

    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::initLimits() {
        Inform msg("AdaptBins");  // INFORM_ALL_NODES
            
        //static IpplTimings::TimerRef histoLimits = IpplTimings::getTimer("initHistoLimits");
        //IpplTimings::startTimer(histoLimits);

        var_selector_m.updateDataArr(bunch_m); // update needed if bunch->create() is called between binnings!
        BinningSelector var_selector = var_selector_m;  
        size_type nlocal             = bunch_m->getLocalNum();

        IpplTimings::startTimer(bInitLimitsT);
        if (nlocal <= 0) {
            msg << "Particles in the bunch = " << nlocal << ". Overwriting limits manually." << endl;
            xMin_m = xMax_m = (nlocal == 0) ? 0 : 0; // var_selector(0); 
        } else {
            Kokkos::MinMaxScalar<value_type> localMinMax;
            // Sadly this is necessary, since Kokkos seems to have a problem when nlocal == 1 where it does not update localMinMax...
            if (nlocal == 1) {
                //std::cout << "-------------->>>> " << var_selector(0) << " - " << bunch_m->R(0) << " - " << bunch_m->P(0) << std::endl;
                Kokkos::View<value_type, Kokkos::HostSpace> host_scalar("host_scalar"); 
                Kokkos::View<value_type> tmp_dvalue("tmp_dvalue");
                Kokkos::parallel_for("tmp_dvalue", 1, KOKKOS_LAMBDA(const size_type) { tmp_dvalue() = var_selector(0); });
                Kokkos::deep_copy(host_scalar, tmp_dvalue);
                localMinMax.max_val = localMinMax.min_val = host_scalar();
                //std::cout << "Meh " << host_scalar() << std::endl;
            } else {
                Kokkos::parallel_reduce("localBinLimitReduction", nlocal, KOKKOS_LAMBDA(const size_type i, Kokkos::MinMaxScalar<value_type>& update) {
                    value_type val = var_selector(i); // localData(i)[2]; // use z axis for binning!
                    update.min_val = Kokkos::min(update.min_val, val);
                    update.max_val = Kokkos::max(update.max_val, val);
                }, Kokkos::MinMax<value_type>(localMinMax));
            }
            xMin_m = localMinMax.min_val;
            xMax_m = localMinMax.max_val;
        }
        IpplTimings::stopTimer(bInitLimitsT);

        // Putting the same to-reduce variable as an argument ensures that every node gets the correct min/max and not just the root node!
        // Note: boradcast does not exist, use allreduce for reduce+broadcast together!
        IpplTimings::startTimer(bAllReduceLimitsT);
        ippl::Comm->allreduce(xMax_m, 1, std::greater<value_type>());
        ippl::Comm->allreduce(xMin_m, 1, std::less<value_type>());
        IpplTimings::stopTimer(bAllReduceLimitsT);

        //IpplTimings::stopTimer(histoLimits);

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
            /*Kokkos::parallel_for("initHistogram", numBins, KOKKOS_LAMBDA(const bin_index_type i) {
                device_histo(i) = 0;
            });*/
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

        // Ensure x is within bounds (clamp it between xMin and xMax --> this is only for bin assignment)
        // x = (x < xMin) ? xMin : ((x > xMax) ? xMax : x);
        x += (x < xMin) * (xMin - x) + (x > xMax) * (xMax - x); // puts x in the bin or nearest bin if out of bounds

        bin_index_type bin = (x - xMin) * binWidthInv; // multiply with inverse of binwidth
        return (bin >= numBins) ? (numBins - 1) : bin;  // Clamp to the maximum bin
    }

    template <typename BunchType, typename BinningSelector>    
    void AdaptBins<BunchType, BinningSelector>::assignBinsToParticles() {
        // Set the bin attribute for the given particle
        Inform msg("AdaptBins");

        // position_view_type localData = bunch_m->R.getView();
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
        // Alternatively explicit capture: [xMin = xMin_m, xMax = xMax_m, binWidth = binWidth_m, numBins = currentBins_m, localData = localData, binIndex = binIndex]

        //static IpplTimings::TimerRef assignParticleBins = IpplTimings::getTimer("assignParticleBins");
        //IpplTimings::startTimer(assignParticleBins);
        Kokkos::parallel_for("assignParticleBinsConst", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const size_type& i) {
                // Access the z-axis position of the i-th particle
                value_type v = var_selector(i); // localData(i)[2];  
                
                // Assign the bin index to the particle (directly on device)
                bin_index_type bin = getBin(v, xMin, xMax, binWidthInv, numBins);
                binIndex(i)        = bin;
        });
        IpplTimings::stopTimer(bAssignUniformBinsT);
        //IpplTimings::stopTimer(assignParticleBins);
        msg << "All bins assigned." << endl; 
    }

    template<typename BunchType, typename BinningSelector>
    template<typename ReducerType>
    void AdaptBins<BunchType, BinningSelector>::executeInitLocalHistoReduction(ReducerType& to_reduce) {
        bin_view_type binIndex        = getBinView(); 
        dview_type device_histo       = localBinHisto_m.template getDeviceView<dview_type>(localBinHisto_m.getHistogram());
        bin_index_type binCount       = getCurrentBinCount();

        //static IpplTimings::TimerRef initLocalHisto = IpplTimings::getTimer("initLocalHistoParallelReduce");
        //IpplTimings::startTimer(initLocalHisto);
        
        //auto start = std::chrono::high_resolution_clock::now(); // TODO: remove

        Kokkos::parallel_reduce("initLocalHist", bunch_m->getLocalNum(), 
            KOKKOS_LAMBDA(const size_type& i, ReducerType& update) {
                bin_index_type ndx = binIndex(i);  // Determine the bin index for this particle
                update.the_array[ndx]++;           // Increment the corresponding bin count in the reduction array
            }, Kokkos::Sum<ReducerType>(to_reduce)
        );
        //IpplTimings::stopTimer(initLocalHisto);

        //auto end = std::chrono::high_resolution_clock::now(); // TODO: remove
        //long long duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); // TODO: remove
        //std::cout << "executeInitLocalHistoReduction;" << bunch_m->getLocalNum() << ";" << binCount << ";" << duration << std::endl; // TODO: remove

        // Copy the reduced results to the final histogram
        Kokkos::parallel_for("finalize_histogram", binCount, 
            KOKKOS_LAMBDA(const bin_index_type& i) {
                device_histo(i) = to_reduce.the_array[i];
            }
        );
        

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

        //static IpplTimings::TimerRef initLocalHisto = IpplTimings::getTimer("initLocalHistoTeamBased");
        //IpplTimings::startTimer(initLocalHisto);
        //auto start = std::chrono::high_resolution_clock::now(); // TODO: remove

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
                if (ndx < binCount) Kokkos::atomic_increment(&team_local_hist(ndx)); // Kokkos::atomic_fetch_add(&team_local_hist(ndx), 1); // Atomic within shared memory
            });
            teamMember.team_barrier();

            // Reduce the team-local histogram into global memory
            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, binCount), [&](const bin_index_type i) {
                Kokkos::atomic_add(&device_histo(i), team_local_hist(i));
            });
        });

        //auto end = std::chrono::high_resolution_clock::now(); // TODO: remove
        //long long duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); // TODO: remove
        //std::cout << "executeInitLocalHistoReductionTeamFor;" << bunch_m->getLocalNum() << ";" << binCount << ";" << duration << std::endl; // TODO: remove
        
        //IpplTimings::stopTimer(initLocalHisto);
        localBinHisto_m.modify_device();
    }

    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::initLocalHisto(HistoReductionMode modePreference) {
        Inform msg("AdaptBins");
        
        bin_index_type binCount = getCurrentBinCount();

        // Determine the execution method based on the bin count and the mode...
        HistoReductionMode mode = determineHistoReductionMode(modePreference, binCount); // modePreference;

        IpplTimings::startTimer(bExecuteHistoReductionT);
        if (mode == HistoReductionMode::HostOnly) {
            msg << "Using host-only parallel_reduce reduction." << endl;
            HostArrayReduction<size_type, bin_index_type>::binCountStatic = binCount; // set size of the histogram 
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
            ippl::Comm->abort(); // Exit, since error!
        }
        IpplTimings::stopTimer(bExecuteHistoReductionT);

        msg << "Reducer ran without error." << endl;
        
        localBinHisto_m.sync(); // since all reductions happen on device --> marked as modified 
    }

    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::initGlobalHistogram() {
        Inform msg("AdaptBins");

        // Get the current number of bins
        bin_index_type numBins = getCurrentBinCount(); // number of local bins = number of global bins!
        
        // Create a view to hold the global histogram on all ranks
        //bin_host_histo_type globalBinHisto("globalBinHistoHost", numBins);
        globalBinHisto_m = h_histo_type_g("globalBinHisto_m", numBins, xMax_m - xMin_m,
                                          binningAlpha_m, binningBeta_m, desiredWidth_m);

        
        // Need host mirror, otherwise the data is not available when the histogram is created using CUDA
        hview_type localBinHistoHost    = localBinHisto_m.template getHostView<hview_type>(localBinHisto_m.getHistogram()); 
        hview_type_g globalBinHistoHost = globalBinHisto_m.template getHostView<hview_type_g>(globalBinHisto_m.getHistogram()); 

        //static IpplTimings::TimerRef globalHistoReduce = IpplTimings::getTimer("allReduceGlobalHisto");
        //IpplTimings::startTimer(globalHistoReduce);

        /*
         * Note: The allreduce also works when the .data() returns a CUDA space pointer.
         *       However, for some reason, copying manually to host and then allreducing is faster. 
         */
        IpplTimings::startTimer(bAllReduceGlobalHistoT);
        ippl::Comm->allreduce(
            localBinHistoHost.data(),           // Pointer to local data
            globalBinHistoHost.data(),              // Pointer to global data
            numBins,                            // Number of elements to reduce
            std::plus<size_type>()              // Reduction operation
        );
        IpplTimings::stopTimer(bAllReduceGlobalHistoT);
        //IpplTimings::stopTimer(globalHistoReduce);

        // The global histogram is currently on host, but can be saved on device
        globalBinHisto_m.modify_host();
        globalBinHisto_m.sync(); 
        //globalBinHisto_m.init(); // syncs and inits the initial postSum/widths array
 
        msg << "Global histogram created." << endl;
    }

    template <typename BunchType, typename BinningSelector>
    template <typename T, unsigned Dim>
    VField_t<T, Dim>& AdaptBins<BunchType, BinningSelector>::LTrans(VField_t<T, Dim>& field, const bin_index_type& currentBin) {
        Inform m("AdaptBins");
        //bin_view_type binIndex            = getBinView();
        //const size_type localNumParticles = bunch_m->getLocalNum(); 
        position_view_type P = bunch_m->P.getView();
        hash_type indices    = sortedIndexArr_m;

        // TODO: remove once in OPAL, since it shoud already exist over there!
        // constexpr double c2 = 299792458.0*299792458.0; // Speed of light in m/s
        // Note: not needed, since P is a normalized value "p/mc"

        // Calculate gamma factor for field back transformation --> TODO: change iteration if decide to use sorted particles!
        // Note that P is saved normalized in OPAL, so technically p/mc
        Vector<T, Dim> gamma_bin2(0.0);
        Kokkos::parallel_reduce("CalculateGammaFactor", getBinIterationPolicy(currentBin), // localNumParticles, 
            KOKKOS_LAMBDA(const size_type& i, Vector<double, 3>& v) {
                Vector<double, 3> v_comp = P(indices(i)); 
                v                       += v_comp; // like this for elemntwise multiplication (not .dot...) // v_comp.dot(v_comp) * (binIndex(i) == currentBin); 
            }, Kokkos::Sum<Vector<T, Dim>>(gamma_bin2));
        bin_index_type npart_bin = getNPartInBin(currentBin);
        /**
         * 
         * TODO Note: when the load balancer is not called often enough, then the adaptive bin
         * can lead to a phenomenon where the number of particles in a bin is zero, which leads to
         * a division by zero. 
         * So: either check if the number is 0 or make sure ranks always have enough particles!
         */

        gamma_bin2  = (npart_bin == 0) ? Vector<double, 3>(0.0) : gamma_bin2/npart_bin; // Now we have <P> for this bin
        gamma_bin2  = -sqrt(1.0 + gamma_bin2*gamma_bin2); // in these units: gamma=sqrt(1 + <P>^2), assuming <P^2>~0 (since bunch per bin should be "considered constant") // -1.0 / sqrt(1.0 - gamma_bin2 / c2); // negative sign, since we want the inverse transformation
        // std::cout << "Gamma factor calculated = " << gamma_bin2 << std::endl;
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
         */
        Inform msg("AdaptBins");

        //static IpplTimings::TimerRef argSortBins      = IpplTimings::getTimer("argSortBins");
        //static IpplTimings::TimerRef permutationTimer = IpplTimings::getTimer("sortPermutationTimer");
        //static IpplTimings::TimerRef isSortedCheck    = IpplTimings::getTimer("isSortedCheck");
        //static IpplTimings::TimerRef binSortingAndScatterT = IpplTimings::getTimer("binSortingAndScatter");

        bin_view_type bins          = getBinView();
        size_type localNumParticles = bunch_m->getLocalNum();
        size_type numBins           = getCurrentBinCount();
        dview_type bin_counts       = localBinHisto_m.template getDeviceView<dview_type>(localBinHisto_m.getHistogram());

        //IpplTimings::startTimer(argSortBins);
        // Get post sum (already calculated with histogram and saved inside local_bin_histo_post_sum_m), use copy to not modify the original
        Kokkos::View<size_type*> bin_offsets("bin_offsets", numBins + 1);
        // typename d_histo_type::dview_type postSumView = localBinHisto_m.view_device(HistoTypeIdentifier::PostSum);
        Kokkos::deep_copy(bin_offsets, localBinHisto_m.template getDeviceView<dview_type>(localBinHisto_m.getPostSum()));

        sortedIndexArr_m  = hash_type("indices", localNumParticles);
        hash_type indices = sortedIndexArr_m;
        
        /*
        TODO: maybe change value_type of hash_type to size_type instead int of at some point???
        */
        //IpplTimings::startTimer(binSortingAndScatterT);
        //auto start = std::chrono::high_resolution_clock::now(); // TODO: remove
        IpplTimings::startTimer(bSortContainerByBinT);
        Kokkos::parallel_for("InPlaceSortIndices", localNumParticles, KOKKOS_LAMBDA(const size_type& i) {
            size_type target_bin = bins(i);
            size_type target_pos = Kokkos::atomic_fetch_add(&bin_offsets(target_bin), 1);

            // Place the current particle directly into its target position
            indices(target_pos) = i;
        });
        IpplTimings::stopTimer(bSortContainerByBinT);
        //auto end = std::chrono::high_resolution_clock::now(); // TODO: remove
        //long long duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); // TODO: remove
        //std::cout << "sortIndices;" << localNumParticles << ";" << numBins << ";" << duration << std::endl; // TODO: remove
        
        //IpplTimings::stopTimer(argSortBins);
        //msg << "Argsort on bin index completed." << endl;
        //Kokkos::fence();

        /*IpplTimings::startTimer(permutationTimer);
        bunch_m->template forAllAttributes([&]<typename Attribute>(Attribute*& attribute) {
            using memory_space    = typename Attribute::memory_space;

            // Ensure indices are in the correct memory space --> copies data ONLY when different memory spaces, so should be efficient
            auto indices_device = Kokkos::create_mirror_view_and_copy(memory_space{}, indices);

            attribute->pack(indices_device);
            attribute->unpack(localNumParticles, true);
        });
        IpplTimings::stopTimer(permutationTimer);*/
        //IpplTimings::stopTimer(binSortingAndScatterT);
        //msg << "Permutation of particle attributes completed." << endl;

        // TODO: remove, just for testing purposes (get new bin view, since the old memory address might be overwritten by this action...)
        IpplTimings::startTimer(bVerifySortingT);
        if (localNumParticles > 1 && !viewIsSorted<bin_index_type>(getBinView(), indices, localNumParticles)) {
            msg << "Sorting failed." << endl;
            ippl::Comm->abort();
        } 
        IpplTimings::stopTimer(bVerifySortingT);
    }

}

#endif // ADAPT_BINS_HPP


