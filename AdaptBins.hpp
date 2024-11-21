#ifndef ADAPT_BINS_HPP
#define ADAPT_BINS_HPP

// #include "AdaptBins.h"

namespace ParticleBinning {

    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::initLimits() {
        Inform msg("AdaptBins");  // INFORM_ALL_NODES

        Kokkos::MinMaxScalar<value_type> localMinMax;
        // position_view_type localData = bunch_m->R.getView();
        var_selector_m.updateDataArr(bunch_m); // update needed if bunch->create() is called between binnings!
        BinningSelector var_selector = var_selector_m;  
        
        static IpplTimings::TimerRef histoLimits = IpplTimings::getTimer("initHistoLimits");
        IpplTimings::startTimer(histoLimits);

        Kokkos::parallel_reduce("localBinLimitReduction", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const size_type i, Kokkos::MinMaxScalar<value_type>& update) {
            value_type val = var_selector(i); // localData(i)[2]; // use z axis for binning!
            update.min_val = Kokkos::min(update.min_val, val);
            update.max_val = Kokkos::max(update.max_val, val);
        }, Kokkos::MinMax<value_type>(localMinMax));
        xMin_m = localMinMax.min_val;
        xMax_m = localMinMax.max_val;

        // Putting the same to-reduce variable as an argument ensures that every node gets the correct min/max and not just the root node!
        // Note: boradcast does not exist, use allreduce for reduce+broadcast together!
        ippl::Comm->allreduce(xMax_m, 1, std::greater<value_type>());
        ippl::Comm->allreduce(xMin_m, 1, std::less<value_type>());

        IpplTimings::stopTimer(histoLimits);

        binWidth_m = (xMax_m - xMin_m) / currentBins_m;
        
        msg << "Initialized limits. Min: " << xMin_m << ", max: " << xMax_m << ", binWidth: " << binWidth_m << endl;
    }

    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::instantiateHistogram(bool setToZero) {
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
        bin_view_type binIndex        = bunch_m->bin.getView();  

        // Declare the variables locally before the Kokkos::parallel_for (to avoid implicit this capture in Kokkos lambda)
        value_type xMin = xMin_m, xMax = xMax_m, binWidthInv = 1.0/binWidth_m;
        bin_index_type numBins = currentBins_m;
        // Alternatively explicit capture: [xMin = xMin_m, xMax = xMax_m, binWidth = binWidth_m, numBins = currentBins_m, localData = localData, binIndex = binIndex]

        static IpplTimings::TimerRef assignParticleBins = IpplTimings::getTimer("assignParticleBins");
        IpplTimings::startTimer(assignParticleBins);

        Kokkos::parallel_for("assignParticleBins", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const size_type i) {
                // Access the z-axis position of the i-th particle
                value_type v = var_selector(i); // localData(i)[2];  
                
                // Assign the bin index to the particle (directly on device)
                bin_index_type bin = getBin(v, xMin, xMax, binWidthInv, numBins);
                binIndex(i)        = bin;
        });
        IpplTimings::stopTimer(assignParticleBins);
        msg << "All bins assigned." << endl; 
    }

    template<typename BunchType, typename BinningSelector>
    template<typename ReducerType>
    void AdaptBins<BunchType, BinningSelector>::executeInitLocalHistoReduction(ReducerType& to_reduce) {
        bin_view_type binIndex        = bunch_m->bin.getView();  
        bin_histo_type localBinHisto  = localBinHisto_m;
        bin_index_type binCount       = getCurrentBinCount();

        static IpplTimings::TimerRef initLocalHisto = IpplTimings::getTimer("initLocalHistoParallelReduce");
        IpplTimings::startTimer(initLocalHisto);
        Kokkos::parallel_reduce("initLocalHist", bunch_m->getLocalNum(), 
            KOKKOS_LAMBDA(const size_type& i, ReducerType& update) {
                bin_index_type ndx = binIndex(i);  // Determine the bin index for this particle
                update.the_array[ndx]++;           // Increment the corresponding bin count in the reduction array
            }, Kokkos::Sum<ReducerType>(to_reduce)
        );
        IpplTimings::stopTimer(initLocalHisto);

        // Copy the reduced results to the final histogram
        Kokkos::parallel_for("finalize_histogram", binCount, 
            KOKKOS_LAMBDA(const bin_index_type& i) {
                localBinHisto(i) = to_reduce.the_array[i];
            }
        );
    }

    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::executeInitLocalHistoReductionTeamFor() {
        bin_view_type binIndex            = bunch_m->bin.getView();
        bin_histo_type localBinHisto      = localBinHisto_m;
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

        static IpplTimings::TimerRef initLocalHisto = IpplTimings::getTimer("initLocalHistoTeamBased");
        IpplTimings::startTimer(initLocalHisto);

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
                Kokkos::atomic_add(&localBinHisto(i), team_local_hist(i));
            });
        });
        
        IpplTimings::stopTimer(initLocalHisto);
    }

    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::initLocalHisto(HistoReductionMode modePreference) {
        Inform msg("AdaptBins");
        
        bin_index_type binCount = getCurrentBinCount();

        // Determine the execution method based on the bin count and the mode...
        HistoReductionMode mode = determineHistoReductionMode(modePreference, binCount);

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

        msg << "Reducer ran without error." << endl;
        
        localBinHistoHost_m = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), localBinHisto_m);
    }

    template <typename BunchType, typename BinningSelector>
    void AdaptBins<BunchType, BinningSelector>::initGlobalHistogram() {
        Inform msg("GetGlobalHistogram");

        // Get the current number of bins
        bin_index_type numBins = getCurrentBinCount(); // number of local bins = number of global bins!
        
        // Create a view to hold the global histogram on all ranks
        bin_host_histo_type globalBinHisto("globalBinHistoHost", numBins);

        static IpplTimings::TimerRef globalHistoReduce = IpplTimings::getTimer("allReduceGlobalHisto");
        IpplTimings::startTimer(globalHistoReduce);
        
        // Need host mirror, otherwise the data is not available when the histogram is created using CUDA
        auto localBinHistoHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), localBinHisto_m);
        
        /*
         * Note: The allreduce also works when the .data() returns a CUDA space pointer.
         *       However, for some reason, copying manually to host and then allreducing is faster. 
         */
        ippl::Comm->allreduce(
            localBinHistoHost.data(),           // Pointer to local data
            globalBinHisto.data(),              // Pointer to global data
            numBins,                            // Number of elements to reduce
            std::plus<size_type>()              // Reduction operation
        );
        IpplTimings::stopTimer(globalHistoReduce);

        // The global histogram is currently on host, but can be saved on device
        globalBinHisto_m = bin_histo_type("globalBinHisto", numBins);
        Kokkos::deep_copy(globalBinHisto_m, globalBinHisto); // Copy to device view

        msg << "Global histogram created." << endl;
    }

    template <typename BunchType, typename BinningSelector>
    template <typename T, unsigned Dim>
    VField_t<T, Dim>& AdaptBins<BunchType, BinningSelector>::LTrans(VField_t<T, Dim>& field) {
        bin_view_type binIndex            = bunch_m->bin.getView();
        const size_type localNumParticles = bunch_m->getLocalNum(); 
        position_view_type P              = bunch_m->P.getView();

        // Calculate gamma factor for field back transformation --> TODO: change iteration if decide to use sorted particles!
        Vector<T, Dim> gamma_bin(0.0);
        Kokkos::parallel_reduce("CalculateGammaFactor", localNumParticles, 
            KOKKOS_LAMBDA(const size_type& i, Vector<double, 3>& local_sum_speed) {
                Vector<double, 3> v_comp = P(i); 
                local_sum_speed         += v_comp.dot(v_comp) * (binIndex(i) == i); 
            }, Kokkos::Sum<Vector<T, Dim>>(gamma_bin));
        gamma_bin = 1.0 / sqrt(1.0 + gamma_bin);
        std::cout << "Gamma factor calculated = " << gamma_bin << std::endl;

        // Next apply the transformation --> do it manually, since fc->E*gamma does not exist in IPPL...
        ippl::parallel_for("TransformFieldWithVelocity", field.getFieldRangePolicy(), 
            KOKKOS_LAMBDA(const ippl::RangePolicy<Dim>::index_array_type& idx) {
                apply(field, idx) *= gamma_bin;
            });

        return field;
    }

}

#endif // ADAPT_BINS_HPP


