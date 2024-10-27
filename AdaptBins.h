/**
 * @file AdaptBins.h
 * @brief Defines a structure to hold particles in energy bins and their associated data.
 * 
 * "AdaptBins" is a improved version of "PartBins" from the old OPAL. AdaptBins uses
 * data structures from IPPL and Kokkos to make improve performance when using MPI and Cuda.
 * In contrast to the old PartBin, this class allows for the use of rebinning during runtime. 
 */

/**
 * @copyright Copyright (c) 2007-2020, Paul Scherrer Institut, Villigen PSI, Switzerland
 * All rights reserved
 *
 * This file is part of OPALX.
 *
 * OPAL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU General Public License
 * along with OPAL. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef ADAPT_BINS_H
#define ADAPT_BINS_H

#include <memory>

#include <iostream>
#include "Ippl.h"
#include <Kokkos_ScatterView.hpp>



// To save the bin count statically
//template<typename SizeType, typename IndexType>
//IndexType ViewReducer<SizeType, IndexType>::numBins = 0;

template <typename BunchType>
class AdaptBins {
public:
    using value_type             = typename BunchType::Layout_t::value_type;
    using particle_position_type = typename BunchType::particle_position_type;
    using position_view_type     = typename particle_position_type::view_type;
    using size_type              = typename BunchType::size_type;
    //using particle_id_type = typename BunchType::particle_index_type;
    //using id_vector_type = typename BunchType::vector_type<size_type, 1>;  // PLayout type with T=size_type and dim=3
    //using value_view_type        = typename ippl::detail::ViewType<ippl::Vector<double, Dim>, 1>::view_type;

    // New bin_type representing the bin index data type
    using bin_index_type         = typename BunchType::bin_index_type; // has to be initialized over there...
    using bin_type               = typename ippl::ParticleAttrib<bin_index_type>;
    using bin_view_type          = typename bin_type::view_type; 

    using bin_histo_type         = typename Kokkos::View<size_type*>; // ippl::Vector<size_type*, 1>;
    using bin_host_histo_type    = typename Kokkos::View<size_type*, Kokkos::HostSpace>; // ippl::Vector<size_type*, 1>;
    
    //using reducer_type =  HistogramReducer<size_type, bin_index_type>;
    //using reducer_type = ViewReducer<size_type, bin_index_type>;


    //template<typename SizeType, typename IndexType>
    //bin_index_type reducer_type::numBins = 0; // save it here once to use it as static

    //using bin_histo_view_type    = typename bin_histo_type::view_type;

    std::shared_ptr<BunchType> bunch_m;  // Shared pointer to the particle container
    bin_index_type maxBins_m;  // Maximum number of bins
    bin_index_type currentBins_m; // Current number of bins
    value_type xMin_m;
    value_type xMax_m;
    value_type binWidth_m;
    bin_histo_type localBinHisto_m;

    // Constructor to initialize with a shared pointer to the particle container
    AdaptBins(std::shared_ptr<BunchType> bunch, bin_index_type maxBins = 10)
        : bunch_m(bunch)
        , maxBins_m(maxBins) {
        
        // Register the bin attribute in the particle bunch
        // bunch_m->addAttribute(bin_type);

        // * (bin == 4) ? q : 0

        currentBins_m = maxBins; // TODO for now...
    }

    // Function to initialize the bin limits (can be extended as needed)
    void initLimits() {
        Inform msg("AdaptBins");  // INFORM_ALL_NODES

        Kokkos::MinMaxScalar<value_type> localMinMax;
        //localMinMax.min_val = std::numeric_limits<value_type>::max();
        //localMinMax.max_val = std::numeric_limits<value_type>::lowest();
    
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

    // Get the current number of bins
    bin_index_type getCurrentBinCount() const {
        return currentBins_m;
    }

    void setCurrentBinCount(bin_index_type nBins) {
        // Set the number of current bins to the new value!
        currentBins_m = (nBins > maxBins_m) ? maxBins_m : nBins; 
        binWidth_m    = (xMax_m - xMin_m) / currentBins_m; // assuming particles did not change!
    }

    void initializeHistogram(bool setToZero = false) {
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

    KOKKOS_INLINE_FUNCTION
    static bin_index_type getBin(value_type x, value_type xMin, value_type xMax, value_type binWidthInv, bin_index_type numBins) {
        // Explanation: Don't access xMin, binWidth, ... through the members to avoid implicit
        // variable capture by Kokkos and potential copying overhead. Instead, pass them as an 
        // argument, s.t. Kokkos can capture them explicitly!
        // Make it static to avoid implicit capture of this inside Kokkos lambda! 

        // Ensure x is within bounds (clamp it between xMin and xMax --> this is only for bin assignment)
        // x = (x < xMin) ? xMin : ((x > xMax) ? xMax : x);
        x += (x < xMin) * (xMin - x) + (x > xMax) * (xMax - x); // puts x in the bin or nearest bin if out of bounds
        
        // Compute the bin index (ensuring it's within [0, numBins-1])
        //bin_index_type bin = (x - xMin) / binWidth;
        //return (bin >= numBins) ? (numBins - 1) : bin;  // Clamp to the maximum bin

        bin_index_type bin = (x - xMin) * binWidthInv; // multiply with inverse of binwidth
        return (bin >= numBins) ? (numBins - 1) : bin;  // Clamp to the maximum bin
    }

    void assignBinsToParticles() {
        // Set the bin attribute for the given particle
        //bunch_m.bin(particleIndex) = bin;
        Inform msg("AdaptBins");

        position_view_type localData = bunch_m->R.getView();
        bin_view_type binIndex       = bunch_m->bin.getView();  
        // Use Kokkos::parallel_for to iterate over the particles and assign them to bins
        /*Kokkos::parallel_for("assignParticleBins", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const int i) {
            value_type v = localData(i)[2];  // Assuming z-axis is the third dimension
            binIndex(i) = getBin(v, xMin, xMax, binWidth, numBins);
        });*/
        
        // This means that new particles were added, not that there was a rebin (yet!)

        // Declare the variables locally before the Kokkos::parallel_for (to avoid implicit this capture in Kokkos lambda)
        value_type xMin = xMin_m, xMax = xMax_m, binWidthInv = 1.0/binWidth_m;
        bin_index_type numBins = currentBins_m;
        // Alternatively explicit capture: [xMin = xMin_m, xMax = xMax_m, binWidth = binWidth_m, numBins = currentBins_m, localData = localData, binIndex = binIndex]

        // Capture variables explicitly...
        static IpplTimings::TimerRef assignParticleBins = IpplTimings::getTimer("assignParticleBins");
        IpplTimings::startTimer(assignParticleBins);
        //initializeHistogram();
        //bin_histo_type localBinHisto = localBinHisto_m; // use copy of histo... view to avoid implicit "this" capture in Kokkos:parallel_for

        Kokkos::parallel_for("assignParticleBins", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const size_type i) {
                // Access the z-axis position of the i-th particle
                value_type v = localData(i)[2];  
                
                // Assign the bin index to the particle (directly on device)
                bin_index_type bin = getBin(v, xMin, xMax, binWidthInv, numBins);
                binIndex(i)        = bin;
                //localBinHisto(bin) += 1;
                //Kokkos::atomic_fetch_add(&localBinHisto(bin), static_cast<bin_index_type>(1)); // since "Kokkos::atomic_add" might not support unsigned long, only int?!
                //Kokkos::atomic_add(&localBinHisto(bin), 1);
                //Kokkos::atomic_increment(&localBinHisto(bin)); // Not needed with parallel_reduce
        });
        IpplTimings::stopTimer(assignParticleBins);
        msg << "All bins assigned." << endl; 

        initLocalHisto();
        msg << "Local Histogram initialized." << endl;
        // Ensure all bins are assigned before further calculations!
        // Kokkos::fence();   
    }

    void initLocalHisto() {
        Inform msg("AdaptBins");
        initializeHistogram(false); // Init histogram (no need to set to 0, since contribute_into overwrites values...)
        msg << "Histogram initialized to 0" << endl;

        bin_view_type binIndex       = bunch_m->bin.getView();
        bin_histo_type localBinHisto = localBinHisto_m;
        bin_index_type binCount      = getCurrentBinCount();
        
        msg << "Starting reducer." << endl;

        static IpplTimings::TimerRef initLocalHisto = IpplTimings::getTimer("initLocalHisto");
        IpplTimings::startTimer(initLocalHisto);

        Kokkos::Experimental::ScatterView<size_type*> scatter_view("scatter_view", binCount);
        Kokkos::parallel_for("initLocalHist", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const size_type i) {
                auto scatter = scatter_view.access();
                bin_index_type ndx = binIndex(i);
                scatter(ndx)++;
            });

        msg << "Reduced, now contributing." << endl;
        scatter_view.contribute_into(localBinHisto);
        
        IpplTimings::stopTimer(initLocalHisto);

        msg << "Reducer ran without error." << endl;
    }

    void doFullRebin(bin_index_type nBins) {
        // Sets number of bins, assigns bin numbers and initializes histogram
        setCurrentBinCount(nBins);
        assignBinsToParticles();
    }

    bin_host_histo_type getGlobalHistogram() {
        Inform msg("GetGlobalHistogram");
        // Get the current number of bins
        bin_index_type numBins = getCurrentBinCount(); // number of local bins = number of global bins!
        //msg << "1" << endl;
        // Create a view to hold the global histogram on all ranks
        bin_host_histo_type globalBinHisto("globalBinHisto", numBins);

        // Need host mirror, otherwise the data is not available when the histogram is created using CUDA
        auto localBinHistoHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), localBinHisto_m);
        // Kokkos::deep_copy(localBinHistoHost, localBinHisto_m);
        
        // Loop through each bin and reduce its value across all nodes
        for (bin_index_type i = 0; i < numBins; ++i) {
            //msg << "3 - " << i << endl;
            size_type localValue = localBinHistoHost(i);
            size_type globalValue = 0;
            //msg << "4 - " << i << endl;
            // Perform the reduction to sum the values across all nodes
            ippl::Comm->allreduce(localValue, globalValue, 1, std::plus<size_type>());  // Reduce to root node
            //msg << "4 - " << i << endl;
            // Store the global value only on rank 0
            //msg << "5 - " << i << endl;
            globalBinHisto(i) = globalValue; // Saves global histogram on ALL ranks...
            //msg << "5 - " << i << endl;
        }

        // Return the global histogram (note: only meaningful on rank 0)
        return globalBinHisto;
    }

    void print() {
        Inform os("PHisto");
        // Only works correct for rank 0
        os << "-----------------------------------------" << endl;
        os << "     Output Global Binning Structure     " << endl;

        bin_index_type numBins = getCurrentBinCount();
        os << "Bins = " << numBins << " hBin = " << binWidth_m << endl;
        os << "Bin #;Val" << endl;

        // Get the global histogram (reduced across all nodes)
        bin_host_histo_type globalBinHisto = getGlobalHistogram();

        // Only rank 0 prints the global histogram
        size_type total = 0;
        for (bin_index_type i = 0; i < numBins; ++i) {
            size_type val = globalBinHisto(i);
            os << i << ";" << val << endl;
            total += val;
        }   
        os << "Total = " << total << endl;
        os << "-----------------------------------------" << endl;
    }

    void debug() {
        Inform msg("KOKKOS DEBUG", INFORM_ALL_NODES);

        int rank = ippl::Comm->rank();

        // Check number of CPU threads
        #ifdef KOKKOS_ENABLE_OPENMP
        int num_threads = Kokkos::OpenMP::concurrency();
        msg << "Number of CPU threads: " << num_threads << endl;
        #endif  

        // Check number of GPUs (CUDA devices)
        #ifdef KOKKOS_ENABLE_CUDA
        int num_gpus = Kokkos::Cuda::detect_device_count();
        msg << "Rank " << rank << " sees " << num_gpus << " GPUs available." << endl;
        #else
        msg << "Rank " << rank << ": Kokkos is not using CUDA (GPU support disabled)." << endl;
        #endif
    
    }
};



#endif 


