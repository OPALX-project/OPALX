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

//#include <memory>
//#include <iostream>
#include "Ippl.h"
#include <Kokkos_ScatterView.hpp> // necessary, since not a standard import of Kokkos (experimental)

namespace ParticleBinning {
    template <typename BunchType>
    class AdaptBins {
    public:
        using value_type             = typename BunchType::Layout_t::value_type;
        using particle_position_type = typename BunchType::particle_position_type;
        using position_view_type     = typename particle_position_type::view_type;
        using size_type              = typename BunchType::size_type;

        // New bin_type representing the bin index data type
        using bin_index_type         = typename BunchType::bin_index_type; // has to be initialized over there...
        using bin_type               = typename ippl::ParticleAttrib<bin_index_type>;
        using bin_view_type          = typename bin_type::view_type; 

        using bin_histo_type         = typename Kokkos::View<size_type*>; // ippl::Vector<size_type*, 1>;
        using bin_host_histo_type    = typename Kokkos::View<size_type*, Kokkos::HostSpace>; // ippl::Vector<size_type*, 1>;

        // Constructor to initialize with a shared pointer to the particle container
        AdaptBins(std::shared_ptr<BunchType> bunch, bin_index_type maxBins = 10)
            : bunch_m(bunch)
            , maxBins_m(maxBins) {

            currentBins_m = maxBins; // TODO for now...
        }

        bin_index_type getCurrentBinCount() const { return currentBins_m; }

        void setCurrentBinCount(bin_index_type nBins) {
            currentBins_m = (nBins > maxBins_m) ? maxBins_m : nBins; 
            binWidth_m    = (xMax_m - xMin_m) / currentBins_m; // assuming particles did not change!
        }

        KOKKOS_INLINE_FUNCTION
        static bin_index_type getBin(value_type x, value_type xMin, value_type xMax, value_type binWidthInv, bin_index_type numBins); 

        void initLimits();
        void initializeHistogram(bool setToZero = false);
        void assignBinsToParticles();
        void initLocalHisto();
        bin_host_histo_type getGlobalHistogram();
        
        void doFullRebin(bin_index_type nBins) {
            setCurrentBinCount(nBins);
            assignBinsToParticles();
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

    private:
        std::shared_ptr<BunchType> bunch_m;  // Shared pointer to the particle container
        bin_index_type maxBins_m;  // Maximum number of bins
        bin_index_type currentBins_m; // Current number of bins
        value_type xMin_m;
        value_type xMax_m;
        value_type binWidth_m;
        bin_histo_type localBinHisto_m;
    };
}

#include "AdaptBins.hpp"

#endif  // ADAPT_BINS_H




