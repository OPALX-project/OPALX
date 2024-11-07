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

namespace ParticleBinning {

    /**
     * @class AdaptBins
     * @brief A class that bins particles in energy bins and allows for runtime rebinning.
     * 
     * This class provides functionality to group particles into energy bins, initialize and reinitialize
     * bin structures, and update bin contents using data structures from IPPL and Kokkos. It is optimized
     * for usage with MPI and CUDA backends.
     *
     * @tparam BunchType The type of particle bunch (container) used in the binning process.
     */
    template <typename BunchType>
    class AdaptBins {
    public:
        using value_type             = typename BunchType::Layout_t::value_type;
        using particle_position_type = typename BunchType::particle_position_type;
        using position_view_type     = typename particle_position_type::view_type;
        using size_type              = typename BunchType::size_type;
        using bin_index_type         = typename BunchType::bin_index_type;
        using bin_type               = typename ippl::ParticleAttrib<bin_index_type>;
        using bin_view_type          = typename bin_type::view_type;
        using bin_histo_type         = Kokkos::View<size_type*>;
        using bin_host_histo_type    = Kokkos::View<size_type*, Kokkos::HostSpace>;

        /**
         * @brief Constructs an AdaptBins object with a specified maximum number of bins.
         * 
         * @param bunch A shared pointer to the particle container.
         * @param maxBins The maximum number of bins to initialize with (default is 10).
         */
        AdaptBins(std::shared_ptr<BunchType> bunch, bin_index_type maxBins = 10)
            : bunch_m(bunch)
            , maxBins_m(maxBins) {

            currentBins_m = maxBins; // TODO for now...
        }

        /**
         * @brief Gets the current number of bins.
         * 
         * @return The current bin count.
         */
        bin_index_type getCurrentBinCount() const { return currentBins_m; }

        /**
         * @brief Sets the current number of bins and adjusts the bin width.
         * 
         * @param nBins The new number of bins.
         */
        void setCurrentBinCount(bin_index_type nBins) {
            currentBins_m = (nBins > maxBins_m) ? maxBins_m : nBins; 
            binWidth_m    = (xMax_m - xMin_m) / currentBins_m; // assuming particles did not change!
        }

        /**
         * @brief Calculates the bin index for a given position value.
         * 
         * This static method calculates which bin a position value falls into based on the bin
         * boundaries and bin width.
         * 
         * @param x The "position" value.
         * @param xMin Minimum boundary for the bins.
         * @param xMax Maximum boundary for the bins.
         * @param binWidthInv Inverse of the bin width for efficiency.
         * @param numBins The total number of bins.
         * @return The calculated bin index for the position value.
         */
        KOKKOS_INLINE_FUNCTION
        static bin_index_type getBin(value_type x, value_type xMin, value_type xMax, value_type binWidthInv, bin_index_type numBins);

        /**
         * @brief Initializes the limits for binning based on the particle data.
         * 
         * This function calculates the minimum and maximum limits (xMin and xMax) from the
         * particle positions, which are then used to define bin boundaries.
         */
        void initLimits();

        /**
         * @brief Initializes the histogram view for binning and optionally sets it to zero.
         * 
         * @param setToZero If true, initializes the histogram view to zero. Default is false. The 0 initialization is not needed if it is overwritten anyways.
         */
        void initializeHistogram(bool setToZero = false);

        /**
         * @brief Assigns each particle in the bunch to a bin based on its position.
         * 
         * This function iterates over all particles in the bunch, calculates their bin
         * index, and updates the bin structure accordingly.
         */
        void assignBinsToParticles();

        /**
         * @brief Initializes a local histogram view for particle binning.
         * 
         * This function prepares a local histogram, enabling binning of particles within
         * local data, which can then be reduced into a global histogram.
         */
        void initLocalHisto();

        /**
         * @brief Retrieves the global histogram across all processes.
         * 
         * This function reduces the local histograms across all MPI processes into
         * a single global histogram view.
         * 
         * @return A view of the global histogram in host space (used for debugging output).
         */
        bin_host_histo_type getGlobalHistogram();

        /**
         * @brief Performs a full rebinning of particles with a specified number of bins.
         * 
         * @param nBins The new number of bins to use for rebinning.
         */
        void doFullRebin(bin_index_type nBins) {
            setCurrentBinCount(nBins);
            assignBinsToParticles();
        }


        /**
         * @brief Prints the current global histogram to the output stream.
         * 
         * This function outputs the global histogram data (bin counts) to the standard output.
         * Note: Only works correctly for rank 0 in an MPI environment.
         */
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

        /**
         * @brief Outputs debug information related to Kokkos and MPI configurations.
         * 
         * This function prints information about the number of threads (in OpenMP) or GPUs
         * (in CUDA) available on the current MPI rank, along with other debug information.
         */
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

        //template<bin_index_type N = 1>
        //void performLocalHistoReduction(bin_index_type binCount);


    private:
        std::shared_ptr<BunchType> bunch_m;    ///< Shared pointer to the particle container.
        bin_index_type maxBins_m;              ///< Maximum number of bins.
        bin_index_type currentBins_m;          ///< Current number of bins in use.
        value_type xMin_m;                     ///< Minimum boundary for bins.
        value_type xMax_m;                     ///< Maximum boundary for bins.
        value_type binWidth_m;                 ///< Width of each bin.
        bin_histo_type localBinHisto_m;        ///< Local histogram view for bin counts.

    };

}

#include "ParallelReduceTools.h"
#include "AdaptBins.hpp"

#endif  // ADAPT_BINS_H


