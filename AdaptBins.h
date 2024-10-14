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

#include "Ippl.h"

/*
// Define a functor for finding min and max values of the third component
struct MinMaxFunctor {
    Kokkos::View<double**> particles; // Assuming particles is a 2D view [N][3]

    MinMaxFunctor(Kokkos::View<double**> particles) : particles(particles) {}

    // The functor operator for the reduction
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, double& minValue, double& maxValue) const {
        // Update the min and max with the third component
        double component = particles(i, 2); // Accessing the third component (index 2)
        minValue = Kokkos::min(minValue, component);
        maxValue = Kokkos::max(maxValue, component);
    }
};*/

// The following custom struct is used to find the Nd minimum and maximum simultaneously in a Kokkos:view
/*template<typename T, int Dim>
struct MinMaxVectorReducer {
    Vector_t<T, Dim> min_val;
    Vector_t<T, Dim> max_val;

    // Constructor initializes the min and max values
    KOKKOS_INLINE_FUNCTION
    MinMaxReducer() {
        for (int d = 0; d < Dim; ++d) {
            min_val[d] = std::numeric_limits<T>::max();
            max_val[d] = std::numeric_limits<T>::lowest();
        }
    }

    // Update min/max values during the reduction operation
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, MinMaxReducer& update) const {
        for (int d = 0; d < Dim; ++d) {
            // Compare the current thread's min/max with the global update
            update.min_val[d] = std::min(update.min_val[d], min_val[d]);
            update.max_val[d] = std::max(update.max_val[d], min_val[d]);
        }
    }

    // Combine the results from different threads
    KOKKOS_INLINE_FUNCTION
    void join(const MinMaxReducer& update) {
        for (int d = 0; d < Dim; ++d) {
            min_val[d] = std::min(min_val[d], update.min_val[d]);
            max_val[d] = std::max(max_val[d], update.max_val[d]);
        }
    }

    // For OpenMP or reduction operations, we overload the += operator
    KOKKOS_INLINE_FUNCTION
    void operator+=(const MinMaxReducer& update) {
        join(update);
    }
};*/

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
    using bin_type               = typename ippl::ParticleAttrib<uint8_t>;
    using bin_view_type          = typename bin_type::view_type; 

    std::shared_ptr<BunchType> bunch_m;  // Shared pointer to the particle container
    size_type maxBins_m;  // Maximum number of bins
    size_type currentBins_m; // Current number of bins
    value_type xMin_m;
    value_type xMax_m;
    value_type binWidth_m;
    ippl::Vector<size_type*, 1> binHisto_m;

    // Constructor to initialize with a shared pointer to the particle container
    AdaptBins(std::shared_ptr<BunchType> bunch, size_type maxBins = 10)
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
            update.min_val = std::min(update.min_val, val);
            update.max_val = std::max(update.max_val, val);
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

    void initializeHistogram(size_type numBins) {
        // Reinitialize the histogram view with the new size (numBins)
        binHisto_m = ippl::Vector<size_type*, 1>(numBins);

        // Optionally, initialize the histogram to zero
        Kokkos::parallel_for("initHistogram", numBins, KOKKOS_LAMBDA(const int i) {
            binHisto_m(i) = 0;
        });
        Kokkos::fence();  // Ensure initialization is done before proceeding
    }

    KOKKOS_INLINE_FUNCTION
    static int getBin(value_type x, value_type xMin, value_type xMax, value_type binWidth, int numBins) {
        // Explanation: Don't access xMin, binWidth, ... through the members to avoid implicit
        // variable capture by Kokkos and potential copying overhead. Instead, pass them as an 
        // argument, s.t. Kokkos can capture them explicitly!
        // Make it static to avoid implicit capture of this inside Kokkos lambda! 

        // Ensure x is within bounds (clamp it between xMin and xMax --> this is only for bin assignment)
        x = (x < xMin) ? xMin : ((x > xMax) ? xMax : x);
        
        // Compute the bin index (ensuring it's within [0, numBins-1])
        int bin = (x - xMin) / binWidth;
        return (bin >= numBins) ? (numBins - 1) : bin;  // Clamp to the maximum bin
    }

    void assignBinsToParticles() {
        // Set the bin attribute for the given particle
        //bunch_m.bin(particleIndex) = bin;
        Inform msg("AdaptBins");

        position_view_type localData = bunch_m->R.getView();
        bin_view_type binIndex = bunch_m->bin.getView();  
        // Use Kokkos::parallel_for to iterate over the particles and assign them to bins
        /*Kokkos::parallel_for("assignParticleBins", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const int i) {
            value_type v = localData(i)[2];  // Assuming z-axis is the third dimension
            binIndex(i) = getBin(v, xMin, xMax, binWidth, numBins);
        });*/
        initializeHistogram(currentBins_m);

        // Declare the variables locally before the Kokkos::parallel_for (to avoid implicit this capture in Kokkos lambda)
        value_type xMin = xMin_m;
        value_type xMax = xMax_m;
        value_type binWidth = binWidth_m;
        int numBins = currentBins_m;
        // Alternatively explicit capture: [xMin = xMin_m, xMax = xMax_m, binWidth = binWidth_m, numBins = currentBins_m, localData = localData, binIndex = binIndex]

        // Capture variables explicitly...
        Kokkos::parallel_for("assignParticleBins", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const int i) {
            // Access the z-axis position of the i-th particle
            value_type v = localData(i)[2];  
            
            // Assign the bin index to the particle (directly on device)
            int bin     = getBin(v, xMin, xMax, binWidth, numBins);
            binIndex(i) = bin;
            Kokkos::atomic_increment(&binHisto_m(bin));
        });

        // Ensure all bins are assigned before further calculations!
        Kokkos::fence();   
        msg << "All bins assigned." << endl; 
    }

    // Get the current number of bins
    size_type getCurrentBinCount() const {
        return currentBins_m;
    }

    // Access bin attribute for a specific particle
    /*bin_type getBinForParticle(size_type particleIndex) const {
        return bunch_m.bin(particleIndex);
    }*/

    Inform& print(Inform& os) {
        os << "-----------------------------------------" << endl;
        os << "        Output Binning Structure         " << endl;

        os << "Bins = " << getCurrentBinCount() << " hBin = " << binWidth_m << " Particle vector length " << endl;
        os << "Bin #;Val" << endl;

        for (size_type i = 0; i < getCurrentBinCount(); ++i) {
            os << i << ";" << binHisto_m(i) << endl;
        }
    }
};



#endif 


