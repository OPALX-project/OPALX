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

/*
template<typename SizeType, typename IndexType>
struct HistogramReducer {
    using bin_histo_type = typename Kokkos::View<SizeType*>; // long int
    using bin_idx_type   = typename Kokkos::View<IndexType*>; // e.g. uint8_t

    bin_histo_type global_histogram;   
    //bin_idx_type binIndex;  // Bin index for each particle
    IndexType numBins;

    // Default constructor (needed for Kokkos) (I don't know why, but it works...)
    HistogramReducer() : global_histogram() { init(); }

    // Constructor to initialize the views
    HistogramReducer(bin_histo_type global_histogram_, const SizeType& numBins_) 
        : global_histogram(global_histogram_)
        , numBins(numBins_) {} //  

    // Initialize the local histogram in each thread
    KOKKOS_INLINE_FUNCTION
    void init() {
        global_histogram = bin_histo_type("global_histogram", 10);
        Kokkos::deep_copy(global_histogram, 0);  // Initialize local histograms to 0
    }

    // Join the reduction result
    KOKKOS_INLINE_FUNCTION
    void join(bin_histo_type& local_histogram, const bin_histo_type& other_histogram) const {
        for (IndexType i = 0; i < numBins; ++i) {
            local_histogram(i) += other_histogram(i);
        }
    }

    // Final reduction step that adds the local histogram into the global histogram
    KOKKOS_INLINE_FUNCTION
    void final(const bin_histo_type& local_histogram) const {
        for (IndexType i = 0; i < numBins; ++i) {
            Kokkos::atomic_add(&global_histogram(i), local_histogram(i));
        }
    }

    // Calls join method to combine results
    KOKKOS_INLINE_FUNCTION
    void operator+=(const HistogramReducer& other) {
        join(global_histogram, other.global_histogram);
    }

    // This operator is used during parallel_reduce to update the histogram
    //KOKKOS_INLINE_FUNCTION
    //void operator()(const SizeType i, HistogramReducer& update) const {
        // Increment the count in the histogram for that bin

    //    global_histogram(binIndex(i))++;
    //}
};*/
/*
// Define the custom scalar type for the reduction
template<typename SizeType, typename IndexType>
struct HistogramStruct {
    using bin_histo_type = typename Kokkos::View<SizeType*>;

    bin_histo_type bins;  // Histogram with bin counts

    // Default constructor (required for Kokkos)
    KOKKOS_INLINE_FUNCTION
    HistogramStruct();

    // Constructor to allocate histogram with a given number of bins
    HistogramStruct(IndexType num_bins) : bins("histogram", num_bins) {}

    // Merge two histograms (for parallel reduction)
    KOKKOS_INLINE_FUNCTION
    void operator+=(const HistogramStruct& other) {
        for (IndexType i = 0; i < bins.extent(0); ++i) {
            bins(i) += other.bins(i);
        }
    }
};

// Custom reducer to manage the histogram reductions
template<typename SizeType, typename IndexType>
struct HistogramReducer {
    using value_type = HistogramStruct<SizeType, IndexType>;

    IndexType numBins;  // Number of bins (runtime value)

    // Constructor: set the number of bins (provided at runtime)
    HistogramReducer(IndexType num_bins) : numBins(num_bins) {}

    // Initializes the local histogram to zero
    KOKKOS_INLINE_FUNCTION
    void init(value_type& histogram) const {
        histogram = value_type(numBins);  // Allocate bins and initialize to 0
        Kokkos::deep_copy(histogram.bins, 0);  // Initialize to 0
    }

    // Reduces (joins) two histograms by adding their bin counts
    KOKKOS_INLINE_FUNCTION
    void join(value_type& dest, const value_type& src) const {
        dest += src;
    }

    // Finalize the reduction: accumulate local histogram to the global histogram
    KOKKOS_INLINE_FUNCTION
    void final(value_type& local, value_type& global) const {
        for (IndexType i = 0; i < numBins; ++i) {
            Kokkos::atomic_add(&global.bins(i), local.bins(i));
        }
    }
};*/

/*template<typename SizeType, typename IndexType>
struct ViewReducer {
    //using value_type = Kokkos::View<SizeType*>;
    //SizeType the_array[N];
    SizeType* bins;

    static IndexType numBins; // changed later
    static void setNumBins(IndexType count) { numBins = count; }

    KOKKOS_INLINE_FUNCTION   // Called only inside "static ViewReducer<SizeType, IndexType> sum()" 
    ViewReducer() { 
        init();
    }

    KOKKOS_INLINE_FUNCTION
    ~ViewReducer() { delete[] bins; }

    KOKKOS_INLINE_FUNCTION
    void init() {
        bins = new SizeType[numBins];
        for (IndexType i = 0; i < numBins; i++) { bins[i] = 0; }
    }
    
    KOKKOS_INLINE_FUNCTION   // add operator
    ViewReducer& operator+= (const ViewReducer& src) {
        for (IndexType i = 0; i < numBins; i++) { bins[i] += src.bins[i]; }
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    ViewReducer& operator=(const ViewReducer& rhs) {
        if (this != &rhs) { 
            if (bins) delete[] bins;
            numBins = rhs.numBins;
            bins = new SizeType[numBins];
            for (IndexType i = 0; i < numBins; i++) { bins[i] = rhs.bins[i]; }
        }
        return *this;
    }
};*/

/*template<typename SizeType, typename IndexType>
struct ViewReducer {
    using value_type = Kokkos::View<SizeType*>;
    value_type bins;  // Kokkos View to store bin values

    static IndexType numBins;  // Static member to store the number of bins

    static void setNumBins(IndexType count) { numBins = count; 
        std::cout << "Hey setNumBins" << std::endl; }

    KOKKOS_INLINE_FUNCTION   // Default constructor
    ViewReducer() {  // Initialize Kokkos View with numBins size
    
        std::cout << "Hey Empty" << std::endl;
        bins = value_type("local_histogram", numBins);
        Kokkos::deep_copy(bins, 0);  // Initialize bins to 0
    }

    KOKKOS_INLINE_FUNCTION   // Copy Constructor // Maybe delete this???????
    ViewReducer(const ViewReducer& rhs) { //}: bins("local_histogram", numBins) {
        
        std::cout << "Hey" << std::endl;
        bins = value_type("local_histogram", numBins);
        Kokkos::deep_copy(bins, rhs.bins);  // Deep copy the content
    }

    KOKKOS_INLINE_FUNCTION   // Destructor (Kokkos View handles memory cleanup automatically)
    ~ViewReducer() { std::cout << "Hey Destructor" << std::endl;}

    KOKKOS_INLINE_FUNCTION   // Add operator
    ViewReducer& operator+= (const ViewReducer& src) {
        std::cout << "Hey+=" << std::endl;
        for (IndexType i = 0; i < numBins; i++) {
            bins(i) += src.bins(i);  // Accumulate the bin values
        }
        return *this;
    }

    KOKKOS_INLINE_FUNCTION   // Assignment operator
    ViewReducer& operator=(const ViewReducer& rhs) {
        std::cout << "Hey=" << std::endl;
        if (this != &rhs) {
            Kokkos::deep_copy(bins, rhs.bins);  // Deep copy the bin values
        }
        return *this;
    }
};*/

/*template<typename SizeType, typename IndexType>
struct ViewReducer {
    using value_type = Kokkos::View<SizeType*>;
    value_type bins;  // Kokkos View to store bin values

    // Constructor to initialize the view with the given number of bins
    KOKKOS_INLINE_FUNCTION
    ViewReducer(IndexType numBins) : bins("local_histogram", numBins) {
        Kokkos::deep_copy(bins, 0);  // Initialize bins to 0
    }

    // Copy Constructor
    KOKKOS_INLINE_FUNCTION
    ViewReducer(const ViewReducer& rhs) : bins("local_histogram", rhs.bins.extent(0)) {
        Kokkos::deep_copy(bins, rhs.bins);  // Deep copy the content
    }

    // Destructor
    KOKKOS_INLINE_FUNCTION
    ~ViewReducer() {}

    // Add operator
    KOKKOS_INLINE_FUNCTION
    ViewReducer& operator+= (const ViewReducer& src) {
        //auto numBins = bins.extent(0);  // Get the number of bins from the current view
        for (IndexType i = 0; i < bins.extent(0); i++) {
            bins(i) += src.bins(i);  // Accumulate the bin values
        }
        return *this;
    }

    // Assignment operator
    KOKKOS_INLINE_FUNCTION
    ViewReducer& operator=(const ViewReducer& rhs) {
        if (this != &rhs) {
            Kokkos::deep_copy(bins, rhs.bins);  // Deep copy the bin values
        }
        return *this;
    }
};*/



/*struct ViewReducer {
    using value_type = typename Kokkos::View<SizeType*>;

    value_type bins = value_type("local_histogram", N);

    KOKKOS_INLINE_FUNCTION   // Default constructor - Initialize to 0's
    ViewReducer() {
        for (IndexType i = 0; i < N; ++i) {
            bins(i) = 0.0;
        }
    }

    KOKKOS_INLINE_FUNCTION   // Copy Constructor
    ViewReducer(const ViewReducer& rhs) { 
        bins = rhs.bins;
    }

    KOKKOS_INLINE_FUNCTION   // add operator
    ViewReducer& operator+=(const ViewReducer& src) {
        for (IndexType i = 0; i < N; ++i) {
            bins(i) += src.bins(i);
        }
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    ViewReducer& operator=(const ViewReducer& rhs) {
        if (this != &rhs) {  // Prevent self-assignment
            bins = rhs.bins;
            //Kokkos::deep_copy(bins, rhs.bins);  // Deep copy the bins
            //Kokkos::fence();
        }
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const IndexType& idx) {
        ++bins(idx);
    }
};*/
//typedef array_type<int,4> ValueType;  // used to simplify code below

// Specialize the reduction identity for ViewReducer
/*namespace Kokkos {
    template<typename SizeType, typename IndexType>
    struct reduction_identity<ViewReducer<SizeType, IndexType>> {
        KOKKOS_FORCEINLINE_FUNCTION
        static ViewReducer<SizeType, IndexType> sum() {
            std::cout << "Hey reduction_identity" << std::endl;
            return ViewReducer<SizeType, IndexType>(); // Pass N for initialization
        }
    };
}*/
/*namespace Kokkos {
    template<typename SizeType, typename IndexType>
    struct reduction_identity<ViewReducer<SizeType, IndexType>> {
        static IndexType numBins;  // Static variable to hold the number of bins

        // Function to set the number of bins before parallel_reduce
        static void setNumBins(IndexType count) {
            std::cout << "Hey setNumBins" << count << std::endl;
            numBins = count;
        }

        // Sum initialization with dynamic numBins
        KOKKOS_FORCEINLINE_FUNCTION
        static ViewReducer<SizeType, IndexType> sum() {
            std::cout << "Hey reduction_identity" << std::endl;
            return ViewReducer<SizeType, IndexType>(numBins);  // Initialize with numBins
        }
    };

    // Define the static member
    template<typename SizeType, typename IndexType>
    IndexType reduction_identity<ViewReducer<SizeType, IndexType>>::numBins = 0;
}*/

/*
template<typename SizeType, typename IndexType>
struct ViewReducer {
    SizeType* bins;
    IndexType numBins;

    KOKKOS_INLINE_FUNCTION
    ViewReducer(IndexType binsCount) 
        : numBins(binsCount) {
        bins = new SizeType[binsCount]();
    }

    KOKKOS_INLINE_FUNCTION
    ViewReducer(const ViewReducer& rhs)
        : numBins(rhs.numBins) {
        bins = new SizeType[numBins];
        for (IndexType i = 0; i < numBins; ++i) {
            bins[i] = rhs.bins[i];
        }
    }

    KOKKOS_INLINE_FUNCTION
    ~ViewReducer() {
        delete[] bins;
    }

    KOKKOS_INLINE_FUNCTION
    void operator+=(const ViewReducer& src) {
        for (IndexType i = 0; i < numBins; ++i) {
            bins[i] += src.bins[i];
        }
    }

    KOKKOS_INLINE_FUNCTION
    void join(ViewReducer& dest, const ViewReducer& src) const {
        for (IndexType i = 0; i < src.numBins; ++i) {
            Kokkos::atomic_add(&dest.bins[i], src.bins[i]);
        }
    }
};

// Define reduction_identity for the custom reducer
namespace Kokkos {
    template<typename SizeType, typename IndexType>
    struct reduction_identity<ViewReducer<SizeType, IndexType>> {
        KOKKOS_INLINE_FUNCTION static ViewReducer<SizeType, IndexType> sum() {
            return ViewReducer<SizeType, IndexType>(0);
        }
    };
}*/









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

    // Get the current number of bins
    bin_index_type getCurrentBinCount() const {
        return currentBins_m;
    }

    void setCurrentBinCount(bin_index_type nBins) {
        // Set the number of current bins to the new value!
        currentBins_m = (nBins > maxBins_m) ? maxBins_m : nBins; 
        binWidth_m    = (xMax_m - xMin_m) / currentBins_m; // assuming particles did not change!
    }

    void initializeHistogram() {
        // Reinitialize the histogram view with the new size (numBins)
        const bin_index_type numBins = getCurrentBinCount();
        localBinHisto_m = bin_histo_type("binHisto_m", numBins);
        bin_histo_type localBinHisto = localBinHisto_m; // avoid implicit "this" capture

        // Optionally, initialize the histogram to zero
        Kokkos::parallel_for("initHistogram", numBins, KOKKOS_LAMBDA(const bin_index_type i) {
            localBinHisto(i) = 0;
        });
        
        //Kokkos::deep_copy(localBinHisto_m, 0);
        //Kokkos::fence();  // Ensure initialization is done before proceeding
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
        //bin_histo_type localBinHisto = localBinHisto_m; // use copy of histo... view to avoid implicit "this" capture in Kokkos:parallel_for

        // Declare the variables locally before the Kokkos::parallel_for (to avoid implicit this capture in Kokkos lambda)
        value_type xMin = xMin_m, xMax = xMax_m, binWidthInv = 1.0/binWidth_m;
        bin_index_type numBins = currentBins_m;
        // Alternatively explicit capture: [xMin = xMin_m, xMax = xMax_m, binWidth = binWidth_m, numBins = currentBins_m, localData = localData, binIndex = binIndex]

        // Capture variables explicitly...
        static IpplTimings::TimerRef assignParticleBins = IpplTimings::getTimer("assignParticleBins");
        IpplTimings::startTimer(assignParticleBins);
        Kokkos::parallel_for("assignParticleBins", bunch_m->getLocalNum(), KOKKOS_LAMBDA(const size_type i) {
                // Access the z-axis position of the i-th particle
                value_type v = localData(i)[2];  
                
                // Assign the bin index to the particle (directly on device)
                bin_index_type bin = getBin(v, xMin, xMax, binWidthInv, numBins);
                binIndex(i)        = bin;
                //localBinHisto(bin) += 1;
                //Kokkos::atomic_fetch_add(&localBinHisto(bin), static_cast<bin_index_type>(1)); // since "Kokkos::atomic_add" might not support unsigned long, only int?!
                //Kokkos::atomic_add(&localBinHisto(bin), 1);
                // Kokkos::atomic_increment(&localBinHisto(bin)); // Not needed with parallel_reduce
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
        initializeHistogram(); // Initialize all bin counts to 0 --> not necessary when using parallel_reduce!
        msg << "Histogram initialized to 0" << endl;
        //for (size_type i = 0; i < 10; ++i) {
            //msg << "Hey" << endl;
        //    msg << localBinHisto_m(i) << " " << endl;
        //}

        bin_view_type binIndex       = bunch_m->bin.getView();
        bin_histo_type localBinHisto = localBinHisto_m;
        bin_index_type binCount      = getCurrentBinCount();
        //Kokkos::reduction_identity<reducer_type>::setBinCount(binCount);

        //reducer_type::setNumBins(binCount);
        //Kokkos::reduction_identity<reducer_type>::setNumBins(binCount);
        //reducer_type to_reduce(binCount); 
        
        //reducer_type resultReducer(binCount);
        
        //resultReducer.init();  
        //ViewReducerFunctor<size_type, bin_index_type> functor(resultReducer, binCount);
        // msg << "Initializing histogram of size " << to_reduce.numBins << endl;

        //reducer_type reducer = reducer_type();
        //msg << "HEY1" << endl;
        /*for (int i = 0; i < 10; ++i) {
            msg << to_reduce.bins(i) << endl;
        }*/
        //static IpplTimings::TimerRef initLocalHisto = IpplTimings::getTimer("initLocalHisto");
        //IpplTimings::startTimer(initLocalHisto);
        msg << "Starting reducer." << endl;

        //Kokkos::parallel_reduce("custom_reduce", Kokkos::RangePolicy<>(0, bunch_m->getLocalNum()), functor, resultReducer);

        /*Kokkos::parallel_reduce("initLocalHist", bunch_m->getLocalNum(),
            KOKKOS_LAMBDA(const size_type i, reducer_type& update) {
                bin_index_type ndx = binIndex(i);
                Kokkos::atomic_increment(&update.bins[ndx]);
                //update.bins(ndx)++;
            }, Kokkos::Sum<reducer_type>(resultReducer));*/

        Kokkos::Experimental::ScatterView<size_type*> scatter_view("scatter_view", binCount);

        // Usage with parallel_for
        Kokkos::parallel_for("initLocalHist", bunch_m->getLocalNum(), 
            KOKKOS_LAMBDA(const size_type i) {
                auto scatter = scatter_view.access();
                bin_index_type ndx = binIndex(i);
                scatter(ndx)++;
            });

        auto hist_view = scatter_view.subview();
        Kokkos::parallel_for("finalize_histogram", binCount, KOKKOS_LAMBDA(const size_type i) {
            localBinHisto(i) = hist_view(i);
        });
        Kokkos::fence();



        /*Kokkos::parallel_reduce("initLocalHist", bunch_m->getLocalNum(), 
            KOKKOS_LAMBDA(const size_type i, reducer_type& update) {
                std::cout << i << std::endl;
                bin_index_type ndx = binIndex(i);  // Get the bin index for the particle
                update.bins(ndx)++;  // Update the appropriate bin count
            }, Kokkos::Sum<reducer_type>(to_reduce)  // Pass the custom reducer
        );*/
        //Kokkos::fence();
        msg << "Reducer ran without error." << endl;

        //Kokkos::deep_copy(localBinHisto_m, tr.the_array);
        //msg << "HEY2" << endl;
        /*for (bin_index_type i = 0; i < binCount; ++i) {
            msg << localBinHisto(i) << endl;
            //localBinHisto_m(i) = resultReducer.bins(i); // [i];
        }*/
        //IpplTimings::stopTimer(initLocalHisto);
    }

    void doFullRebin(bin_index_type nBins) {
        // Sets number of bins, assigns bin numbers and initializes histogram
        setCurrentBinCount(nBins);
        assignBinsToParticles();
    }

    bin_histo_type getGlobalHistogram() {
        Inform msg("GetGlobalHistogram");
        // Get the current number of bins
        bin_index_type numBins = getCurrentBinCount(); // number of local bins = number of global bins!
        //msg << "1" << endl;
        // Create a view to hold the global histogram on all ranks
        bin_histo_type globalBinHisto("globalBinHisto", numBins);
        //msg << "2" << endl;
        // Loop through each bin and reduce its value across all nodes
        for (bin_index_type i = 0; i < numBins; ++i) {
            //msg << "3 - " << i << endl;
            size_type localValue = localBinHisto_m(i);
            size_type globalValue = 0;

            // Perform the reduction to sum the values across all nodes
            ippl::Comm->allreduce(localValue, globalValue, 1, std::plus<size_type>());  // Reduce to root node
            //msg << "4 - " << i << endl;
            // Store the global value only on rank 0
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

        os << "Bins = " << getCurrentBinCount() << " hBin = " << binWidth_m << endl;
        os << "Bin #;Val" << endl;

        // Get the global histogram (reduced across all nodes)
        bin_histo_type globalBinHisto = getGlobalHistogram();

        // Only rank 0 prints the global histogram
        for (bin_index_type i = 0; i < getCurrentBinCount(); ++i) {
            os << i << ";" << globalBinHisto(i) << endl;
        }

        size_type total = 0;
        // Compute the sum of all elements in the histogram
        Kokkos::parallel_reduce("sum_histogram", globalBinHisto.extent(0), 
            KOKKOS_LAMBDA(const size_type i, size_type& localSum) {
                localSum += globalBinHisto(i);
            }, total);

        os << "Total = " << total << endl;
        os << "-----------------------------------------" << endl;
    }
};



#endif 


