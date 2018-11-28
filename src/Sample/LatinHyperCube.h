#ifndef OPAL_LATIN_HYPERCUBE_H
#define OPAL_LATIN_HYPERCUBE_H

#include "Sample/SamplingMethod.h"
#include "Sample/RNGStream.h"

#include <algorithm>

class LatinHyperCube : public SamplingMethod
{

public:
    typedef typename std::uniform_real_distribution<double> dist_t;

    LatinHyperCube(double lower, double upper, std::size_t n)
        : binsize_m((upper - lower) / double(n))
        , dist_m(0.0, 1.0)
    {
        this->fillBins_m(n);
        
        RNGInstance_m = RNGStream::getInstance();
    }
    
    LatinHyperCube(double lower, double upper, int seed, std::size_t n)
        : binsize_m((upper - lower) / double(n))
        , dist_m(0.0, 1.0)
    {
        this->fillBins_m(n);
        
        RNGInstance_m = RNGStream::getInstance(seed);
    }

    ~LatinHyperCube() {
        RNGStream::deleteInstance(RNGInstance_m);
    }

    void create(boost::shared_ptr<SampleIndividual>& ind, std::size_t i) {
        /* values are created within [0, 1], thus, they need to be mapped
         * the domain [lower, upper]
         */
        ind->genes[i] = map2domain_m(RNGInstance_m->getNext(dist_m));
    }
    
private:
    double map2domain_m(double val) {
        /* y = mx + q
         * 
         * [0, 1] --> [a, b]
         * 
         * y = (b - a) * x + a
         * 
         * where a and b are the lower, respectively, upper
         * bound of the current bin.
         */
        
        std::size_t bin = bin_m.back();
        bin_m.pop_back();
        
        if ( bin_m.size() + 1e4 < bin_m.capacity() )
            bin_m.shrink_to_fit();
        
        double high = (bin + 1) * binsize_m;
        double low  = bin * binsize_m;
        
        return  (high - low) * val + low;
    }
    
    void fillBins_m(std::size_t n) {
        bin_m.reserve(n);
        std::iota(bin_m.begin(), bin_m.end(), 0);
        
        std::random_device rd;
        std::mt19937_64 eng(rd());
        std::shuffle(bin_m.begin(), bin_m.end(), eng);
    }

private:
    RNGStream *RNGInstance_m;
    
    std::vector<std::size_t> bin_m;
    double binsize_m;
    
    dist_t dist_m;
};

#endif
