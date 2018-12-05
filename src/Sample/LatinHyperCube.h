#ifndef OPAL_LATIN_HYPERCUBE_H
#define OPAL_LATIN_HYPERCUBE_H

#include "Sample/SamplingMethod.h"
#include "Sample/RNGStream.h"

#include <algorithm>
#include <deque>

class LatinHyperCube : public SamplingMethod
{

public:
    typedef typename std::uniform_real_distribution<double> dist_t;

    LatinHyperCube(double lower, double upper)
        : binsize_m(0.0)
        , upper_m(upper)
        , lower_m(lower)
        , dist_m(0.0, 1.0)
    {
        RNGInstance_m = RNGStream::getInstance();
    }
    
    LatinHyperCube(double lower, double upper, int seed)
        : binsize_m(0.0)
        , upper_m(upper)
        , lower_m(lower)
        , dist_m(0.0, 1.0)
    {
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
    
    void allocate(std::size_t n) {
        
        binsize_m = ( upper_m - lower_m ) / double(n);
        
        this->fillBins_m(n);
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
        
        return  binsize_m * (val + bin) + lower_m;
    }
    
    void fillBins_m(std::size_t n) {
        bin_m.resize(n);
        std::iota(bin_m.begin(), bin_m.end(), 0);
        
        std::random_device rd;
        std::mt19937_64 eng(rd());
        std::shuffle(bin_m.begin(), bin_m.end(), eng);
    }

private:
    RNGStream *RNGInstance_m;
    
    std::deque<std::size_t> bin_m;
    double binsize_m;
    
    double upper_m;
    double lower_m;
    
    dist_t dist_m;
};

#endif
