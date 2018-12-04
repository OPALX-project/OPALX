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
        , RNGInstance_m(RNGStream::getInstance())
        , seed_m(RNGStream::getGlobalSeed())
    {}
    
    LatinHyperCube(double lower, double upper, int seed)
        : binsize_m(0.0)
        , upper_m(upper)
        , lower_m(lower)
        , dist_m(0.0, 1.0)
        , RNGInstance_m(nullptr)
        , seed_m(seed)
    {}

    ~LatinHyperCube() {
        if ( RNGInstance_m )
            RNGStream::deleteInstance(RNGInstance_m);
    }

    void create(boost::shared_ptr<SampleIndividual>& ind, std::size_t i) {
        /* values are created within [0, 1], thus, they need to be mapped
         * the domain [lower, upper]
         */
        ind->genes[i] = map2domain_m(RNGInstance_m->getNext(dist_m));
    }
    
    void allocate(const CmdArguments_t& args, const Comm::Bundle_t& comm) {
        int id = comm.island_id;
        
        if ( !RNGInstance_m )
            RNGInstance_m = RNGStream::getInstance(seed_m + id);
        
        int nSamples = args->getArg<int>("nsamples", true);
        int nMasters = args->getArg<int>("num-masters", true);
        
        int nLocSamples = nSamples / nMasters;
        int rest = nSamples - nMasters * nLocSamples;
        
        if ( id < rest )
            nLocSamples++;
        
        int startBin = 0;
        
        if ( rest == 0 )
            startBin = nLocSamples * id;
        else {
            if ( id < rest ) {
                startBin = nLocSamples * id;
            } else {
                startBin = (nLocSamples + 1) * rest + (id - rest) * nLocSamples;
            }
        }
        
        binsize_m = ( upper_m - lower_m ) / double(nSamples);
        
        this->fillBins_m(nSamples, startBin, seed_m + id);
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
    
    void fillBins_m(std::size_t n, int startBin, std::size_t seed) {
        bin_m.resize(n);
        std::iota(bin_m.begin(), bin_m.end(), startBin);
        
        std::mt19937_64 eng(seed);
        std::shuffle(bin_m.begin(), bin_m.end(), eng);
    }

private:
    std::deque<std::size_t> bin_m;
    double binsize_m;
    
    double upper_m;
    double lower_m;
    
    dist_t dist_m;
    
    RNGStream *RNGInstance_m;
    
    std::size_t seed_m;
};

#endif
