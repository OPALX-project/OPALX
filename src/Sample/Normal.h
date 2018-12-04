#ifndef OPAL_NORMAL_RANDOM_SAMPLING_H
#define OPAL_NORMAL_RANDOM_SAMPLING_H

#include "Sample/SamplingMethod.h"
#include "Sample/RNGStream.h"

#include <type_traits>

class Normal : public SamplingMethod
{

public:
    typedef std::normal_distribution<double> dist_t;


    Normal(double lower, double upper)
        : dist_m(0.5 * (lower + upper), (upper - lower) / 10)
        , RNGInstance_m(nullptr)

    {}

    Normal(double lower, double upper, double seed)
        : dist_m(0.5 * (lower + upper), (upper - lower) / 10)
        , RNGInstance_m(nullptr)

    {}

    ~Normal() {
        if ( RNGInstance_m)
            RNGStream::deleteInstance(RNGInstance_m);
    }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = RNGInstance_m->getNext(dist_m);
    }
    
    void allocate(std::size_t n, std::size_t seed) {
        RNGInstance_m = RNGStream::getInstance(seed);
    }

private:
    dist_t dist_m;

    RNGStream *RNGInstance_m;
};

#endif