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

    {
        RNGInstance_m = RNGStream::getInstance();
    }

    Normal(double lower, double upper, double seed)
        : dist_m(0.5 * (lower + upper), (upper - lower) / 10)

    {
        RNGInstance_m = RNGStream::getInstance(seed);
    }

    ~Normal() {
        RNGStream::deleteInstance(RNGInstance_m);
    }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = RNGInstance_m->getNext(dist_m);
    }

private:
    RNGStream *RNGInstance_m;

    dist_t dist_m;
};

#endif