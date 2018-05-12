#ifndef OPAL_NORMAL_RANDOM_SAMPLING_H
#define OPAL_NORMAL_RANDOM_SAMPLING_H

#include "Sample/SamplingMethod.h"

#include <random>
#include <type_traits>

class Normal : public SamplingMethod
{

public:
    typedef std::normal_distribution<double> dist_t;


    Normal(double lower, double upper, double seed)
        : eng_m(seed)
        , dist_m(0.5 * (lower + upper), (upper - lower) / 10)

    { }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = dist_m(eng_m);
    }

private:
    std::mt19937_64 eng_m;

    dist_t dist_m;
};

#endif