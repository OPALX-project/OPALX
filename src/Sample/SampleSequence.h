#ifndef OPAL_SAMPLE_SEQUENCE_H
#define OPAL_SAMPLE_SEQUENCE_H

#include "Sample/SamplingMethod.h"

template <typename T>
class SampleSequence : public SamplingMethod
{

public:

    SampleSequence(T lower, T upper, size_t modulo, int nSample)
        : lower_m(lower)
        , step_m( (upper - lower) / double(nSample - 1) )
        , n_m(0)
        , size_m(nSample)
        , counter_m(0)
        , mod_m(modulo)

    { }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = static_cast<T>(lower_m + step_m * n_m);
        incrementCounter();
    }

private:
    void incrementCounter() {
        ++ counter_m;
        if (counter_m % mod_m == 0)
            ++ n_m;
        if (n_m >= size_m)
            n_m = 0;
    }

    T lower_m;
    double step_m;
    unsigned int n_m;
    int size_m;
    size_t counter_m;
    size_t mod_m;
};

#endif