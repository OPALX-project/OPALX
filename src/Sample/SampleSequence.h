#ifndef OPAL_SAMPLE_SEQUENCE_H
#define OPAL_SAMPLE_SEQUENCE_H

#include "Sample/SamplingMethod.h"

template <typename T>
class SampleSequence : public SamplingMethod
{
    // provides a sequence of equidistant sampling points. It
    // can't be garanteed that the sampling is equidistant if
    // an integer type is chosen and the difference between
    // the upper and lower limit isn't divisible by the number
    // of sampling points.

public:

    SampleSequence(T lower, T upper, size_t modulo, int nSample)
        : lowerLimit_m(lower)
        , stepSize_m( (upper - lower) / double(nSample - 1) )
        , sampleNr_m(0)
        , numSamples_m(nSample)
        , volumeLowerDimensions_m(modulo)
        , individualCounter_m(0)
    { }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = static_cast<T>(lowerLimit_m + stepSize_m * sampleNr_m);
        incrementCounter();
    }

private:
    void incrementCounter() {
        ++ individualCounter_m;
        if (individualCounter_m % volumeLowerDimensions_m == 0)
            ++ sampleNr_m;

        sampleNr_m = sampleNr_m % numSamples_m;
    }

    T lowerLimit_m;
    double stepSize_m;
    unsigned int sampleNr_m;
    unsigned int numSamples_m; // size of this "dimension"
    size_t volumeLowerDimensions_m; // the "volume" of the sampling space of the lower "dimensions"
    size_t individualCounter_m; // counts how many "individuals" have been created
};

#endif