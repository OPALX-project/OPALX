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
        , numSamples_m(nSample)
        , volumeLowerDimensions_m(modulo)
    { }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        
        unsigned int id = ind->id;
        
        int bin = int(id / volumeLowerDimensions_m) % numSamples_m;
        
        ind->genes[i] = static_cast<T>(lowerLimit_m + stepSize_m * bin);
    }
    
private:
    T lowerLimit_m;
    double stepSize_m;
    unsigned int numSamples_m; // size of this "dimension"
    size_t volumeLowerDimensions_m; // the "volume" of the sampling space of the lower "dimensions"
};

#endif