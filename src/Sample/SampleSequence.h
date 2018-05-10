#ifndef OPAL_SAMPLE_SEQUENCE_H
#define OPAL_SAMPLE_SEQUENCE_H

#include "Sample/SamplingMethod.h"

class SampleSequence : public SamplingMethod
{

public:
    
    SampleSequence(double lower, double upper, int nSample)
        : lower_m(lower)
        , step_m( (upper - lower) / double(nSample - 1) )
        , n_m(0)
        
    { }
    
    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = lower_m + step_m * n_m;
        ++n_m;
    }
    
private:
    double lower_m;
    double step_m;
    int n_m;
};

#endif
