#ifndef OPAL_SAMPLING_METHOD_H
#define OPAL_SAMPLING_METHOD_H


#include "Sample/SampleIndividual.h"

#include <boost/smart_ptr.hpp>

class SamplingMethod
{
    
public:
    virtual void create(boost::shared_ptr<SampleIndividual>& ind, int i) = 0;
};

#endif
