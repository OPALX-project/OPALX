#ifndef OPAL_SAMPLING_OPERATOR_H
#define OPAL_SAMPLING_OPERATOR_H


#include "Sample/Individual.h"

#include <boost/smart_ptr.hpp>

class SamplingOperator
{
    
public:
    virtual void create(boost::shared_ptr<Individual>& ind, int i) = 0;
};

#endif
