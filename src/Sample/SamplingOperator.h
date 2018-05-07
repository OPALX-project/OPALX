#ifndef OPAL_SAMPLING_OPERATOR_H
#define OPAL_SAMPLING_OPERATOR_H


#include "Sample/Individual.h"

class SamplingOperator
{
    
public:
    virtual void create(Individual& ind, int i) = 0;
};

#endif
