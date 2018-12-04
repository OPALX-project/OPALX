#ifndef OPAL_SAMPLING_METHOD_H
#define OPAL_SAMPLING_METHOD_H


#include "Sample/SampleIndividual.h"

#include <boost/smart_ptr.hpp>

class SamplingMethod
{

public:
    virtual void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) = 0;
    
    /*!
     * Allocate memory for sampling. Not every sampling method
     * requires that.
     * 
     * This function is used to reduce memory since only the
     * sampler ranks need these sampling methods.
     * 
     * @param n number of samples
     */
    virtual void allocate(std::size_t n, std::size_t seed) {
        /* Some sampling methods require a container.
         * In order to reduce memory only samplers should allocate
         * the memory
         */
    }
};

#endif