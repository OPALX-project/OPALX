#ifndef OPAL_UNIFORM_H
#define OPAL_UNIFORM_H

#include "Sample/SamplingMethod.h"
#include "Sample/RNGStream.h"

#include <type_traits>

template <typename T>
class Uniform : public SamplingMethod
{

public:
    typedef typename std::conditional<
                        std::is_integral<T>::value,
                        std::uniform_int_distribution<T>,
                        std::uniform_real_distribution<T>
                     >::type dist_t;

    Uniform(T lower, T upper)
        : dist_m(lower, upper)
        , RNGInstance_m(nullptr)
    {
        
    }

    Uniform(T lower, T upper, int seed)
        : dist_m(lower, upper)
        , RNGInstance_m(nullptr)

    {
        
    }

    ~Uniform() {
        if ( RNGInstance_m )
            RNGStream::deleteInstance(RNGInstance_m);
    }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = RNGInstance_m->getNext(dist_m);
    }
    
    void allocate(std::size_t n, std::size_t seed) {
        RNGInstance_m = RNGStream::getInstance(seed);
    }

private:
    dist_t dist_m;
    
    RNGStream *RNGInstance_m;
};

#endif