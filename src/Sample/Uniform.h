#ifndef OPAL_UNIFORM_H
#define OPAL_UNIFORM_H

#include "Sample/SamplingMethod.h"

#include <random>
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
    
    
    Uniform(T lower, T upper, int seed)
        : eng_m(seed)
        , dist_m(lower, upper)
        
    { }
    
    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = dist_m(eng_m);
    }
    
private:
    std::mt19937_64 eng_m;
    
    dist_t dist_m;
};

#endif
