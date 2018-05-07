#ifndef OPAL_UNIFORM_H
#define OPAL_UNIFORM_H

#include "Sample/SamplingOperator.h"

#include <random>
#include <type_traits>

template <typename T>
class Uniform : public SamplingOperator
{

public:
    typedef std::conditional<
                std::is_integral<T>::value,
                std::uniform_int_distribution<>,
                std::uniform_real_distribution<>
            > dist_t;
    
    
    Uniform(int seed, T lower, T upper)
        : eng_m(seed)
        , dist_m(lower, upper)
        
    { }
    
    
    void create(boost::shared_ptr<Individual>& ind, int i) {
        ind->genes[i] = dist_m(eng_m);
    }
    
private:
    std::mt19937_64 eng_m;
    
    dist_t dist_m;
};

#endif
