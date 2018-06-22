#ifndef RNGSTREAM_H
#define RNGSTREAM_H

#include <random>

class RNGStream
{
public:
    static RNGStream* getInstance();
    static RNGStream* getInstance(unsigned int seed);
    static void deleteInstance(RNGStream* & generator);

    static void setGlobalSeed(unsigned int seed);

    std::mt19937_64 & getGenerator();

    template <class DISTR>
    typename DISTR::result_type getNext(DISTR & RNGDist) {
        return RNGDist(RNGenerator_m);
    }
private:
    RNGStream():
        RNGenerator_m(globalSeed_sm),
        isGlobal_m(true)
    { }

    RNGStream(unsigned int seed):
        RNGenerator_m(seed),
        isGlobal_m(false)
    { }

    ~RNGStream()
    { }

    static RNGStream *globalInstance_sm;
    static unsigned int globalSeed_sm;
    static unsigned int numGlobalInstances_sm;
    std::mt19937_64 RNGenerator_m;
    bool isGlobal_m;

};
#endif