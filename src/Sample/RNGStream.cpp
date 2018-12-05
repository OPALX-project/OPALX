#include "Sample/RNGStream.h"

#include <iostream>

RNGStream * RNGStream::globalInstance_sm = NULL;
unsigned int RNGStream::globalSeed_sm = 42;
unsigned int RNGStream::numGlobalInstances_sm = 0;

RNGStream* RNGStream::getInstance() {
    if (globalInstance_sm == NULL)
        globalInstance_sm = new RNGStream();

    ++ numGlobalInstances_sm;
    return globalInstance_sm;
}

RNGStream* RNGStream::getInstance(unsigned int seed) {
    return new RNGStream(seed);
}

void RNGStream::deleteInstance(RNGStream* & generator) {
    if (generator->isGlobal_m) {
        -- numGlobalInstances_sm;

        if (numGlobalInstances_sm == 0) {
            delete generator;
        }
    } else {
        delete generator;
    }

    generator = NULL;
    return;
}

void RNGStream::setGlobalSeed(unsigned int seed) {
    globalSeed_sm = seed;

    if (globalInstance_sm != NULL)
        globalInstance_sm->RNGenerator_m.seed(seed);
}

unsigned int RNGStream::getGlobalSeed() {
    return globalSeed_sm;
}

std::mt19937_64 & RNGStream::getGenerator() {
    return RNGenerator_m;
}