#ifndef OPAL_SAMPLE_GAUSSIAN_SEQUENCE_H
#define OPAL_SAMPLE_GAUSSIAN_SEQUENCE_H

#include "Sample/SamplingMethod.h"
#include "Utilities/Util.h"

#ifdef WITH_UNIT_TESTS
#include <gtest/gtest_prod.h>
#endif

class SampleGaussianSequence : public SamplingMethod
{
    // provides a sequence of sampling points that have a Gaussian distribution
    // with
    //      mean = 0.5 * (upper + lower)
    //      sigma = (upper - lower) / 10
    // This can be achieved if the integral of the Gaussian between the sampling
    // points are all equal. The sampling points are therefore computed using
    // the inverse error function at equally distributed arguments between
    // -1 and 1.

public:

    SampleGaussianSequence(double lower, double upper, size_t modulo, int nSample)
        : sampleNr_m(0)
        , numSamples_m(nSample)
        , volumeLowerDimensions_m(modulo)
        , individualCounter_m(0)
    {
        double mean = 0.5 * (lower + upper);
        double sigma = (upper - lower) / 10; // +- 5 sigma
        double factor = sigma / sqrt(2);
        double dx = 2.0 / nSample;
        for (long i = 0; i < nSample; ++ i) {
            double x = -1.0 + (i + 0.5) * dx;
            double y = Util::erfinv(x);
            sampleChain_m.push_back(mean + factor * y);
        }
    }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = getNext();
    }

    double getNext() {
        double sample = sampleChain_m[sampleNr_m];
        incrementCounter();

        return sample;
    }

private:
#ifdef WITH_UNIT_TESTS
    FRIEND_TEST(GaussianSampleTest, ChainTest);
#endif
    std::vector<double> sampleChain_m;
    unsigned int sampleNr_m;
    unsigned int numSamples_m; // size of this "dimension"
    size_t volumeLowerDimensions_m; // the "volume" of the sampling space of the lower "dimensions"
    size_t individualCounter_m; // counts how many "individuals" have been created

    void incrementCounter() {
        ++ individualCounter_m;
        if (individualCounter_m % volumeLowerDimensions_m == 0)
            ++ sampleNr_m;

        sampleNr_m = sampleNr_m % numSamples_m;
    }
};

#endif