#ifndef OPAL_SAMPLE_GAUSSIAN_SEQUENCE_H
#define OPAL_SAMPLE_GAUSSIAN_SEQUENCE_H

#include "Sample/SamplingMethod.h"

#ifdef WITH_UNIT_TESTS
#include <gtest/gtest_prod.h>
#endif

class SampleGaussianSequence : public SamplingMethod
{

public:

    SampleGaussianSequence(double lower, double upper, int nSample)
        : n_m(0)
    {
        double mean = 0.5 * (lower + upper);
        double sigma = (upper - lower) / 10; // +- 5 sigma
        double factor = sigma / sqrt(2);
        double dx = 2.0 / nSample;
        double oldY = -2.5;
        for (long i = 1; i < nSample; ++ i) {
            double x = -1.0 + i * dx;
            double y = erfinv(x);
            chain_m.push_back(mean + factor * (y + oldY));
            oldY = y;
        }

        chain_m.push_back(mean + factor * (2.5 + oldY));

    }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = getNext();
    }

    double getNext() {
        double sample = chain_m[n_m];
        incrementCounter();

        return sample;
    }

private:
#ifdef WITH_UNIT_TESTS
    FRIEND_TEST(GaussianSampleTest, ChainTest);
#endif
    std::vector<double> chain_m;
    unsigned int n_m;

    void incrementCounter() {
        ++ n_m;
        if (n_m >= chain_m.size())
            n_m = 0;
    }

#define erfinv_a3 -0.140543331
#define erfinv_a2 0.914624893
#define erfinv_a1 -1.645349621
#define erfinv_a0 0.886226899

#define erfinv_b4 0.012229801
#define erfinv_b3 -0.329097515
#define erfinv_b2 1.442710462
#define erfinv_b1 -2.118377725
#define erfinv_b0 1

#define erfinv_c3 1.641345311
#define erfinv_c2 3.429567803
#define erfinv_c1 -1.62490649
#define erfinv_c0 -1.970840454

#define erfinv_d2 1.637067800
#define erfinv_d1 3.543889200
#define erfinv_d0 1

double erfinv (double x)
{
    double x2, r, y;
    int  sign_x;

    if (x < -1 || x > 1)
        return NAN;

    if (x == 0)
        return 0;

    if (x > 0)
        sign_x = 1;
    else {
        sign_x = -1;
        x = -x;
    }

    if (x <= 0.7) {

        x2 = x * x;
        r =
            x * (((erfinv_a3 * x2 + erfinv_a2) * x2 + erfinv_a1) * x2 + erfinv_a0);
        r /= (((erfinv_b4 * x2 + erfinv_b3) * x2 + erfinv_b2) * x2 +
              erfinv_b1) * x2 + erfinv_b0;
    }
    else {
        y = sqrt (-log ((1 - x) / 2));
        r = (((erfinv_c3 * y + erfinv_c2) * y + erfinv_c1) * y + erfinv_c0);
        r /= ((erfinv_d2 * y + erfinv_d1) * y + erfinv_d0);
    }

    r = r * sign_x;
    x = x * sign_x;

    r -= (erf (r) - x) / (2 / sqrt (M_PI) * exp (-r * r));
    r -= (erf (r) - x) / (2 / sqrt (M_PI) * exp (-r * r));

    return r;
}

#undef erfinv_a3
#undef erfinv_a2
#undef erfinv_a1
#undef erfinv_a0

#undef erfinv_b4
#undef erfinv_b3
#undef erfinv_b2
#undef erfinv_b1
#undef erfinv_b0

#undef erfinv_c3
#undef erfinv_c2
#undef erfinv_c1
#undef erfinv_c0

#undef erfinv_d2
#undef erfinv_d1
#undef erfinv_d0
};

#endif