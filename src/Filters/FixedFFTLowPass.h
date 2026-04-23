#ifndef OPALX_FIXED_LOW_PASS_FFT_HH
#define OPALX_FIXED_LOW_PASS_FFT_HH

#include "Filters/Filter.h"

class FixedFFTLowPassFilter: public Filter {
public:
    FixedFFTLowPassFilter(const int &N);
    void apply(std::vector<double> &LineDensity);
    void calc_derivative(std::vector<double> &histogram, const double &h);

private:
    int number_frequencies_m;
};

#endif // OPALX_FIXED_LOW_PASS_FFT_HH
