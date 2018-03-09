#ifndef CLASSIC_STENCIL_HH
#define CLASSIC_STENCIL_HH

#include "Filters/Filter.h"

class IlyaPogorelovFilter: public Filter {
public:
    IlyaPogorelovFilter() { ;}
    void apply(std::vector<double> &histogram);
    void calc_derivative(std::vector<double> &histogram, const double &h);
};

#endif // CLASSIC_STENCIL_HH
