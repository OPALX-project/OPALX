#ifndef OPALX_STENCIL_HH
#define OPALX_STENCIL_HH

#include "Filters/Filter.h"

class StencilFilter : public Filter {
public:
    StencilFilter() { ; }
    void apply(std::vector<double>& histogram);
    void calc_derivative(std::vector<double>& histogram, const double& h);
};

#endif  // OPALX_STENCIL_HH