#include "TrimCoils/TrimCoilFit.h"

#include <cmath>

TrimCoilFit::TrimCoilFit(double bmax,
                         double rmin,
                         double rmax,
                         const std::vector<double>& coefnum,
                         const std::vector<double>& coefdenom):
    TrimCoil(),
    bmax_m(bmax),
    coefnum_m(coefnum),
    coefdenom_m(coefdenom)
{
    // convert to m
    const double mm2m = 0.001;
    rmin_m = rmin * mm2m;
    rmax_m = rmax * mm2m;
}

void TrimCoilFit::doApplyField(const double r, const double z, double *br, double *bz)
{
    // FIXME

    if (std::abs(bmax_m) < 1e-20) return;
    // check range
    if (r < rmin_m || r > rmax_m) return;

    double num = 0.0;
    for (std::size_t i = 0; i < coefnum_m.size(); ++i) {
        num += std::pow(coefnum_m[i], i);
    }
    double denom = 0.0;
    for (std::size_t i = 0; i < coefdenom_m.size(); ++i) {
        denom += std::pow(coefdenom_m[i], i);
    }

    double dr = num / denom;

    *bz -= 0.0;
    *br -= dr*z;
}
