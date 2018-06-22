#include "TrimCoil.h"

TrimCoil::TrimCoil(double bmax,
                   double rmin,
                   double rmax)
{
    // convert to m
    const double mm2m = 0.001;
    rmin_m = rmin * mm2m;
    rmax_m = rmax * mm2m;
    // convert to kG
    bmax_m = bmax * 10.0;
}

void TrimCoil::applyField(const double r, const double z, double *br, double *bz)
{
    doApplyField(r,z,br,bz);
}