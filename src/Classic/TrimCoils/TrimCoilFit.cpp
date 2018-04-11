#include "TrimCoils/TrimCoilFit.h"

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
}
