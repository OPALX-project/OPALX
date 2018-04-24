#ifndef TRIM_COILMIRRORED_H
#define TRIM_COILMIRRORED_H

#include "TrimCoils/TrimCoil.h"

/// TrimCoilMirrored class
/// Shape mirrored from TC-15 shape
/// http://accelconf.web.cern.ch/AccelConf/ipac2017/papers/thpab077.pdf

class TrimCoilMirrored : public TrimCoil{

public:
    TrimCoilMirrored(double bmax,
                     double rmin,
                     double rmax,
                     double slope);

    virtual ~TrimCoilMirrored() { };

private:
    /// Maximum B field (kG)
    double bmax_m;
    /// Minimum radius (m)
    double rmin_m;
    /// Maximum radius (m)
    double rmax_m;
    /// Slope in (1 / mm)
    double bslope_m;

    /// @copydoc TrimCoil::doApplyField
    virtual void doApplyField(const double r, const double z, double *br, double *bz);

    TrimCoilMirrored() = delete;
};

#endif //TRIM_COILMIRRORED_H
