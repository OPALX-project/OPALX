#ifndef TRIM_COILPHASEFIT_H
#define TRIM_COILPHASEFIT_H

#include "TrimCoils/TrimCoil.h"

#include <vector>

/// TrimCoilPhaseFit class
/// General rational function fit of the phase shift

class TrimCoilPhaseFit : public TrimCoil {

public:
    TrimCoilPhaseFit(double bmax,
                     double rmin,
                     double rmax,
                     const std::vector<double>& coefnum,
                     const std::vector<double>& coefdenom);

    virtual ~TrimCoilPhaseFit() { };

private:
    TrimCoilPhaseFit() = delete;

    std::vector<double> coefnum_m;
    std::vector<double> coefdenom_m;

    /// @copydoc TrimCoil::doApplyField
    virtual void doApplyField(const double r, const double z, double *br, double *bz);
};

#endif //TRIM_COILPHASEFIT_H
