#ifndef TRIM_COILBFIT_H
#define TRIM_COILBFIT_H

#include "TrimCoils/TrimCoil.h"

#include <vector>

/// TrimCoilBFit class
/// General rational function fit
/// https://gitlab.psi.ch/OPAL/src/issues/157

class TrimCoilBFit : public TrimCoil {

public:
    TrimCoilBFit(double bmax,
                 double rmin,
                 double rmax,
                 const std::vector<double>& coefnum,
                 const std::vector<double>& coefdenom);

    virtual ~TrimCoilBFit() { };

private:
    TrimCoilBFit() = delete;

    std::vector<double> coefnum_m;
    std::vector<double> coefdenom_m;

    /// @copydoc TrimCoil::doApplyField
    virtual void doApplyField(const double r, const double z, double *br, double *bz);
};

#endif //TRIM_COILBFIT_H
