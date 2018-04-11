#ifndef TRIM_COILFIT_H
#define TRIM_COILFIT_H

#include "TrimCoils/TrimCoil.h"

#include <vector>

/// TrimCoilFit class
/// General rational function fit
/// https://gitlab.psi.ch/OPAL/src/issues/157

class TrimCoilFit : public TrimCoil {

public:
    TrimCoilFit(double bmax,
                double rmin,
                double rmax,
                const std::vector<double>& coefnum,
                const std::vector<double>& coefdenom);

    virtual ~TrimCoilFit() { };

private:
    TrimCoilFit() = delete;
    /// Maximum B field (T)
    double bmax_m;
    /// Minimum radius (m)
    double rmin_m;
    /// Maximum radius (m)
    double rmax_m;

    std::vector<double> coefnum_m;
    std::vector<double> coefdenom_m;

    /// @copydoc TrimCoil::doApplyField
    virtual void doApplyField(const double r, const double z, double *br, double *bz);
};

#endif //TRIM_COILFIT_H
