#include "TrimCoils/TrimCoilPhaseFit.h"

#include <cmath>

TrimCoilPhaseFit::TrimCoilPhaseFit(double bmax,
                                   double rmin,
                                   double rmax,
                                   const std::vector<double>& coefnum,
                                   const std::vector<double>& coefdenom):
    TrimCoil(bmax, rmin, rmax),
    coefnum_m(coefnum),
    coefdenom_m(coefdenom)
{
    // normal polynom if no denominator coefficients (denominator = 1)
    if (coefdenom_m.empty())
      coefdenom_m.push_back(1.0);
}

void TrimCoilPhaseFit::doApplyField(const double r, const double z, double *br, double *bz)
{
    if (std::abs(bmax_m) < 1e-20) return;
    // check range
    if (r < rmin_m || r > rmax_m) return;

    double num     = 0.0; // numerator
    double d_num   = 0.0; // derivative of numerator
    double d2_num  = 0.0; // second derivative of numerator
    double powr    = 1.0; // power of r
    unsigned int order = coefnum_m.size();

    // add constant and first term
    num += coefnum_m[0];
    if (order > 1) {
        num   += coefnum_m[1] * r;
        d_num += coefnum_m[1];
    }
    for (std::size_t i = 2; i < coefnum_m.size(); ++i) {
        d2_num += coefnum_m[i] * powr * i * (i-1);
        powr   *= r; // r^(i-1)
        d_num  += coefnum_m[i] * powr * i;
        num    += coefnum_m[i] * powr * r;
    }

    double denom    = 0.0; // denominator
    double d_denom  = 0.0; // derivative of denominator
    double d2_denom = 0.0; // derivative of denominator
    powr            = 1.0; // power of r
    order           = coefdenom_m.size();

    // add constant
    denom += coefdenom_m[0];
    if (order > 1) {
        denom   += coefdenom_m[1] * r;
        d_denom += coefdenom_m[1];
    }
    for (std::size_t i = 2; i < coefdenom_m.size(); ++i) {
        d2_denom += coefdenom_m[i] * powr * i * (i-1);
        powr     *= r;
        d_denom  += coefdenom_m[i] * powr * i;
        denom    += coefdenom_m[i] * powr * r;
    }

    double phase = num / denom;

    // derivative of phase with quotient rule (B - field)
    double d_phase = (d_num * denom - d_denom * num) / (denom*denom);

    // second derivitive of phase (dB/dr) with quotient rule
    // (d2_num - 2*(num/denom)' * d_denom - (num/denom) * d2_denom) / denom
    double d2_phase = (d2_num - 2*d_phase*d_denom - phase * d2_denom) / denom;

    //std::cout << "r " << r << " dr " <<  dr << std::endl;

    *bz += - bmax_m * d_phase;
    *br += - bmax_m * d2_phase * z;
}