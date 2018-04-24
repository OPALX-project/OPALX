#include "TrimCoil.h"

TrimCoil::TrimCoil(){}

void TrimCoil::applyField(const double r, const double z, double *br, double *bz)
{
    doApplyField(r,z,br,bz);
}