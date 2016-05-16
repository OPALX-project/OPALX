#include "Algorithms/PolynomialTimeDependence.h"

Inform &PolynomialTimeDependence::print(Inform &os) {
  Inform::FmtFlags_t ff = os.flags();
  os << std::scientific;
  for (unsigned int i=0; i<this->coeffs.size(); i++)
    os << "c" << i << "= " << this->coeffs[i] << " ";
  os << endl;
  os.flags(ff);
  return os;
}
