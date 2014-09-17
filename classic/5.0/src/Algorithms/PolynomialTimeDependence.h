#ifndef _CLASSIC_SRC_ALGORITHMS_POLYNOMIALTIMEDEPENDENCE_H_
#define _CLASSIC_SRC_ALGORITHMS_POLYNOMIALTIMEDEPENDENCE_H_

#include <vector>
#include <iostream>

#include "Algorithms/AbstractTimeDependence.h"

class PolynomialTimeDependence : public AbstractTimeDependence {
  public:
    PolynomialTimeDependence(std::vector<double> ptd) : coeffs(ptd) {}
    PolynomialTimeDependence() {}
    ~PolynomialTimeDependence() {}
    inline double getValue(double time);
    PolynomialTimeDependence* clone() {
      std::vector<double> temp(coeffs);
      PolynomialTimeDependence* d = new PolynomialTimeDependence(temp);
      return d;
    }

  private:
    std::vector<double> coeffs;
};

double PolynomialTimeDependence::getValue(double time) {
    double x = 0.;
    double t_power = 1.;
    for (std::size_t i = 0; i < coeffs.size() ; ++i) {
        x += coeffs[i]*t_power;
        t_power *= time;
    }
    return x;
}

#endif

