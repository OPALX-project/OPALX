#ifndef PHYSICS_H
#define PHYSICS_H

#include <cmath>
#include <functional>

/*!
 *  \addtogroup physics
 *  @{
 */

/// Defines the typical constants in physics
namespace physics {
  /// Vacuum permittivity \f$ \left[\frac{As}{Vm} = \frac{C^{2}}{Nm^{2}} = \frac{F}{m}\right] \f$ (deu: elektrische Feldkonstante)
  long double eps0 = 8.854187817e-12;
  /// Vacuum permeability \f$ \left[\frac{N}{A^{2}} = \frac{Vs}{Am}\right] \f$ (deu: magnetische Feldkonstante)
  long double mu0 = 4.0e-7*M_PI;
  /// Elementary charge \f$ \left[C\right] \f$
  long double q0 = 1.0;	//1.60217733e-19;
  /// Speed of light \f$ \left[\frac{m}{s}\right] \f$
  long double c = 2.99792458e8;
  /// Proton mass \f$ \left[\frac{MeV}{c^{2}}\right] \f$
  long double E0=938.272;
  /// Cyclotron unit length \f$ \left[m\right] \f$
  /*!
   * This quantity is defined in the paper "Transverse-Longitudinal Coupling by Space Charge in Cyclotrons" 
   * of Dr. Christian Baumgarten (2012)
   * The lambda function takes the orbital frequency \f$ \omega_{o} \f$ (also defined in paper) as input argument.
   */
  std::function<double(double)> acon = [](double wo) { return c/wo; };
  /// Cyclotron unit \f$ \left[T\right] \f$ (Tesla)
  /*!
   * The lambda function takes the orbital frequency \f$ \omega_{o} \f$ as input argument.
   */
  std::function<double(double)> bcon = [](double wo) { return E0*1.0e7/(q0*c*acon(wo)); };
};

/*! @} End of Doxygen Groups*/

#endif