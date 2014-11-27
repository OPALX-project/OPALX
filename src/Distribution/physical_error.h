#ifndef PHYSICAL_ERROR_H
#define PHYSICAL_ERROR_H

#include "error.h"

/// This class writes unphysical errors to the terminal and stops execution.
/*!
 * It inherits from the base class Error.
 */
class PhysicalError : public Error
{
public:
  /// No stationary distribution can be found due to instability of system
  static std::string match;
  
  /// Unphysical value
  static std::string negative;
  
  /// Undefined value
  static std::string undefined;
  
  /// Unstable solution
  static std::string unstable;
  
  /// Imaginary solution
  static std::string imag;
  
};

std::string PhysicalError::match = "MATCHING NOT POSSIBLE";
std::string PhysicalError::negative = "NEGATIVE VALUE";
std::string PhysicalError::undefined = "UNDEFINED PHYSICAL VALUE";
std::string PhysicalError::unstable = "NOT STABLE SOLUTION";
std::string PhysicalError::imag = "IMAGINARY VALUE";

#endif