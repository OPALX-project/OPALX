#ifndef MATH_H
#define MATH_H

/*!
 *  \addtogroup matt
 *  @{
 */

/// Defines additional functions
namespace matt {
  /// Computes the sign of the input argument
  template<typename T>
  T sign(const T val) {
    return (val<0) ? T(-1) : T(1);
  }
}

/*! @} End of Doxygen Groups*/

#endif