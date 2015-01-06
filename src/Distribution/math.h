/**
 * @file math.h
 * Mathematical functions that aren't part of the standard of C++
 *
 * @author Matthias Frey
 * @version 1.0
 */

#ifndef MATH_H
#define MATH_H

#include <cmath>

/*!
 *  \addtogroup matt
 *  @{
 */

/// @brief Defines additional mathematical functions
namespace matt {
    /// Computes the sign of the input argument
    template<typename T>
        T sign(const T val) {
            return (std::signbit(val)) ? T(-1) : T(1);
        }
}

/*! @} End of Doxygen Groups*/

#endif
