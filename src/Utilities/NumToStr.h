#ifndef OPAL_NumToStr_HH
#define OPAL_NumToStr_HH 1

// ------------------------------------------------------------------------
// $RCSfile: NumToStr.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:12 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include <sstream>
#include <string>

/// Convert number to string.
//  Return the string corresponding to the numeric argument.

template <typename T>
std::string NumToStr(T t) {
    std::ostringstream oss;
    oss << t;
    return oss.str();
}

#endif // OPAL_NumToStr_HH
