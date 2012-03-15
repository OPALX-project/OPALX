// ------------------------------------------------------------------------
// $RCSfile: NumToStr.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Functions: std::string NumToStr(short)
//            std::string NumToStr(unsigned short)
//            std::string NumToStr(int)
//            std::string NumToStr(unsigned int)
//            std::string NumToStr(long)
//            std::string NumToStr(unsigned long)
//            std::string NumToStr(float)
//            std::string NumToStr(double)
//            std::string NumToStr(long double)
//   Convert number to string.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:12 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Utilities/NumToStr.h"


std::string NumToStr(short s) {
    std::ostringstream oss;
    oss << s;
    return oss.str();
}

std::string NumToStr(unsigned short s) {
    std::ostringstream oss;
    oss << s;
    return oss.str();
}

std::string NumToStr(int i) {
    std::ostringstream oss;
    oss << i;
    return oss.str();
}

std::string NumToStr(unsigned int i) {
    std::ostringstream oss;
    oss << i;
    return oss.str();
}

std::string NumToStr(long l) {
    std::ostringstream oss;
    oss << l;
    return oss.str();
}

std::string NumToStr(unsigned long l) {
    std::ostringstream oss;
    oss << l;
    return oss.str();
}

std::string NumToStr(float f) {
    std::ostringstream oss;
    oss << f;
    return oss.str();
}

std::string NumToStr(double d) {
    std::ostringstream oss;
    oss << d;
    return oss.str();
}

std::string NumToStr(long double d) {
    std::ostringstream oss;
    oss << d;
    return oss.str();
}

