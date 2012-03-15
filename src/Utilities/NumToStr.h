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
// Functions: string NumToStr(short)
//            string NumToStr(unsigned short)
//            string NumToStr(int)
//            string NumToStr(unsigned int)
//            string NumToStr(long)
//            string NumToStr(unsigned long)
//            string NumToStr(float)
//            string NumToStr(double)
//            string NumToStr(long double)
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:12 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include <sstream>

/// Convert number to string.
//  Return the string corresponding to the numeric argument.

std::string NumToStr(short s);

std::string NumToStr(unsigned short s);

std::string NumToStr(int i);

std::string NumToStr(unsigned int i);

std::string NumToStr(long l);

std::string NumToStr(unsigned long l);

std::string NumToStr(float f);

std::string NumToStr(double d);

std::string NumToStr(long double d);


#endif // OPAL_NumToStr_HH
