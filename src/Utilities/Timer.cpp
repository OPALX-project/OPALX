// ------------------------------------------------------------------------
// $RCSfile: Timer.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Timer
//   Get Calendar date and time.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:48 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Utilities/Timer.h"
#include <ctime>

namespace OPALTimer {
    // Class Timer
    // ------------------------------------------------------------------------

    Timer::Timer() {
        ::time(&timer);
    }


    Timer::~Timer()
    {}


    string Timer::date() const {
        char buffer[12];
        strftime(buffer, 12, "%d/%m/%Y", localtime(&timer));
        return string(buffer, 10);
    }


    string Timer::time() const {
        char buffer[12];
        strftime(buffer, 12, "%H.%M.%S", localtime(&timer));
        return string(buffer, 8);
    }
}
