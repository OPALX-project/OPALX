#ifndef STEPPERS_H
#define STEPPERS_H

#include "RK4.h"
#include "LF2.h"
#include "MTS.h"


namespace stepper {
    
    enum INTEGRATOR {
        UNDEFINED   = -1,
        RK4         = 0,
        LF2         = 1,
        MTS         = 2
    };
};

#endif
