/** slice.h
   slice class definition

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   07-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: slice.h 110 2007-05-15 08:52:36Z bakker $
*/


#ifndef _SLICE_DEF
#define _SLICE_DEF

#include <stdio.h>

/*
 #include "element.h"
*/

/** Index for parameters of each slice: see class below 
   beta normalized velocity (total) [-]
   z    slice position [m]
   x    beam size x (rms) [m]   
   y    beam size y (rms) [m]
   px   beam divergence x [rad]
   py   beam divergence y [rad]
   X0   position centroid x [m]
   Y0   position centroid y [m]
   pX0  angular deflection centriod x
   pY0  angular deflection centroid y

   N.B. INDEX MUST MATCH DEFINITION IN RK ROUTINE IN bunch.C !!!!!!!
*/

#define SLNPAR 10

#define SLI_z     0
#define SLI_beta  1
#define SLI_x     2
#define SLI_px    3
#define SLI_y     4
#define SLI_py    5
#define SLI_x0    6
#define SLI_px0   7
#define SLI_y0    8
#define SLI_py0   9

class EnvelopeSlice {
    friend class EnvelopeBunch;
    /**
      ElementList 
     *eList;    // list of active elements
     */

    double *p;        // parameters: array[SLNPAR]
    double *p_old;    // parameters before last integration step: array [SLNPAR]
    double *p_scr;    // parameters at screen position (= NULL if not defined)

    int p_local;   /* p-array (see above)  defined locally ?
                      (may be set to external to speed up MPI) */

    int valid;     /* slice is valid ?
                      (flag cleared on orrurance of beta <= 0) */
    int was_valid; // backup copy of valid

public:
    EnvelopeSlice(FILE *f = NULL);
    ~EnvelopeSlice();

    /// sets the memory buffer to a new location
    /// pointer to buffer
    void setBuffer(double *);

    /// read slice parameters from file
    int read(FILE *f);
    /// write slice parameters to file
    void write(FILE *f = stdout);
    /// write header of a slice file
    void writeHeader(FILE *f = stdout);

    /// get gamma value
    double gamma();

    /**
      void findActive(           // mark all elements that act between
      double,                     // z and z+c*dt (time step input)
      Element *);                 // first element in system
    */

    /// backup p to p_old
    void backup();
    /// restore the previous state
    void restore();

    // mark parameters at screen position, returns 1 if set
    // position of the screen
    // system time
    int mark(double, double);
    
    /// clear marking
    void clear_mark();

    /// check validity of slice, return true on change
    int  check();
    /// request validity of slice
    int  is_valid();
};

typedef EnvelopeSlice *EnvelopeSliceP;

#endif
