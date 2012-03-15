/* integrate.h
   integration routines

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   09-03-06      Created                                         Rene Bakker
*/


#ifndef _INTEGRATE_DEF
#define _INTEGRATE_DEF

/* qromb()
   Returns the integral of the function func from a to b
   with error tolerance eps. 
   Integration is performed by Romberg's method. The order
   is set in in the program source code (default 5).
*/
double qromb(
  double (*)(double), // function 
  double, double,     // boundaries
  double = 1.0e-4);   // error (eps)

#endif
