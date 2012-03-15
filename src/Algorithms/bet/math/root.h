/* root.h
   root finding functions

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   09-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: root.h 29 2007-04-14 17:03:18Z l_bakker $
*/


#ifndef _ROOT_DEF
#define _ROOT_DEF

/* findRoot()
   Originally rtsafe()
   Using a combination of Newton-Raphson and bisection, find the root of a
   function bracketed between x1 and x2. The root, returned as the function
   value rtsafe, will be refined until its accuracy is known within ±xacc.
   funcd is a user-supplied routine that returns both the function value
   and the first derivative of the function. */
double findRoot
(
    void (*funcd)(double, double *, double *), // function to find root
    double, double,                            // boundaries
    double);                                   // accuracy

#endif
