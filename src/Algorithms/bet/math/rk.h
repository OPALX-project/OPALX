/* rk.h
   Runge-Kutta integration implementation

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   07-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: rk.h 29 2007-04-14 17:03:18Z l_bakker $
*/


#ifndef _RK_DEF
#define _RK_DEF

#include <stdio.h>

/* rk4()
   Given values for the variables y[0..n-1] use the fourth-order Runge-Kutta 
   method to advance the solution over an interval h and return the 
   incremented variables in y[0..n-1]. The user supplies the routine 
   derivs(x,y,dydx), which returns derivatives dydx at x.t
*/
void rk4(
  double y[],    // initial condition         
  int n, 	 // number of equations  
  double x, 	 // start	       
  double h, 	 // interval
  void (*derivs)(double, double [], double []));

/* odeint
   Runge-Kutta driver with adaptive stepsize control. 

   Integrate starting values ystart[0..nvar-1] 
   from x1 to x2 with accuracy eps, 

   storing intermediate results in global variables. 
   h1 should be set as a guessed first stepsize, 
   hmin as the minimum allowed stepsize (can be zero). 

   On output:
   - nok and nbad are the number of good and bad (but retried and fixed) 
     steps taken
   - ystart is replaced by values at the end of the integration interval. 
   - RETURNS 0 on success and 1 on failure
   
   derivs is the user-supplied routine for calculating the right-hand side 
   derivative.

   rkqs is the name of the stepper routine (see above).
*/
int odeint(
 double ystart[],  // initial condition  
 int    nvar,      // number of equations 
 double x1,        // start		 
 double x2,        // end		 
 double eps,       // accuracy		 
 double h1,        // guessed first step-size
 double hmin,      // minimum step-size	 
 int    *nok,      // number of good steps
 int    *nbad,     // number of bad steps
 void (*derivs)(double, double [], double []));

void rkActivateBuffer(int);             // > 0, 0 - deactivates
void rkPrintBuffer(FILE *f = stdout);

#endif
