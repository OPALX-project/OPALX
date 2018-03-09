/* functions.h
   implementation of special functions

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   09-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: functions.h 29 2007-04-14 17:03:18Z l_bakker $
*/


#ifndef _FUNCTIONS_DEF
#define _FUNCTIONS_DEF

/* besselj()
   Interger bessel function of the first kind
   n = 0, 1, ......
*/
double bessj(int n, double x);

/* complementary error function */
double errfc(double x);

#endif
