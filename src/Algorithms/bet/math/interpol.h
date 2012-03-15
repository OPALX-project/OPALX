/* interpol.h
   interpolation function definition

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   09-03-06      Created                                         Rene Bakker
*/


#ifndef _INTERPOL_DEF
#define _INTERPOL_DEF

/* spline()
   Given arrays x[0..n-1] and y[0..n-1] containing a tabulated function,
   i.e., y[i] = f(x[i]), with x[1] < x[2] < .. . < x[N-1], this routine
   returns an array y2[0..n-1] that contains the second derivatives of the
   interpolating function at the tabulated points x[i].
   Derrivate at bounderies set to zero !
*/
void spline
(
    double x[],
    double y[],
    int n,
    double y2[]);

/* splint()
   Given the arrays xa[0..n-1] and ya[0..n-1], which tabulate a function
   (with the xa[i] in order), and given the array y2a[0..n-1], which
   is the output from spline above, and given a value of x, this routine
   returns a cubic-spline interpolated value y.
*/
void splint
(
    double xa[],
    double ya[],
    double y2a[],
    int n,
    double x,
    double *y);

/* lsplint()
   Given the arrays xa[0..n-1] and ya[0..n-1], which tabulate a function
   (with the xa[i] in order), and given the array y2a[0..n-1], which
   is the output from spline above, and given a value of x, this routine
   returns a cubic-spline interpolated value y, which is limited between
   the y-values of the adjacent points.
*/
void lsplint
(
    double xa[],
    double ya[],
    double y2a[],
    int n,
    double x,
    double *y);

#endif
