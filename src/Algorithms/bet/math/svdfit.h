/* svdfit.h
   SVD fitting routines

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   09-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: svdfit.h 29 2007-04-14 17:03:18Z l_bakker $
*/


#ifndef _SVDFIT_DEF
#define _SVDFIT_DEF

/* svdfit 
   Given a set of data points x[0..ndata-1],y[0..ndata-1], with individual standard 
   deviations sig[1..ndata], use X^2 minimization to determine the coefficients 
   a[0..ma-1] of the fitting function y = Sum a[i]*afunc[i](x). 
   Here we solve the fitting equations using singular value decomposition of the 
   ndata by ma matrix, as in paragraph 2.6. 

   The program returns values for the  ma fit parameters a, and X^2, chisq. 
   The user supplies a routine funcs(x,afunc,ma) that returns the ma basis functions 
   evaluated at x = x in the array afunc[0..ma-1].
*/

void svdfit
(
 double *, // x[]
 double *, // y[]
 double *, // sig[]
 int,      // ndata
 double *, // a[]
 int,      // ma
 double *, // error[]
 double *, // chisq
 void (*)(double, double *, int));

void svdfit // like previous with sig[i] = 0.0, i [0..ndata-1]
(
 double *, // x[]
 double *, // y[]
 int,      // ndata
 double *, // a[]
 int,      // ma
 double *, // error[]
 double *, // chisq
 void (*)(double, double *, int));

void svdfitP // Polynomial fit of order ma
(
 double *,  // x[]
 double *,  // y[]
 int,       // ndata
 double *,  // a[]
 int,       // ma
 double *,  // error[]
 double *); // chisq

#endif
