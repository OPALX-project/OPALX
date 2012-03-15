/* linfit.h
   linear fitting routine

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   09-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: root.h 29 2007-04-14 17:03:18Z l_bakker $
*/


#ifndef _LINFIT_DEF
#define _LINFIT_DEF


/* linfit() Given a set of data points x[0..ndata-1],y[0..ndata-1] with
   individual standard deviations sig[0..ndata-1], fit them to a
   straight line y = a + bx by minimizing chi2. Returned are a,b and
   their respective probable uncertainties siga and sigb, the
   chi-square chi2, and the goodness-of-fit probability q (that the
   fit would have chi2 this large or larger). If mwt=0 on input, then
   the standard deviations are assumed to be unavailable: q is
   returned as 1.0 and the normalization of chi2 is to unit standard
   deviation on all points.  
*/

void linfit(double x[], double y[], int ndata, 
	    double sig[], int mwt, double *a, double *b, double *siga,
	    double *sigb, double *chi2, double *q);


/* linfit() 
   Given a set of data points x[0..ndata-1],y[0..ndata-1], fit them to
   a straight line y = a + bx by minimizing chi2. Returned are a,b and
   their respective probable uncertainties siga and sigb, and the
   chi-square chi2.
*/
void linfit(double x[], double y[], int ndata, 
	    double *a, double *b, double *siga,
	    double *sigb, double *chi2);

#endif
