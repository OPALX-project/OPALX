/* savgol.h
   Savitzky-Golay Smoothing Filters

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   09-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: savgol.h 29 2007-04-14 17:03:18Z l_bakker $
*/


#ifndef _SAVGOL_DEF
#define _SAVGOL_DEF



/* savgol()
   Returns in c[1..np]
   consistent with the argument respns in
   routine convlv, a set of Savitzky-Golay filter coefficients. 
   nl is the number of leftward (past) data points used, while 
   nr is the number of rightward (future) data points, making 
      the total number of data points used nl+nr+1. 
   ld is the order of the derivative desired 
      (e.g., ld = 0 for smoothed function). 
   m  is the order of the smoothing polynomial, 
      also equal to the highest conserved moment; 
      usual values are m = 2 or m = 4.
*/
void savgol(
  double *,          // c (data array)
  int,               // number of points [0..np-1]
  int,int,           // nl, nr, 
  int,int);          // ld, m


/* convlv()
   Convolves or deconvolves a real data set data[1..n] (including any 
   user-supplied zero padding) with a response function respns[1..n]. 
   The response function must be stored in wrap-around order in the 
   first m elements of respns, where m is an odd integer <=n. Wrap-
   around order means that the first half of the array respns contains 
   the impulse response function at positive times, while the second 
   half of the array contains the impulse response function at negative 
   times, counting down from the highest element respns[m]. On input 
   isign is +1 for convolution, -1 for deconvolution. The answer is 
   returned in the first n components of ans. However, ans must be 
   supplied in the calling program with dimensions [1..2*n], for 
   consistency with twofft. n MUST be an integer power of two.
*/
void convlv(
  double *,      //data[], 
  int,           //  n, 
  double *,      //respns[]
  int,           // m,
  int,           // isign, 
  double *);     // ans[]


/* sgSmooth()
   Smoothes c[0..n-1] using a Savitzky-Golay filter.
   nl is the number of leftward (past) data points used, while 
   nr is the number of rightward (future) data points, making 
      the total number of data points used nl+nr+1. 
   ld is the order of the derivative desired 
      (e.g., ld = 0 for smoothed function). 
   m  is the order of the smoothing polynomial, 
      also equal to the highest conserved moment; 
      usual values are m = 2 or m = 4.
*/
void sgSmooth (
  double *,          // c (data array)
  int,               // number of points [0..n-1]
  int,int,           // nl, nr, 
  int,int);          // ld, m

#endif
