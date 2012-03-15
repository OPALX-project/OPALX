/* Driver for routine linfit */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>

#include "error.h"
#include "linfit.h"

#define NPT 100
#define SPREAD 0.5


/* vector
   allocate an array of doubles [0..n-1] and check memory */
static double *vector(int n) {
  double *b;

  b = (double *) malloc(sizeof(double)*n);
  if (!b) {
    writeError(errModeAll,errFatal,"Insufficient memory malloc %d bytes (svdfit.C)",sizeof(double)*n);
  }
  return b;
} /* vector */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

static double ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;
  
  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

static double gasdev(long *idum)
{
  double ran1(long *idum);
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}

int main(void)
{
  long idum=(-117);
  int i,mwt;
  double a,b,chi2,q,siga,sigb,*x,*y,*sig;
  
  x=vector(NPT);
  y=vector(NPT);
  sig=vector(NPT);
  for (i=0;i<NPT;i++) {
    x[i]=0.1*i;
    y[i] = -2.0*x[i]+1.0+SPREAD*gasdev(&idum);
    sig[i]=SPREAD;
  }
  for (mwt=-1;mwt<2;mwt++) {
    if (mwt < 0) {
      linfit(x,y,NPT,&a,&b,&siga,&sigb,&chi2);
      printf("\nIgnoring standard deviations (short)\n");
      q = 1.0;
    } else {
      linfit(x,y,NPT,sig,mwt,&a,&b,&siga,&sigb,&chi2,&q);
      if (mwt == 0)
	printf("\nIgnoring standard deviations\n");
      else
	printf("\nIncluding standard deviations\n");
    }
    printf("%12s %9.6f %18s %9.6f \n",
	   "a  =  ",a,"uncertainty:",siga);
    printf("%12s %9.6f %18s %9.6f \n",
	   "b  =  ",b,"uncertainty:",sigb);
    printf("%19s %14.6f \n","chi-squared: ",chi2);
    printf("%23s %10.6f \n","goodness-of-fit: ",q);
  }
  free(sig);
  free(y);
  free(x);
  return 0;
}

