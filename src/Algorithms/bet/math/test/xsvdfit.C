/* Driver for routine svdfit */


#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>

#include "error.h"
#include "svdfit.h"

void fpoly(double x, double p[], int np)
{
  int j;
  
  p[0]=1.0;
  for (j=1;j<np;j++) p[j]=p[j-1]*x;
}

double poly(double x,double p[], int np) {
  double sum = p[0];
  int    i;

  for (i=1; i<np; i++) sum += (p[i]*pow(x,i));
  return sum;
}

void fleg(double x, double pl[], int nl)
{
  int j;
  double twox,f2,f1,d;
  
  pl[0]=1.0;
  pl[1]=x;
  if (nl > 2) {
    twox=2.0*x;
    f2=x;
    d=1.0;
    for (j=2;j<nl;j++) {
      f1=d++;
      f2 += twox;
      pl[j]=(f2*pl[j-1]-f1*pl[j-2])/d;
    }
  }
}

double leg(double x,double p[], int np) {
  return (p[0] + p[1]*x + p[2]*(3.0*x*x-1.0)/2.0 +
	  p[3]*(5.0*pow(x,3)-3.0*x)/2.0 +
	  p[4]*(35.0*pow(x,4)-30.0*x*x +3.0)/8.0);
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum)
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

#define NPT 100
#define SPREAD 0.2
#define NPOL 5

double gasdev(long *idum)
{
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
  long idum=(-911);
  int i;
  double chisq,*x,*y,*sig,*a,*e;
  FILE *f;
  
  initErrorMsg();

  x  =(double *) malloc(sizeof(double)*NPT);
  y  =(double *) malloc(sizeof(double)*NPT);
  sig=(double *) malloc(sizeof(double)*NPT);
  a=(double *) malloc(sizeof(double)*NPOL);
  e=(double *) malloc(sizeof(double)*NPOL);

  for (i=0;i<NPT;i++) {
    x[i]=0.02*i;
    y[i]=1.0+x[i]*(2.0+x[i]*(3.0+x[i]*(4.0+x[i]*5.0)));
    y[i] *= (1.0+SPREAD*gasdev(&idum));
    sig[i]=y[i]*SPREAD;
  }

  for (i=0; i<NPOL; i++) a[i] = 1.0;
  svdfit(x,y,sig,NPT,a,NPOL,e,&chisq,fpoly);
  printf("\npolynomial fit (1):\n\n");
  for (i=1;i<=NPOL;i++)
    printf("%12.6f %s %10.6f\n",a[i],"  +-",e[i]);
  printf("\nChi-squared %12.6f (%12.6f)\n",chisq,sqrt(chisq)/NPT);
  f = fopen("poly1.dat","w");
  for (i=0; i<NPT; i++) 
    fprintf(f,"%12.6f \t %12.6f \t %12.6f\n",
	    x[i],y[i],poly(x[i],a,NPOL));
  fclose(f);

  for (i=0; i<NPOL; i++) a[i] = 1.0;
  svdfitP(x,y,NPT,a,NPOL,e,&chisq);
  printf("\npolynomial fit (2):\n\n");
  for (i=1;i<=NPOL;i++)
    printf("%12.6f %s %10.6f\n",a[i],"  +-",e[i]);
  printf("\nChi-squared %12.6f (%12.6f)\n",chisq,sqrt(chisq)/NPT);
  f = fopen("poly2.dat","w");
  for (i=0; i<NPT; i++) 
    fprintf(f,"%12.6f \t %12.6f \t %12.6f\n",
	    x[i],y[i],poly(x[i],a,NPOL));
  fclose(f);

  for (i=0; i<NPOL; i++) a[i] = 1.0;
  svdfit(x,y,sig,NPT,a,NPOL,e,&chisq,fleg);
  printf("\nLegendre polynomial fit:\n\n");
  for (i=1;i<=NPOL;i++)
    printf("%12.6f %s %10.6f\n",a[i],"  +-",e[i]);
  printf("\nChi-squared %12.6f (%12.6f)\n",chisq,sqrt(chisq)/NPT);
  f = fopen("legp.dat","w");
  for (i=0; i<NPT; i++) 
    fprintf(f,"%12.6f \t %12.6f \t %12.6f\n",
	    x[i],y[i],leg(x[i],a,NPOL));
  fclose(f);

  free(a); free(e);
  free(y);
  free(x);
  return 0;
}
