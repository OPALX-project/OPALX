/* Driver for routine rtsafe */

#include <stdlib.h>
#include <stdio.h>

#include "bessel.h"
#include "error.h"
#include "root.h"

#define N 100
#define NBMAX 20
#define X1 1.0
#define X2 50.0

static double fx(double x)
{
  return bessj(0,x);
}

static void funcd(double x,double *fn, double *df)
{
  *fn=bessj(0,x);
  *df = -bessj(1,x);
}

/* vector
   allocate an array of doubles [0..n-1] and check memory */
static double *vector(int n) {
  double *b;

  b = (double *) malloc(sizeof(double)*n);
  if (!b) {
    writeError(errModeAll,errFatal,"Insufficient memory malloc %d bytes (rk.C)",sizeof(double)*n);
  }
  return b;
} /* vector */

static void zbrak
(
 double (*fx)(double), 
 double x1, double x2, int n, double xb1[],
 double xb2[], int *nb)
{
  int nbb,i;
  double x,fp,fc,dx;
  
  nbb=0;
  dx=(x2-x1)/n;
  fp=(*fx)(x=x1);
  for (i=1;i<=n;i++) {
    fc=(*fx)(x += dx);
    if (fc*fp < 0.0) {
      xb1[++nbb]=x-dx;
      xb2[nbb]=x;
      if(*nb == nbb) return;
      
    }
    fp=fc;
  }
  *nb = nbb;
}

int main(void)
{
  int i,nb=NBMAX;
  double xacc,root,*xb1,*xb2;
  
  initErrorMsg();

  xb1=vector(NBMAX+1);
  xb2=vector(NBMAX+1);
  zbrak(fx,X1,X2,N,xb1,xb2,&nb);
  printf("\nRoots of bessj0:\n");
  printf("%21s %15s\n","x","f(x)");
  for (i=1;i<=nb;i++) {
    xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
    root=findRoot(funcd,xb1[i],xb2[i],xacc);
    printf("root %3d %14.6f %14.6f\n",i,root,fx(root));
  }
  free(xb2);
  free(xb1);
  return 0;
}
