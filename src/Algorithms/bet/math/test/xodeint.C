/* Driver for routine odeint */


#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>

#include "error.h"
#include "bessel.h"
#include "rk.h"

#define N 4

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


int nrhs;   /* counts function evaluations */

static void derivs(double x,double y[],double dydx[])
{
  nrhs++;
  dydx[0] = -y[1];
  dydx[1]=y[0]-(1.0/x)*y[1];
  dydx[2]=y[1]-(2.0/x)*y[2];
  dydx[3]=y[2]-(3.0/x)*y[3];
}

int main(int argc,char **argv)
{
  int i,nbad,nok,odeint_result;
  double 
    eps,
    h1=0.1,hmin=0.0,x1=1.0,x2=10.0,*ystart;
  
  initErrorMsg();

  ystart=vector(N);

  eps           = 1.0e-30;
  odeint_result = 1;
  while (odeint_result == 1) {
    printf("eps = %le\n",eps);
    ystart[0]=bessj(0,x1);
    ystart[1]=bessj(1,x1);
    ystart[2]=bessj(2,x1);
    ystart[3]=bessj(3,x1);
    nrhs=0;
    rkActivateBuffer(100);

    printf("Init done.\n");

    odeint_result = odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs);

    if (odeint_result == 0) {
      printf("\n%s %13s %3d\n","successful steps:"," ",nok);
      printf("%s %20s %3d\n","bad steps:"," ",nbad);
      printf("%s %9s %3d\n","function evaluations:"," ",nrhs);
      
      rkPrintBuffer();
      printf("\n");
      for (i=0; i<N; i++) printf("y[%d] = %14.6f\n",i,ystart[i]);
      printf("x1=%10.4f y=%14.6f\n",x1,bessj(3,x1));
      printf("x2=%10.4f y=%14.6f\n",x2,bessj(3,x2));
    }
    eps *= 10.0;
  }
  free(ystart);
  return 0;
}



