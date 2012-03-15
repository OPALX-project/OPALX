/* Driver for routine rk4 */

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

void derivs(double x,double y[],double dydx[])
{
  dydx[0] = -y[1];
  dydx[1]=y[0]-(1.0/x)*y[1];
  dydx[2]=y[1]-(2.0/x)*y[2];
  dydx[3]=y[2]-(3.0/x)*y[3];
}

int main(void)
{
  int i,j;
  double h,x=1.0,*y,*yout;
  
  y=vector(N);
  yout=vector(N);
  y[0]=bessj(0,x);
  y[1]=bessj(1,x);
  y[2]=bessj(2,x);
  y[3]=bessj(3,x);
  printf("\n%16s %5s %12s %12s %12s\n",
	 "Bessel function:","j0","j1","j3","j4");
  for (i=0;i<5;i++) {
    h=0.2*i;
    for (j=0; j<4; j++) yout[j] = y[j];
    rk4(yout,N,x,h,derivs);
    printf("\nfor a step size of: %6.2f\n",h);
    printf("%12s","rk4:");
    for (j=0;j<4;j++) printf(" %12.6f",yout[j]);
    printf("\n%12s %12.6f %12.6f %12.6f %12.6f\n","actual:",
	   bessj(0,x+h),bessj(1,x+h),bessj(2,x+h),bessj(3,x+h));
  }
  free(yout);
  free(y);
  return 0;
}
