/* Driver for routine savgol */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>

#include "error.h"
#include "savgol.h"

#define NMAX 1000
#define NTEST 6

#define N 16        /* data array size */
#define M 9         /* response function dimension - must be odd */
#define N2 (2*N)

static double *vector(long n)
/* allocate a double vector with subscript range v[0..n] */
{
  double *v;
  
  v=(double *)malloc((size_t) ((n+1)*sizeof(double)));
  if (!v) writeError(errModeAll,errFatal,
		     "savgol.C allocation failure in vector()");
  return v;
}

int xsavgol(int argc,char **argv)
{
  int i,j,m,nl,np,nr;
  double *c,sum;
  static int mtest[NTEST]={2,2,2,2,4,4};
  static int nltest[NTEST]={2,3,4,5,4,5};
  static int nrtest[NTEST]={2,1,0,5,4,5};
  static char *ans[NTEST]={
    "                      -0.086  0.343  0.486  0.343 -0.086",
    "               -0.143  0.171  0.343  0.371  0.257",
    "         0.086 -0.143 -0.086  0.257  0.886",
    " -0.084  0.021  0.103  0.161  0.196  0.207  0.196  0.161  0.103  0.021 -0.084",
    "         0.035 -0.128  0.070  0.315  0.417  0.315  0.070 -0.128  0.035",
    "  0.042 -0.105 -0.023  0.140  0.280  0.333  0.280  0.140 -0.023 -0.105  0.042"};
  
  c = (double *) malloc(sizeof(double)*(NMAX-1));
  printf("M nl nr\n");
  printf("\t\t\tSample Savitzky-Golay Coefficients\n");
  for (i=0;i<NTEST;i++) {
    m=mtest[i];
    nl=nltest[i];
    nr=nrtest[i];
    np=nl+nr+1;
    savgol(c,np,nl,nr,0,m);
    for (sum=0.0,j=0;j<np;j++) sum += c[j];
    printf("%1d %1d %1d\n",m,nl,nr);
    for (j=nl;j<5;j++) printf("%7s"," ");
    for (j=nl-1;j>=0;j--) printf("%7.3f",c[j]);
    for (j=0;j<nr;j++) printf("%7.3f",c[np-j-1]);
    printf("\n");
    printf("Sum = %7.3f\n",sum);
    printf("Compare answer:\n%s\n",ans[i]);
    printf("=================================================================\n");
  }
  free(c);
  return 0;
}


int xconvlv(int argc, char **argv)
{
  unsigned long i,j;
  int isign;
  double cmp,*data,*respns,*resp,*ans;
  
  data=vector(N);
  respns=vector(N);
  resp=vector(N);
  ans=vector(N2);
  for (i=1;i<=N;i++)
    if ((i >= N/2-N/8) && (i <= N/2+N/8))
      data[i]=1.0;
    else
      data[i]=0.0;
  for (i=1;i<=M;i++) {
    if ((i > 2) && (i < 7))
      respns[i]=1.0;
    else
      respns[i]=0.0;
    resp[i]=respns[i];
  }
  isign=1;
  convlv(data,N,resp,M,isign,ans);
  /* compare with a direct convolution */
  printf("%3s %14s %13s\n","i","CONVLV","Expected");
  for (i=1;i<=N;i++) {
    cmp=0.0;
    for (j=1;j<=M/2;j++) {
      cmp += data[((i-j-1+N) % N)+1]*respns[j+1];
      cmp += data[((i+j-1) % N)+1]*respns[M-j+1];
    }
    cmp += data[i]*respns[1];
    printf("%3ld %15.6f %12.6f\n",i,ans[i],cmp);
  }
  free(ans);
  free(resp);
  free(respns);
  free(data);
  return 0;
}

int main(int argc,char **argv) {
  xsavgol(argc,argv);
  xconvlv(argc,argv);
  return 0;
}
