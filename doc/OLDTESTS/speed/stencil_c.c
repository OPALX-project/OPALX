/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#include <stdio.h>
#include <time.h>
#include <malloc.h>

/**********************************************************************/

int c1(int n, double *a,
       double *b, double b0)
{
  int i;
  for (i=0; i<n; ++i)
    a[i] = b[i]*b0;
  return n;
}

int c2(int n, double *a,
       double *b, double b0, double b1)
{
  int i;
  for (i=0; i<n; ++i)
    a[i] = b[i]*b0 + b[i+1]*b1;
  return n*3;
}

int c3(int n, double *a,
       double *b, double b0, double b1, double b2)
{
  int i;
  for (i=0; i<n; ++i)
    a[i] = b[i]*b0 + b[i+1]*b1 + b[i+2]*b2;
  return n*5;
}

int c4(int n, double *a,
       double *b, double b0, double b1, double b2, double b3)
{
  int i;
  for (i=0; i<n; ++i)
    a[i] = b[i]*b0 + b[i+1]*b1 + b[i+2]*b2 + b[i+3]*b3;
  return n*7;
}

int c5(int n, double *a,
       double *b, double b0, double b1, double b2, double b3, double b4)
{
  int i;
  for (i=0; i<n; ++i)
    a[i] = b[i]*b0 + b[i+1]*b1 + b[i+2]*b2 + b[i+3]*b3 + b[i+4]*b4;
  return n*9;
}

int c6(int n, double *a,
       double *b,double b0,double b1,double b2,double b3,double b4,double b5)
{
  int i;
  for (i=0; i<n; ++i)
    a[i] = b[i]*b0 + b[i+1]*b1 + b[i+2]*b2 + b[i+3]*b3 + b[i+4]*b4 + b[i+5]*b5;
  return n*11;
}

/**********************************************************************/

int main(int argc, char *argv[])
{
  int total_its = atoi(argv[1]);
  int max_len = atoi(argv[2]);
  int len=1;
  double *mem = malloc(sizeof(double)*(max_len*5 + 10));
  double *a = mem;
  double *b = a+max_len;
  double *c = b+max_len;
  double *d = c+max_len;

  {
    int i;
    for (i=0;i<max_len;i++)
      mem[i] = sin((double)i);
  }

  printf("%8s %10s %10s %10s %10s %10s %10s\n","length",
	 "mflops 1","mflops 2","mflops 3", "mflops 4","mflops 5","mflops 6");

  while (len<max_len) {
    clock_t tic0,tic1,tic2,tic3,tic4,tic5,tic6;
    double t1,t2,t3,t4,t5,t6;
    double mf1,mf2,mf3,mf4,mf5,mf6;
    int op1,op2,op3,op4,op5,op6;
    int its = total_its / len;
    int it,next_len;

    tic0 = clock();
    for (it=0; it<its; ++it)
      op1 = c1(len,a,b,1.0);
    tic1 = clock();
    for (it=0; it<its; ++it)
      op2 = c2(len,a,b,1.0,2.0);
    tic2 = clock();
    for (it=0; it<its; ++it)
      op3 = c3(len,a,b,1.0,2.0,3.0);
    tic3 = clock();
    for (it=0; it<its; ++it)
      op4 = c4(len,a,b,1.0,2.0,3.0,4.0);
    tic4 = clock();
    for (it=0; it<its; ++it)
      op5 = c5(len,a,b,1.0,2.0,3.0,4.0,5.0);
    tic5 = clock();
    for (it=0; it<its; ++it)
      op6 = c6(len,a,b,1.0,2.0,3.0,4.0,5.0,6.0);
    tic6 = clock();

    t1 = 1e6*(tic1-tic0)/(its*(double)CLOCKS_PER_SEC);
    t2 = 1e6*(tic2-tic1)/(its*(double)CLOCKS_PER_SEC);
    t3 = 1e6*(tic3-tic2)/(its*(double)CLOCKS_PER_SEC);
    t4 = 1e6*(tic4-tic3)/(its*(double)CLOCKS_PER_SEC);
    t5 = 1e6*(tic5-tic4)/(its*(double)CLOCKS_PER_SEC);
    t6 = 1e6*(tic6-tic5)/(its*(double)CLOCKS_PER_SEC);
    mf1 = op1/t1;
    mf2 = op2/t2;
    mf3 = op3/t3;
    mf4 = op4/t4;
    mf5 = op5/t5;
    mf6 = op6/t6;

    printf("%8i %10g %10g %10g %10g %10g %10g\n",len,
	   mf1,mf2,mf3,mf4,mf5,mf6);

    next_len = len*1.1;
    if ( next_len == len )
      len += 1;
    else
      len = next_len;
  }
}
/***************************************************************************
 * $RCSfile: stencil_c.c,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:07:07 $
 * IPPL_VERSION_ID: $Id: stencil_c.c,v 1.1.1.1 2000/11/30 21:07:07 adelmann Exp $ 
 ***************************************************************************/
