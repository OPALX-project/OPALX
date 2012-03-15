/* savgol.C
   Savitzky-Golay Smoothing Filters

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    ------------
   23/11/2006    Copied from "numerical recipies"                Rene Bakker

   Last Revision:
   $Id: savgol.C 29 2007-04-14 17:03:18Z l_bakker $
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Algorithms/bet/error.h"
#include "Algorithms/bet/math/savgol.h"

/* internal functions
   ====================================================================== */

static double *vector(long n)
/* allocate a double vector with subscript range v[0..n] */
{
  double *v;
  
  v=(double *)malloc((size_t) ((n+1)*sizeof(double)));
  if (!v) {
    //writeError(errModeAll,errFatal,
    //	     "savgol.C allocation failure in vector()");
    writeError();
  }
    return v;
}

int *ivector(long n)
/* allocate an int vector with subscript range v[nl..n] */
{
  int *v;
  
  v=(int *)malloc((size_t) ((n+1)*sizeof(int)));
  if (!v) {
    //writeError(errModeAll,errFatal,
    //	     "savgol.C allocation failure in ivector()");
    writeError();
  }
  return v;
}

static double **matrix(long nr, long nc)
/* allocate a double matrix with subscript range m[0..nr][0..nc] */
{
  long 
    i, 
    nrow=nr+1,
    ncol=nc+1;
  double 
    **m;
  
  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+1)*sizeof(double*)));
  if (!m) {
    //writeError(errModeAll,errFatal,
    //	     "savgol.C allocation failure 1 in matrix()");
    writeError();
  }
  /* allocate rows and set pointers to them */
  m[0]=(double *) malloc((size_t)((nrow*ncol+1)*sizeof(double)));
  if (!m[0]) {
    //writeError(errModeAll,errFatal,
    //		  "savgol.C allocation failure 2 in matrix()");
    writeError();
  }
  for(i=1;i<=nr;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

static void free_matrix(double **a) {
  if (a) {
    if (a[0]) free(a[0]);
    free(a);
  }
}

static int minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

/* lubksb() 
   Solves the set of n linear equations A·X = B. 
   Here a[1..n][1..n] is input, not as the matrix A but rather as its 
   LU decomposition, determined by the routine ludcmp. indx[1..n] is input 
   as the permutation vector returned by ludcmp. b[1..n] is input as the 
   right-hand side vector B, and returns with the solution vector X. a, n, 
   and indx are not modified by this routine and can be left in place for 
   successive calls with different right-hand sides b. This routine takes 
   into account the possibility that b will begin with many zero elements, 
   so it is efficient for use in matrix inversion.
*/
static void lubksb(double **a, int n, int *indx, double b[])
{
  int i,ii=0,ip,j;
  double sum;
  
  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

#define TINY 1.0e-20;

/* ludcmp() 
   Given a matrix a[1..n][1..n], this routine replaces it by the LU 
   decomposition of a rowwise permutation of itself. a and n are input. 
   a is output, arranged as in equation (2.3.14) above;
   indx[1..n] is an output vector that records the row permutation effected 
   by the partial pivoting; d is output as ±1 depending on whether the number 
   of row interchanges was even or odd, respectively. This routine is used 
   in combination with lubksb to solve linear equations or invert a matrix.
*/
static void ludcmp(double **a, int n, int *indx, double *d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv;
  
  vv=vector(n);
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) {
      //writeError(errModeAll,errFatal,
      //	 "Singular matrix in routine ludcmp (savgol)");
      writeError();
    }
      
      vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
  free(vv);
}
#undef TINY

/* savgol()
   Returns in c[1..np], in wrap-around order (N.B.!) 
   consistent with the argument respons in
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
void savgol(double c[], int np, int nl, int nr, int ld, int m)
{
  int imj,ipj,j,k,kk,mm,*indx;
  double d,fac,sum,**a,*b;
  
  if (np < nl+nr+1 || nl < 0 || nr < 0 || ld > m || nl+nr < m) {
    //writeError(errModeAll,errFatal,"bad args in savgol");
    writeError();
  }

  indx=ivector(m+1);
  a=matrix(m+1,m+1);
  b=vector(m+1);
  for (ipj=0;ipj<=(m << 1);ipj++) {
    sum=(ipj ? 0.0 : 1.0);
    for (k=1;k<=nr;k++) sum += pow((double)k,ipj);
    for (k=1;k<=nl;k++) sum += pow((double)-k,ipj);
    mm=FMIN(ipj,2*m-ipj);
    for (imj = -mm;imj<=mm;imj+=2) a[1+(ipj+imj)/2][1+(ipj-imj)/2]=sum;
  }
  ludcmp(a,m+1,indx,&d);
  for (j=1;j<=m+1;j++) b[j]=0.0;
  b[ld+1]=1.0;
  lubksb(a,m+1,indx,b);
  for (kk=1;kk<=np;kk++) c[kk]=0.0;
  for (k = -nl;k<=nr;k++) {
    sum=b[1];
    fac=1.0;
    for (mm=1;mm<=m;mm++) sum += b[mm+1]*(fac *= k);
    kk=((np-k) % np)+1;
    c[kk]=sum;
  }
  free(b);
  free_matrix(a);
  free(indx);
}

#undef FMIN

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/* four1()
   Replaces data[1..2*nn] by its discrete Fourier transform, if isign 
   is input as 1; or replaces data[1..2*nn] by nn times its inverse 
   discrete Fourier transform, if isign is input as -1. data is a 
   complex array of length nn or, equivalently, a real array of length 
   2*nn. nn MUST be an integer power of 2 (this is not checked for!).
*/
void four1(double data[], int nn, int isign)
{
  int n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;
  
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
#undef SWAP

/* realft()
   Calculates the Fourier transform of a set of n real-valued data points. 
   Replaces this data (which is stored in array data[1..n]) by the positive 
   frequency half of its complex Fourier transform. The real-valued first 
   and last components of the complex transform are returned as elements 
   data[1] and data[2], respectively. n must be a power of 2. This routine 
   also calculates the inverse transform of a complex data array if it is 
   the transform of real data. (Result in this case must be multiplied 
   by 2/n.)static
*/
static void realft(double data[], int n, int isign)
{
  void four1(double data[], int nn, int isign);
  int i,i1,i2,i3,i4,np3;
  double c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;
  
  theta=3.141592653589793/(double) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data,n>>1,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np3=n+3;
  for (i=2;i<=(n>>2);i++) {
    i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r-data[2];
  } else {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    four1(data,n>>1,-1);
  }
}

/* twofft()
   Given two real input arrays data1[1..n] and data2[1..n], 
   this routine calls four1 and returns two complex output arrays, 
   fft1[1..2n] and fft2[1..2n], each of complex length n (i.e., 
   real length 2*n), which contain the discrete Fourier transforms 
   of the respective data arrays. n MUST be an integer power of 2
*/
static void twofft(double data1[], double data2[], double fft1[], double fft2[],
	    int n)
{
  void four1(double data[], int nn, int isign);
  int nn3,nn2,jj,j;
  double rep,rem,aip,aim;
  
  nn3=1+(nn2=2+n+n);
  for (j=1,jj=2;j<=n;j++,jj+=2) {
    fft1[jj-1]=data1[j];
    fft1[jj]=data2[j];
  }
  four1(fft1,n,1);
  fft2[1]=fft1[2];
  fft1[2]=fft2[2]=0.0;
  for (j=3;j<=n+1;j+=2) {
    rep=0.5*(fft1[j]+fft1[nn2-j]);
    rem=0.5*(fft1[j]-fft1[nn2-j]);
    aip=0.5*(fft1[j+1]+fft1[nn3-j]);
    aim=0.5*(fft1[j+1]-fft1[nn3-j]);
    fft1[j]=rep;
    fft1[j+1]=aim;
    fft1[nn2-j]=rep;
    fft1[nn3-j] = -aim;
    fft2[j]=aip;
    fft2[j+1] = -rem;
    fft2[nn2-j]=aip;
    fft2[nn3-j]=rem;
  }
}

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

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
void convlv(double data[], int n, double respns[], int m,
		   int isign, double ans[])
{
  int i,no2;
  double dum,mag2,*fft;
  
  fft=vector(n<<1);

  for (i=1;i<=(m-1)/2;i++)
    respns[n+1-i]=respns[m+1-i];
  for (i=(m+3)/2;i<=n-(m-1)/2;i++)
    respns[i]=0.0;
  twofft(data,respns,fft,ans,n);
  no2=n>>1;
  for (i=2;i<=n+2;i+=2) {
    if (isign == 1) {
      ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
      ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
    } else if (isign == -1) {
      if ((mag2=SQR(ans[i-1])+SQR(ans[i])) == 0.0) {
	//	writeError(errModeAll,errFatal,
	//   "savgol.C Deconvolving at response zero in convlv");
	writeError();
      }

      ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
      ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
    } else writeError(); 

      //writeError(errModeAll,errFatal,
    //	      "No meaning for isign in convlv");
  }
  ans[2]=ans[n+1];
  realft(ans,n,-1);

  free(fft);
}

#undef SQR


/* sgSmooth()
   Smoothes c[0..n-1] with a Savitzky-Golay filter.
   nl is the number of leftward (past) data points used, while 
   nr is the number of rightward (future) data points, making 
      the total number of data points used nl+nr+1. 
   ld is the order of the derivative desired 
      (e.g., ld = 0 for smoothed function). 
   m  is the order of the smoothing polynomial, 
      also equal to the highest conserved moment; 
      usual values are m = 2 or m = 4.
*/
void sgSmooth(double *c, int n, int nl, int nr, int ld, int m)
{
  double 
    *cIn,*cOut,
    *cf;
  int    
    isign,
    nn,np,i;

  // make dimension 2^m with m integer
  nn = (int) pow(2.0,(int) (log(1.0*n)/log(2.0))+1);
  

  // memory allocation
  cf   = vector(nn);
  cIn  = vector(nn);
  cOut = vector(nn*2);

  // fill data array
  cIn[0] = 0.0;
  memcpy(&cIn[1],c,sizeof(double)*n);
  for (i=n+1; i<=nn; i++) cIn[i] = 0.0;

  // create filter coefficients
  np = nl + nr + 1;
  if ((np % 2) == 0) ++np;
  savgol(cf,np,nl,nr,ld,m);

  // filter
  isign=1;
  convlv(cIn,nn,cf,np,isign,cOut);

  // move data back
  memcpy(c,&cOut[1],sizeof(double)*n);  

  free(cIn); free(cOut);
  free(cf);
  //  writeError(errModeMaster,errMessage,"SG n = %d, nn = %d, np = %d",n,nn,np);
}
