#ifndef _RANLIB_CC_
#define _RANLIB_CC_

#include "ranlib.h"

//=======================================================================
// Member functions of RANLIB_class
// --------------------------------
// Here all the member functions are provided which are needed for 
// generating uniform random deviates in ]0,1[. Extra functionality 
// can be incorporated by #define RANLIBX which includes the files 
// ranlibx.h and ranlibx.cc.
//
// Author: Michael Schmelling / May 1999
//=======================================================================

RANLIB_class::~RANLIB_class()
//================================
// destructor - return local store
//================================
{
  free(work);
  free(scal_f);
  free(scal_d);
}

//=============================
// instantiate the constructors
//=============================

RANLIB_class::RANLIB_class() { init(); }
RANLIB_class::RANLIB_class(int seed) { init(); iluxinit(seed, 4, work); }
RANLIB_class::RANLIB_class(int seed, int luxl) { init(); iluxinit(seed, luxl, work); }


void RANLIB_class::init()
//===============================================
// default initialization of the random generator
//===============================================
{
  // create working space for the ilux-kernel

  nwork = 34;
  work  = (int*)malloc(nwork*sizeof(int));

  // initialize padding 

  maxpad = 101;
  scal_f = (float*) malloc(maxpad*sizeof(float));
  scal_d = (double*)malloc(maxpad*sizeof(double));
  for(int ipad=0; ipad<maxpad; ipad++) {
    scal_f[ipad] = pow((double) 2,(double) -24-ipad);
    scal_d[ipad] = pow((double) 2,(double) -48-ipad);
  }
  maxpad--;

  // default initialization of work

  iluxinit(314159265, 4, work); 
}



int RANLIB_class::ilux() 
//==================================================
// return a random integer in the range [0,16777215]
//==================================================
{ 
  return ilux(work); 
}


float RANLIB_class::flux()
//===========================================================
// return a single precision random number in the range ]0,1[
//===========================================================
{
  // generate a random integer

  int irn  = ilux(work);
  int ipad = 0;

  // padding up to 24 significant bits

  while(irn<8388608) {
    if(ipad==maxpad) {
      if(irn==0) irn = 1;
      break;
    }
    if(npad==0) {
      rpad = ilux(work); 
      npad = 24;
    }
    ipad++;
    npad--;
    irn  += irn; 
    rpad += rpad;
    if(rpad>16777215) {
      irn++;
      rpad -= 16777216;
    }
  } 

  // float and scale according to the padding

  return irn*scal_f[ipad];
}


double RANLIB_class::dlux()
//===========================================================
// return a double precision random number in the range ]0,1[
//===========================================================
{
  // generate a 48bit random integer (stored as double)
  // nb: decimal points for the integers below are essential!

  double rn   = ilux(work)*16777216. + ilux(work);
  int    ipad = 0;

  // padding up to 48 significant bits

  while(rn<140737488355328.) {
    if(ipad==maxpad) {
      if(rn==0) rn = 1;
      break;
    }
    if(npad==0) {
      rpad = ilux(work);
      npad = 24; 
    }
    ipad++;
    npad--;
    rn   += rn; 
    rpad += rpad;
    if(rpad>16777215) {
      rn   += 1;
      rpad -= 16777216;
    }
  } 

  // scale according to the padding and return

  return rn*scal_d[ipad];
}


int RANLIB_class::ilux(int *ix) 
//===================================================================
// kernel of the random generator
// ------------------------------
// This function is an optimized version of the kernel of the RANLUX
// random generator by F.James as distributed e.g. with CERNLIB 97A.
// It advances the state vector ix[] according to the specified 
// luxury level, returning the next random integer in the sequence. 
// The algebra is realized entirely in integer arithmetic, which 
// on most processors is still faster than floating point operations. 
// For the same initial seeds and luxury level the return values are 
// bit-identical to RANLUX values multiplied by 2^24.  
//===================================================================
{
  // count calls to the ranlux kernel

  minor++;
  if(minor==1000000000) {
    major++;
    minor = 0;
  }

  // bookkeeping - if possible return a number from the buffer

  if(ix[25]--) {
    if(ix[26]--) return ix[ix[26]];
    ix[26]  = -1;
  } else {
    ix[25]  = 23;
    ix[26] -= ix[27];  
  }

  // entire sweeps of subtract with carry

  while(ix[26] < 0) {
    ix[23] = ix[ 9] - ix[23] - ix[24];
    ix[22] = ix[ 8] - ix[22]; if(ix[23]<0) {ix[23] += 16777216; --ix[22];}
    ix[21] = ix[ 7] - ix[21]; if(ix[22]<0) {ix[22] += 16777216; --ix[21];}
    ix[20] = ix[ 6] - ix[20]; if(ix[21]<0) {ix[21] += 16777216; --ix[20];}
    ix[19] = ix[ 5] - ix[19]; if(ix[20]<0) {ix[20] += 16777216; --ix[19];}
    ix[18] = ix[ 4] - ix[18]; if(ix[19]<0) {ix[19] += 16777216; --ix[18];}
    ix[17] = ix[ 3] - ix[17]; if(ix[18]<0) {ix[18] += 16777216; --ix[17];}
    ix[16] = ix[ 2] - ix[16]; if(ix[17]<0) {ix[17] += 16777216; --ix[16];}
    ix[15] = ix[ 1] - ix[15]; if(ix[16]<0) {ix[16] += 16777216; --ix[15];}
    ix[14] = ix[ 0] - ix[14]; if(ix[15]<0) {ix[15] += 16777216; --ix[14];}
    ix[13] = ix[23] - ix[13]; if(ix[14]<0) {ix[14] += 16777216; --ix[13];}
    ix[12] = ix[22] - ix[12]; if(ix[13]<0) {ix[13] += 16777216; --ix[12];}
    ix[11] = ix[21] - ix[11]; if(ix[12]<0) {ix[12] += 16777216; --ix[11];}
    ix[10] = ix[20] - ix[10]; if(ix[11]<0) {ix[11] += 16777216; --ix[10];}
    ix[ 9] = ix[19] - ix[ 9]; if(ix[10]<0) {ix[10] += 16777216; --ix[ 9];}
    ix[ 8] = ix[18] - ix[ 8]; if(ix[ 9]<0) {ix[ 9] += 16777216; --ix[ 8];}
    ix[ 7] = ix[17] - ix[ 7]; if(ix[ 8]<0) {ix[ 8] += 16777216; --ix[ 7];}
    ix[ 6] = ix[16] - ix[ 6]; if(ix[ 7]<0) {ix[ 7] += 16777216; --ix[ 6];}
    ix[ 5] = ix[15] - ix[ 5]; if(ix[ 6]<0) {ix[ 6] += 16777216; --ix[ 5];}
    ix[ 4] = ix[14] - ix[ 4]; if(ix[ 5]<0) {ix[ 5] += 16777216; --ix[ 4];}
    ix[ 3] = ix[13] - ix[ 3]; if(ix[ 4]<0) {ix[ 4] += 16777216; --ix[ 3];}
    ix[ 2] = ix[12] - ix[ 2]; if(ix[ 3]<0) {ix[ 3] += 16777216; --ix[ 2];}
    ix[ 1] = ix[11] - ix[ 1]; if(ix[ 2]<0) {ix[ 2] += 16777216; --ix[ 1];}
    ix[ 0] = ix[10] - ix[ 0]; if(ix[ 1]<0) {ix[ 1] += 16777216; --ix[ 0];}
    ix[24] = 0;               if(ix[ 0]<0) {ix[ 0] += 16777216; ++ix[24];}
    ix[26] += 24;
  }

  return ix[ix[26]];
}


void RANLIB_class::iluxinit(int seed_arg, int luxl_arg, int *ix)
//==============================================================
// initialize working array ix for the ilux-kernel
//==============================================================
{ 
  // set Luescher's p-23 according to the selected lux-level

  if(luxl_arg<=0) ix[27] =   1;
  if(luxl_arg==1) ix[27] =  25;
  if(luxl_arg==2) ix[27] =  74;
  if(luxl_arg==3) ix[27] = 200;
  if(luxl_arg>=4) ix[27] = 366;

  // initialize bookkeeping information

  ix[26] = 0;   // index of current value
  ix[25] = 24;  // bookkeeping of available numbers 
  ix[24] = 0;   // carry 

  // simple linear congruential for initialization of sequence

  seed = seed_arg;
  for(int i=0; i<24; i++) {
    int k = seed/53668;
    seed = 40014*(seed-k*53668) - k*12211;
    if(seed<0) seed += 2147483563;
    ix[i] = seed%16777216;
  } 

  // store the initial seed and luxlevel

  seed = seed_arg;
  luxl = luxl_arg;

  // zero counters and padding buffers

  major = 0;
  minor = 0;
  rpad  = 0;
  npad  = 0;

  // and return

  return;
}


int RANLIB_class::save(const char *filename)
//====================================
// save the generator status to file
//====================================
{
  FILE *file = fopen(filename,"a");
  if(file==NULL) return -1;

  fprintf(file,"RANLIB_class::begin # saving ranlib generator status\n");
  for(int i=0; i<28; i++) {
     fprintf(file,"%10d",work[i]);
     if(i%6==5) fprintf(file,"\n");
  }
  fprintf(file,"%10d%10d\n",rpad,npad);
  fprintf(file,"%10d %d %d %d # seed,luxlevel,major,minor\n",seed,luxl,major,minor);
  fprintf(file,"RANLIB_class::end \n");
  fclose(file);
  return 0;
}


int RANLIB_class::restore(const char *filename)
//=============================================
// restore the generator status from file
//=============================================
{
  char line[256];

  FILE *file = fopen(filename,"r");
  if(file==NULL) return -1;
  
  int index=-1;
  while(fgets(line,256,file)) {
    // cut line at typical comment-escapes 
    if(strstr(line,"#")) *strstr(line,"#") = '\0';
    if(strstr(line,"/")) *strstr(line,"/") = '\0';
    if(strstr(line,"!")) *strstr(line,"!") = '\0';
    if(strstr(line,"*")) *strstr(line,"*") = '\0';
    // read the line
    if(index>-2) {
      char *token = strtok(line," .,;\n\t");
      while(token) {
        index++;
        if(index<nwork) work[index] = atoi(token);
        token=strtok(NULL," .,;\n\t");
      }
    }
    if(index>nwork) break;
    // look for begin/end markers
    if(strstr(line,"RANLIB_class::begin")) index = -1;
    if(strstr(line,"RANLIB_class::end"  )) index = -2;
  }    
  fclose(file);

  if(index != -2) return -1;

  rpad  = work[28];
  npad  = work[29];
  seed  = work[30];
  luxl  = work[31];
  major = work[32];
  minor = work[33];

  return 0;
}


int RANLIB_class::selftest()
//============================================================
// perform a series of checksum tests 
// ----------------------------------
// 1) for all luxury levels and different seeds the first
//    million random-integers produced by the kernel are 
//    compared to the RANLUX output.
// 2) for all luxury levels the save and restore mechanism
//    is verified, i.e. after initialization the status is 
//    saved, then 10000 random numbers are generated by the 
//    kernel, then the status is restored and the same 
//    numbers are generated again.
// 3) The bit-padding of flux() and dlux() is checked for
//    one million calls each.
//============================================================
{
  int    expo, fail;
  double sum, x, y, sum0, sum1, sum2, fint, frac;
  //---------------------------------------------

  fail = 0;
  printf("\n");
  printf("Start RANLIB selftest\n");
  printf("---------------------\n");

  // verify that first 10^6 calls to the kernel conform to RANLUX

  for(int lux=0; lux<5; lux++) {
    if(lux==0) {iluxinit(926531415, lux, work); sum = 8384860421793.;}
    if(lux==1) {iluxinit(265314159, lux, work); sum = 8379462381030.;}
    if(lux==2) {iluxinit(653141592, lux, work); sum = 8383758430272.;}
    if(lux==3) {iluxinit(531415926, lux, work); sum = 8388528261456.;}
    if(lux==4) {iluxinit(314159265, lux, work); sum = 8382862985048.;}
    for(int i=0; i<1000000; i++) sum -= ilux();
    if(sum) fail++;
    if(sum) printf("ranlib  kernel: luxlevel=%d ... FAILED\n",lux);
    else    printf("ranlib  kernel: luxlevel=%d ... OK\n",lux);
  }

  // then check the save and restore mechanism

  for(int lux=0; lux<5; lux++) {
    iluxinit(314159265, lux, work); sum = 0;
    save("selftest.tmp");
    for(int i=0; i<10000; i++) sum += int(16777216.*flux());
    restore("selftest.tmp");
    for(int i=0; i<10000; i++) sum -= int(16777216.*flux());;
    if(sum) fail++;
    if(sum) printf("save / restore: luxlevel=%d ... FAILED\n",lux);
    else    printf("save / restore: luxlevel=%d ... OK\n",lux);
  }
  system("rm -f selftest.tmp");

  // check bit patterns of the uniform float generator

  sum0 = 1001097.;
  sum1 = 0.;
  sum2 = 12581955507227.;

  iluxinit(314159265, 4, work); 
  for(int i=0; i<1000000; i++) {
    x = flux();
    y = frexp(x, &expo);
    frac  = int(16777216.*modf(y*16777216.,&fint));
    sum0 += expo;
    sum1 -= frac;
    sum2 -= fint;
  }
  if(sum0||sum1||sum2) fail++;
  if(sum0||sum1||sum2) printf("float  padding: .............. FAILED\n");
  else                 printf("float  padding: .............. OK\n");

  // check bit patterns of the uniform double generator

  sum0 = 1000496.; 
  sum1 = 8394467436707.; 
  sum2 = 12579877804190.; 

  iluxinit(314159265, 4, work); 
  for(int i=0; i<1000000; i++) {
    x = dlux();
    y = frexp(x, &expo);
    frac  = int(16777216.*modf(y*16777216.,&fint));
    sum0 += expo;
    sum1 -= frac;
    sum2 -= fint;
  }
  if(sum0||sum1||sum2) fail++;
  if(sum0||sum1||sum2) printf("double padding: .............. FAILED\n");
  else                 printf("double padding: .............. OK\n");

  if(fail) printf("Selftest.................. FAILED\n");
  else     printf("Selftest.................. PASSED\n");
  printf("---------------------------------\n\n");

  return fail;
}


//=================================================================
// Extra member functions RANLIB_class
// -----------------------------------
// Here those member functions of RANLIB_class are provided which
// go beyond the generation of uniform random deviates in ]0,1[.
// The prototypes are defined in ranlibx.h. 
//
// Author: Michael Schmelling / May 1999
//=================================================================


//==========================================================
// uniformly distributed random numbers between low and high 
//==========================================================

double RANLIB_class::uniform(double low, double high) 
       {return low + (high-low)*dlux(); }
float  RANLIB_class::uniform(float  low, float  high) 
       {return low + (high-low)*flux(); }

//===============================================================
// gaussian random numbers with mean av and standard deviation sd  
//===============================================================

double RANLIB_class::gauss(double av, double sd) 
       {return av+sd*sqrt(-2*log(dlux()))*cos(6.283185307179586*dlux());}
float  RANLIB_class::gauss(float  av, float  sd) 
       {return av+sd*sqrt(-2*log(flux()))*cos(6.283185307179586*flux());}


float RANLIB_class::rnlat(float a)
//==================================
// dn/dx = (a-1)*(a-2)*x/(1+x)^a
//==================================
{
  static double aold=0;
  static double x, y, b, g, h, d, c, t, u, z;
  //-----------------------------------------

  // setup some useful constants

  if(aold != a) {
    aold = a;
    b =  1 - a;
    g = -1 - b;
    h = -1/g;
    d = (1-g)/(3*b);
    c = sqrt(-2*g/b);
  }

  // get a random number based on an approximate analytic 
  // inversion of the cumulative distribution

  y = flux();
  z = 1 - y;
  z = 1 - c*sqrt(z) + z*(d-z*(1+d-c));
  t = pow(z,h);
  x = t - 1;

  // improve the result by a Newton step

  u = b*x;
  x = x + (1-u-t*y/z)*h*t/u;

  return float(x);
}


double RANLIB_class::rnlat(double a)
//==================================
// dn/dx = (a-1)*(a-2)*x/(1+x)^a
//==================================
{
  static double aold=0;
  static double x, y, b, g, h, d, c, t, u, z;
  //-----------------------------------------

  // setup some useful constants

  if(aold != a) {
    aold = a;
    b =  1 - a;
    g = -1 - b;
    h = -1/g;
    d = (1-g)/(3*b);
    c = sqrt(-2*g/b);
  }

  // get a random number based on an approximate analytic 
  // inversion of the cumulative distribution

  y = dlux();
  z = 1 - y;
  z = 1 - c*sqrt(z) + z*(d-z*(1+d-c));
  t = pow(z,h);
  x = t - 1;

  u = b*x;
  x = x + (1-u-t*y/z)*h*t/u;

  // improve the result by a Newton step

  return x;
}


int RANLIB_class::nbd(double av, double k) 
//==========================================
// generate a negative binomial distribution
//==========================================
{
  int    n = -1, np, nm;
  double xrn, sum, pn, pnp, pnm ;
  double eps = 1.e-10;

  // check arguments

  if(av<=0 || k<=0) {
    printf("invalid arguments av=%f k=%f to rnNBD \n",av,k); 
    return n;
  }

  // if possible, do normal approximation  (avcut=sqrt(ncall), kcut=4*avcut)

  if(av > 10000 && k > 40000 ) {
    double sigma = sqrt(av*(1+av/k)); 
    while(n<0) n = int(0.5 + gauss(av,sigma));
    return n;
  }
 
  // proper generation...

  // first set up some constants

  double r    = av/k;
  double f    = r/(1+r);
  double fk1  = f*(k-1);
  double lnp0 = -k*log(1+r);
  int    nmax = int(av-r-0.5); 

  xrn = dlux();

  // for small averages  start generation at n=0

  if(nmax < 29) {
    n   = 0;
    sum = exp(lnp0);
    pn  = sum;   
    while(sum<xrn) {
      n++;
      pn  *= (f + fk1/n);
      sum += pn;
      if(pn<eps) break;
    }
    return n;
  }

  // do generation around the maximum of the distribution

  n  = nmax;
  pn = gammln(k+n) - gammln(k) - gammln(1+n);
  pn = exp(pn + n*log(f) + lnp0);

  np  = nm  = n;
  pnp = pnm = sum = pn;

  if(sum >= xrn) return n;

  while(1) {
    np++;
    pnp *= (f+fk1/np);
    sum += pnp;
    if(sum >= xrn || pnp<eps) return np;
    if(nm > 0) {
      pnm /= (f + fk1/nm);
      sum += pnm;
      if(sum >= xrn) return nm;
      nm--;
    } 
  }

} 


double RANLIB_class::gammln(double x) 
//====================================================================
// logarithm of gamma function - derived from NRII p207
//====================================================================
{
  double tmp = (x+0.5)*log(x+5.5) - (x+5.5); 
  double ser = 2.506628275107297465e+00       
    + 1.909551718930763968e+02/(x+1) - 2.168366818437280017e+02/(x+2) 
    + 6.019441764023333263e+01/(x+3) - 3.087513239285458511e+00/(x+4) 
    + 3.029638705253258867e-03/(x+5) - 1.352385959072596015e-05/(x+6);
  return tmp + log(ser/x);
}


double RANLIB_class::lgPowerK(double lgxl, double lgxh, double lgxk, double gl, double gh)
//========================================================================================
// generate a random number according to a power-law with a kink
// -------------------------------------------------------------
// Given a distribution 
//                            dn/dx = nl x^gl for x < xk
//                        and dn/dx = nh x^gh for x > xk
//
// with normalizations nl and nh such that the distribution is continuous
// at the location xk of the kink, and spectral indices gl and gh below 
// and above the kink, respectively, then lgPowerK returns log10(x) with
// x distributed according to dn/dx in the range [xl,xh]. Note that the 
// range related arguments are given as log10, e.g. lgxl=log10(xl).
//========================================================================================
{
  static double rlgxl, rlgxh, rgl, rgh, rlgxk;  
  static double xlPg, xhPg, xkPgl, xkPgh, g, dx, dxl, dxh;
  static double ilow, ihig, lfrac;
  static int    nokink; 
  //------------------------------------------------------

  // setup and store some useful constants for a new set of arguments

  if( (lgxl != rlgxl) || (lgxh != rlgxh) || (lgxk != rlgxk) ||   
      (gl   != rgl  ) || (gh   != rgh  )                     ) {

    // store the current arguments

    rlgxl = lgxl; rlgxh = lgxh; rlgxk = lgxk;
    rgl   = gl;   rgh   = gh;

    // prepare the generator

    nokink = 1; 
    if(     lgxh <= lgxk) g = gl + 1; // all below the knee
    else if(lgxl >= lgxk) g = gh + 1; // all above the knee
    else                  nokink = 0; // knee is within window
    
    if(nokink) {
      xlPg = pow(10,g*lgxl); 
      xhPg = pow(10,g*lgxh);
      dx   = xhPg - xlPg;
    }
    else {
      xlPg  = pow(10,(gl+1)*lgxl); 
      xhPg  = pow(10,(gh+1)*lgxh);
      xkPgl = pow(10,(gl+1)*lgxk);
      xkPgh = pow(10,(gh+1)*lgxk);
      ilow  = ( 1 - xlPg/xkPgl)/(gl+1);
      ihig  = (-1 + xhPg/xkPgh)/(gh+1);
      lfrac = ilow/(ilow+ihig); 
      dxl   = xkPgl - xlPg;
      dxh   = xkPgh - xhPg;
    } 
  }

  // generate an lgx-value within the selected range
   
  double xrn = dlux();

  if(nokink)         return log10(xlPg + xrn*dx )/g;
  else {
    if(dlux()<lfrac) return log10(xlPg + xrn*dxl)/(gl+1);
    else             return log10(xhPg + xrn*dxh)/(gh+1);
  }
}
//=================================================================
#endif
