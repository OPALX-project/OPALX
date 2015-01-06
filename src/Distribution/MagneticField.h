#ifndef MAGNETICFIELD_H
#define MAGNETICFIELD_H

#include <iostream>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <Ippl.h>

#define MALLOC(n,a) ((a*)malloc((n)*sizeof(a)))

/// Temporary class for comparison with reference program of Dr. Christian Baumgarten.
class MagneticField
{
public:
  /// Returns the interpolated values of the magnetic field: B(r,theta), dB/dr, dB/dtheta
  /*!
   * @param *bint stores the interpolated magnetic field
   * @param *brint stores the interpolated radial derivative of the magnetic field
   * @param *btint stores the interpolated "theta" derivative of the magnetic field
   * @param th is the angle theta where to interpolate
   * @param nr is the number of radial steps
   * @param nth is the number of angle steps per sector
   * @param r is radial position where to interpolate
   * @param r0 is the minimal extent of eh cyclotron
   * @param dr is the radial step size
   * @param **b stores the data of the magnetic field (read in by ReadSectorMap())
   */
  static int interpolate(double*, double*, double*, double, int, int, double, double, double, float**);
  
  /// Returns the interpolated values of the magnetic field: B(r,theta), dB/dr, dB/dtheta
  /*!
   * @param *bint stores the interpolated magnetic field
   * @param *brint stores the interpolated radial derivative of the magnetic field
   * @param *btint stores the interpolated "theta" derivative of the magnetic field
   * @param ith is the angle theta where to interpolate
   * @param nr is the number of radial steps
   * @param nth is the number of angle steps per sector
   * @param r is radial position where to interpolate
   * @param r0 is the minimal extent of eh cyclotron
   * @param dr is the radial step size
   * @param **b stores the data of the magnetic field (read in by ReadSectorMap())
   * @param btflag true if dB/dtheta should be computed, else false
   */
  static int interpolate1(double*,double*,double*,int,int,int,double,double, double, float**,int);
  
  /// Reads the magnetic field data from a file
  /*!
   * @param **b stores the read data
   * @param nr is the number of radial steps
   * @param nth is the number of angle steps per sector
   * @param nsc is a symmetry parameter
   * @param fname is the filename
   * @param sum is 0
   */
  static int ReadSectorMap(float**, int, int, int, std::string, double);
  static void MakeNFoldSymmetric(float**, int, int, int, int);
  static float** malloc2df(int, int);
  
private:
  // ipolmethod==0 => lineare Interpolation
  static int ipolmethod;
  
  static void intcf(double, std::vector<double>&);
  static void _intcf(double, std::vector<double>&, int);
  static double norm360(double);
  
};

int MagneticField::interpolate(double* bint,double* brint,double* btint,double th,int nr,int nth,
	   double r,double r0, double dr, float** b) {
    double dth,frb,fth;
    std::vector<double> cofr(8);
    std::vector<double> cofth(8);
    int k0,k1,nrb,i,i1,j,ith,j1,xth,ctrl=0;
    *bint=0.0;
    *brint=0.0;
    *btint=0.0;
    // Map "th" into interval [0..360[: 
    th=norm360(th);
    // Map radius "r" onto radial grid:
    if (r<r0) return 0;
    frb=(r-r0)/dr;
    // Split radius in fractional and integer value for interpolation:
    nrb=(int)frb;
    frb-=(double)nrb;
    dth=360./(double)nth;
    fth=th/dth;
    ith=(int)fth;
    fth-=(double)ith;
    //     r-interpolation for b-field;
    ctrl=0;
    k0=0;
    k1=4;
    if (nrb+1>nr) return 2;
    if (nrb+1==nr) {
	/* Extrapolation: */
	ctrl=2;
	k0=0;
	k1=3;
	nrb--;
    } else
    if (nrb+2==nr) {
	ctrl=1;
	k0=0;
	k1=3;	
    }
    // Calculate coefficients for Laplace-Interpolation:
    _intcf(frb,cofr,ctrl);
    intcf(fth,cofth);
    
    // Sum 'em up:
    for (i=k0;i<k1;i++) {
      xth=0;
      i1=std::abs(nrb+i-1);
      if (i1<0) {
	i1=std::abs(i1);
	xth=nth/2;
      }
      for (j=0;j<4;j++) { 
	j1=(ith+j-1+xth) % nth;
	if (j1<0) j1+=nth;
	*bint+=cofr[i]*cofth[j]*b[j1][i1];	// problem with b, since size not 1440
	*brint+=cofr[i+4]*cofth[j]*b[j1][i1];
	*btint+=cofr[i]*cofth[j+4]*b[j1][i1];
      }
    }
    
    *btint*=((double)nth)/(2.0*M_PI);
    *brint/=dr;
    return 0;
}

int MagneticField::interpolate1(double *bint,double *brint,double *btint,int ith,int nr,int nth,
	   double r,double r0, double dr, float **b,int btflag) {
    double frb,dth;
    
    std::vector<double> cof(8);
    
    int nrb,i,i0,i1,ith1,ith2,ctrl;
    // Map "ith" into interval [0..nth-1]: 
    while (ith>=nth) ith-=nth;
    while (ith<0) ith+=nth;
    // Map radius "r" onto radial grid:
    frb=(r-r0)/dr-1.0;
    // Split radius in fractional and integer value for interpolation:
    nrb=(int)frb;
    frb-=(double)nrb;
    ctrl=0;
    i0=0;
    i1=4;
    if (nrb<0) {
	ctrl=-1;
	i0=1;
	i1=4;
    }
    if (nrb+3==nr) {
	ctrl=1;
	i0=0;
	i1=3;	
    }
    if (nrb+2==nr) {
	ctrl=2;
	i0=0;
	i1=3;
	nrb--;
    }
    if (nrb+2>nr) {
	return 2;
    }
    //     r-interpolation for b-field;
    // Calculate coefficients for Laplace-Interpolation:
    _intcf(frb,cof,ctrl);
    *bint=0.0;
    *brint=0.0;
    // Sum 'em up:
    for (i=i0;i<i1;i++) { 
      *bint+=cof[i]*b[ith][nrb+i];
      *brint+=cof[i+4]*b[ith][nrb+i];
    }
    (*brint)/=dr;
    // Theta derivative?
    if (btflag) {
	*btint=0.0;
	dth=2.0*M_PI/(double)nth;
	ith1=ith-1;
	if (ith1<0) ith1+=nth;
	ith2=ith+1;
	if (ith2>=nth) ith2-=nth;
	for (i=i0;i<i1;i++) {
	  *btint+=cof[i]*(b[ith2][nrb+i]-b[ith1][nrb+i]);
	}
	(*btint)/=(2.0*dth);
    }
    return 0;
}

// Sum=0: Read Field Map with nr,nth entries from file "fname" into array "b"
// Sum=1: Add Field Map with nr,nth entries from file "fname" to array "b"
// b=b[nth*nsc][nr];
// Return value is BOOL(success).
int MagneticField::ReadSectorMap(float** b,int nr,int nth,int nsc, std::string fname, double sum=0) {
    double tmp;
    int i,j,k,err=0;
    FILE *f;
    
    f = fopen(fname.c_str(),"r");
    if (f != NULL) {
	j=0;
	while ((j<nr)&&(err==0)) {		// j: distance (radius)
	    for (i=0;i<nth;i++) {		// i: angles
		err=(fscanf(f,"%lf",&tmp)==EOF);
		if (!err) {
		    if (sum) {
			for (k=0;k<nsc;k++)
			    b[i+nth*k][j]+=sum*tmp;	// b is a (1440 x 8) x 141 matrix
		    } else {
			for (k=0;k<nsc;k++)
			    b[i+nth*k][j]=tmp;
		    }
		} else fprintf(stderr,"Error reading %s at (%d,%d)\n",fname.c_str(),j,i);

	    }
	    if (!err) j++;
	}
	fclose(f);
	if (!err) INFOMSG("Sector map " << fname << " read in" << endl);
	return j;
    } else ERRORMSG("Can not open " << fname << endl);
    return 0;
}

float** MagneticField::malloc2df(int n,int m) {
    float **a;
    int i;
    a=MALLOC(n,float*);
    for (i=0;i<n;i++) {
	a[i]=MALLOC(m,float);
	std::memset(a[i],0,m*sizeof(float));
    }
    return a;
}

int MagneticField::ipolmethod=1;

void MagneticField::intcf(double x, std::vector<double>& cof) {
  //     double 3-pt interpolation for general use;
  //      dimension cof(8);
  double x2,x3;
  
  if (ipolmethod) {
    x2=x*x;
    x3=x2*x;
    cof[0]=x2-.5*(x3+x);
    cof[1]=1.5*x3-2.5*x2+1.;
    cof[2]=2.*x2-1.5*x3+.5*x;
    cof[3]=.5*(x3-x2);
    cof[4]=2.*x-1.5*x2-.5;
    cof[5]=4.5*x2-5.*x;
    cof[6]=4.*x-4.5*x2+.5;
    cof[7]=1.5*x2-x;
  } else {
    cof[0]=0.0;
    cof[1]=1.0-x;
    cof[2]=x;
    cof[3]=0.0;
    cof[4]=0.0;
    cof[5]=-1.0;
    cof[6]=1.0;
    cof[7]=0.0;
  }
}

void MagneticField::_intcf(double x, std::vector<double>& cof, int ctrl) {
  double x2=x*x;
  switch (ctrl) {
    case 0:	intcf(x,cof); break;
    case -1:
      cof[0]=0.0;
      cof[1]=1.-1.5*x+0.5*x2;
      cof[2]=2.*x-x2;
      cof[3]=0.5*(x2-x);
      cof[4]=0.0;
      cof[5]=-1.5+x;
      cof[6]=2.*(1.-x);
      cof[7]=x-0.5;
      break;
    case 1:
      cof[0]=0.5*(x2-x);
      cof[1]=1.-x2;
      cof[2]=0.5*(x+x2);
      cof[3]=0.0;
      cof[4]=x-0.5;
      cof[5]=-2.0*x;
      cof[6]=0.5+x;
      cof[7]=0.0;
      break;
    case 2: /* extrapolating: */
      cof[0]=0.5*(x+x2);
      cof[1]=-2.0*x-x2;
      cof[2]=1.+1.5*x+0.5*x2;
      cof[3]=0.0;
      cof[4]=0.5+x;
      cof[5]=-2.0*(1.+x);
      cof[6]=1.5+x;
      cof[7]=0.0;
      break;
    }
}

double MagneticField::norm360(double phi) {
    double ph=std::fmod(phi,360.0);
    if (ph<0.0) ph+=360.0;
    if (std::fabs(ph-360.0)<1E-3) ph=0.0;
    return ph;
}

void MagneticField::MakeNFoldSymmetric(float** bmag, int N, int nr, int nth, int nsc) {
  double bav;
  for (int i=0;i<nr;i++) {
    for (int j=0;j<nth;j++) {
      bav=0.0;
      for (int k=0;k<nsc;k++) bav+=bmag[(j+k*nth) % N][i];
      bav/=(double)nsc;
      for (int k=0;k<nsc;k++) bmag[(j+k*nth) % N][i]=bav;
    }
  }
}

#endif
