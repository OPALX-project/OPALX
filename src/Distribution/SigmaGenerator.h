#ifndef SIGMAGENERATOR_H
#define SIGMAGENERATOR_H

#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <list>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

#include "error.h"
#include "math.h"
#include "physics.h"
#include "rdm.h"

#include "matrix_vector_operation.h"
#include "ClosedOrbitFinder.h"
#include "MapGenerator.h"

#include <boost/numeric/ublas/io.hpp>

#ifdef PRINT
#include <fstream>
#endif

#ifdef DEBUG
#include <iomanip>
#endif

/// This class computes the matched distribution
/*!
 * It uses the class <b>ClosedOrbitFinder</b> to get the parameters (inverse bending radius, path length
 * field index and tunes) to initialize the sigma matrix. It further has a private class member of type
 * <b>RDM</b> for decoupling.
 * The main function of this class is <b>match(value_type, size_type)</b>, where it iteratively tries to find a matched distribution for given
 * emittances, energy and current. The computation stops when the L2-norm is smaller than a user-defined tolerance.
 */

template<typename Value_type, typename Size_type>
class SigmaGenerator
{
public:
  /// Type of variables
  typedef Value_type value_type;
  /// Type of constant variables
  typedef const value_type const_value_type;
  /// Type for specifying sizes
  typedef Size_type size_type;
  /// Type for storing maps
  typedef boost::numeric::ublas::matrix<value_type> matrix_type;
  /// Type for storing the sparse maps
  typedef boost::numeric::ublas::compressed_matrix<value_type,boost::numeric::ublas::row_major> sparse_matrix_type;
  /// Type for storing vectors
  typedef boost::numeric::ublas::vector<value_type> vector_type;
  /// Container for storing the properties for each angle
  typedef std::vector<value_type> container_type;
  
  /// Constructs an object of type SigmaGenerator
  /*!
   * @param I specifies the current for which a matched distribution should be found, \f$ [I] = mA \f$
   * @param ex is the emittance in x-direction (horizontal), \f$ \left[\varepsilon_{x}\right] = \pi\ mm\ mrad  \f$
   * @param ey is the emittance in y-direction (longitudinal), \f$ \left[\varepsilon_{y}\right] = \pi\ mm\ mrad  \f$
   * @param ez is the emittance in z-direction (vertical), \f$ \left[\varepsilon_{z}\right] = \pi\ mm\ mrad  \f$
   * @param wo is the orbital frequency, \f$ \left[\omega_{o}\right] = \frac{rad}{s} \f$
   * @param E is the energy, \f$ [E] = MeV \f$
   * @param nh is the harmonic number
   * @param m is the mass of the particles \f$ [m] = kg \f$
   * @param Emin is the minimum energy [MeV] needed in cyclotron, \f$ \left[E_{min}\right] = MeV \f$
   * @param Emax is the maximum energy [MeV] reached in cyclotron, \f$ \left[E_{max}\right] = MeV \f$
   * @param nSymmetry is the number of sectors (symmetry assumption)
   * @param N is the number of angle splits
   * @param fieldmap is the location of the file that specifies the magnetic field
   */
  SigmaGenerator(value_type, value_type, value_type, value_type, value_type, value_type, value_type, value_type, value_type, value_type, size_type, size_type, const std::string&);
  
  /// Searches for a matched distribution. Returns "true" if a matched distribution within the accuracy could be found, returns "false" otherwise
  /*!
   * @param accuracy is used for computing the equilibrium orbit (to a certain accuracy)
   * @param maxit is the maximum number of iterations (the program stops if within this value no stationary distribution was found)
   * @param maxitOrbit is the maximum number of iterations for finding closed orbit
   * @param order is the order of the Taylor approximation for constructing the cyclotron map out of the force matrix
   */
  bool match(value_type, size_type, size_type, size_type);
  
  /// Block diagonalizes the symplex part of the one turn transfer matrix
  /** It computes the transformation matrix <b>R</b> and its inverse <b>invR</b>. It returns a vector containing the four eigenvalues
   * (alpha,beta,gamma,delta) (see formula (38), paper: Geometrical method of decoupling)
   */
  /*!
   * @param Mturn is a 6x6 dimensional one turn transfer matrix
   * @param R is the 4x4 dimensional transformation matrix (gets computed)
   * @param invR is the 4x4 dimensional inverse transformation (gets computed)
   */
  vector_type decouple(const matrix_type&, sparse_matrix_type&, sparse_matrix_type&);
  
  /// Checks if the sigma-matrix is an eigenellipse and returns the L2 error (The idea of this function is taken from Dr. Christian Baumgarten's program).
  /*!
   * @param Mturn is the one turn transfer matrix
   * @param sigma is the sigma matrix to be tested
   */
  value_type isEigenEllipse(const matrix_type&, const matrix_type&);
  
  /// Returns the converged stationary distribution
  matrix_type& getSigma();
  
  /// Returns the number of iterations needed for convergence (if not converged, it returns zero)
  size_type getIterations();
  
private:
  /// Stores the value of the current
  value_type I_;
  /// Stores the desired emittances
  std::array<value_type,3> emittance_;
  /// Is the orbital frequency
  value_type wo_;
  /// Stores the user-define energy
  value_type E_;
  /// Relativistic factor (which can be computed out ot the kinetic energy and rest mass (potential energy))
  value_type gamma_;
  /// Relativistic factor squared
  value_type gamma2_;
  /// Harmonic number
  value_type nh_;
  /// Velocity (c/v)
  value_type beta_;
  /// Is the mass of the particles
  value_type m_;
  /// Is the number of iterations needed for convergence
  size_type niterations_;
  /// Is true if converged, false otherwise
  bool converged_;
  /// Minimum energy needed in cyclotron
  value_type Emin_;
  /// Maximum energy reached in cyclotron
  value_type Emax_;
  /// Number of (symmetric) sectors
  size_type nSymmetry_;
  /// Number of angle splits
  size_type N_;
  /// Number of angle splits per sector
  Size_type nSector_;
  /// Location of magnetic fieldmap
  std::string fieldmap_;
  
  /// Stores the stationary distribution (sigma matrix)
  matrix_type sigma_;
  
  /// Vector storing the second moment matrix for each angle
  std::vector<matrix_type> sigmas_;
  
  /// <b>RDM</b>-class member used for decoupling
  RDM<value_type, size_type> rdm_;
  
  /// Initializes a first guess of the sigma matrix with the assumption of a spherical symmetric beam (ex = ey = ez). For each angle split the same initial guess is taken.
  /*!
   * @param nuz is the vertical tune
   * @param ravg is the average radius of the closed orbit
   */
#ifdef DEBUG
public:
#endif
  void initialize(value_type, value_type);
  void chris_initialize(value_type, value_type, value_type);
  
#ifdef DEBUG
public:
#endif
  /// Reduces the 6x6 matrix to a 4x4 matrix (for more explanations consider the source code comments)
  template<class matrix>
  void reduce(matrix&);
  /// Expands the 4x4 matrix back to dimension 6x6 (for more explanations consider the source code comments)
  template<class matrix>
  void expand(matrix&);
  
  /// Computes the new initial sigma matrix
  /*!
   * @param M is the 6x6 one turn transfer matrix
   * @param eigen stores the 4 eigenvalues: alpha, beta, gamma, delta (in this order)
   * @param R is the transformation matrix
   * @param invR is the inverse transformation matrix
   */
  matrix_type updateInitialSigma(const matrix_type&, const vector_type&, sparse_matrix_type&, sparse_matrix_type&);
#ifdef DEBUG
private:
#endif
  /// Computes new sigma matrices (one for each angle)
  /*!
   * Mscs is a vector of all space charge maps
   * Mcycs is a vector of all cyclotron maps
   */
  void updateSigma(const std::vector<matrix_type>&, const std::vector<matrix_type>&);
  
  /// Returns the L2-error norm between the old and new sigma-matrix
  /*!
   * @param oldS is the old sigma matrix
   * @param newS is the new sigma matrix
   */
#ifdef DEBUG
public:
#endif
  value_type L2ErrorNorm(const matrix_type&, const matrix_type&);
};

// -----------------------------------------------------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type>
SigmaGenerator<Value_type, Size_type>::SigmaGenerator(value_type I, value_type ex, value_type ey, value_type ez, value_type wo, 
			       value_type E, value_type nh, value_type m, value_type Emin, value_type Emax, size_type nSymmetry, size_type N, const std::string& fieldmap)
: I_(I), wo_(wo), E_(E),
  gamma_(E/physics::E0+1.0), gamma2_(gamma_*gamma_),
  nh_(nh), beta_(std::sqrt(1.0-1.0/gamma2_)), m_(m), niterations_(0), converged_(false),
  Emin_(Emin), Emax_(Emax), nSymmetry_(nSymmetry), N_(N), nSector_(N/nSymmetry), fieldmap_(fieldmap), sigmas_(nSector_) 
{
  // set emittances (initialization like that due to old compiler version)
  emittance_[0] = ex;
  emittance_[1] = ey;
  emittance_[2] = ez;
}

template<typename Value_type, typename Size_type>
bool SigmaGenerator<Value_type, Size_type>::match(value_type accuracy, size_type maxit, size_type maxitOrbit, size_type order) {
  /* compute the equilibrium orbit for energy E_
   * and get the the following properties:
   * - inverse bending radius h
   * - step sizes of path ds
   * - tune nuz
   */
  
  bool flag = true; // true: my functions; false: Christian's functions
  ClosedOrbitFinder<value_type, size_type, boost::numeric::odeint::runge_kutta4<container_type> > cof(E_,wo_,N_,accuracy,maxitOrbit,Emin_,Emax_,flag,fieldmap_);
  
  container_type h = cof.getInverseBendingRadius();
  container_type r = cof.getOrbit();
  container_type ds = cof.getPathLength();
  container_type fidx = cof.getFieldIndex();
  value_type nur, nuz;
  std::tie(nur,nuz) = cof.getTunes(flag);
  
#ifdef DEBUG
  std::cout << "nur = " << nur << " nuz = " << nuz << " ravg = " << cof.getAverageRadius() << std::endl; std::cin.get();
#endif
  
  double bgam = std::sqrt(cof.getGamma()*cof.getGamma()-1);
  emittance_[0] *= bgam*M_PI;
  emittance_[1] *= bgam*M_PI;
  emittance_[2] *= bgam*M_PI;
#ifdef DEBUG
  std::cout << "Normalized emittances: (" << emittance_[0]/M_PI << ", " << emittance_[1]/M_PI;
  std::cout << ", " << emittance_[2]/M_PI << ")" << std::endl;
#endif
  // initialize sigma matrices (for each angle one) (first guess)
  initialize(nuz,cof.getAverageRadius());
  
#ifdef DEBUG
  std::cout << "Initial Sigma: " << std::endl;
  std::cout << sigmas_[0] << std::endl; std::cin.get();
#endif
  
  // calculate only for a single sector (a nSymmetry_-th) of the whole cyclotron
//   size_type sector = N_/nSymmetry_;	// --> nth = 180
  
#ifdef DEBUG
  std::cout << "Number of sectors: " << nSymmetry_ << std::endl << "Number of angle slices: " << nSector_ << std::endl;
  std::cout << "Set N = nth <--> sector" << std::endl;
#endif
  
  // object for space charge map and cyclotron map
  MapGenerator<value_type, size_type> mapgen(nSector_,I_,gamma_,wo_,nh_);
  
  // compute cyclotron map and space charge map for each angle and store them into a vector
  std::vector<typename MapGenerator<value_type, size_type>::matrix_type> Mcycs(nSector_), Mscs(nSector_);
  
  for(size_type i=0; i<nSector_; ++i) {
    Mcycs[i] = mapgen.generateMcyc(h[i],fidx[i],ds[i],order);
    Mscs[i] = mapgen.generateMsc(sigmas_[i],ds[i]);
  }
  
  // one turn matrix
  mapgen.combine(Mscs,Mcycs);
  matrix_type Mturn = mapgen.getMap();
  
  // (inverse) transformation matrix
  sparse_matrix_type R, invR;
  
  // eigenvalues
  vector_type eigen(4);
  
  matrix_type newSigma(6,6);
  
  bool stop = false;
  
  // initialize the error to be the maximal possible value of datatype "value_type"
  value_type error = std::numeric_limits<value_type>::max();
  
  while(error > accuracy && !stop) {
    // decouple transfer matrix and compute (inverse) tranformation matrix
    eigen = decouple(Mturn,R,invR);
    
    // construct new initial sigma-matrix
    newSigma = updateInitialSigma(Mturn,eigen,R,invR);
#ifdef DEBUG
    std::cout << "newSigma:" << std::endl << newSigma << std::endl; std::cin.get();
#endif
    
    // compute new sigma matrices for all angles (except for initial sigma)
    updateSigma(Mscs,Mcycs);
    
    // compute error
    error = L2ErrorNorm(sigmas_[0],newSigma);
#ifdef DEBUG    
    std::cerr << "error computed: " << error << " " << accuracy << std::endl;
#endif
    // write new initial sigma-matrix into vector
    sigmas_[0] = newSigma;
    
    // compute new space charge maps
    for(size_type i=0; i<nSector_; ++i) {
      Mscs[i] = mapgen.generateMsc(sigmas_[i],ds[i]);
    }
    
    // construct new one turn transfer matrix M
    mapgen.combine(Mscs,Mcycs);
    Mturn = mapgen.getMap();
    
    // check if number of iterations has maxit exceeded.
    stop = (niterations_++ > maxit);
  }
  
  // store converged sigma-matrix
  sigma_.resize(6,6,false);
  sigma_.swap(newSigma);
  
  // returns if the sigma matrix has converged
  converged_ = error<accuracy;
  return converged_;
}

template<typename Value_type, typename Size_type>
typename SigmaGenerator<Value_type, Size_type>::vector_type SigmaGenerator<Value_type, Size_type>::decouple(const matrix_type& Mturn, sparse_matrix_type& R,
										      sparse_matrix_type& invR) {  
  // copy one turn matrix
  matrix_type M(Mturn);
  
  // reduce 6x6 matrix to 4x4 matrix
  reduce<matrix_type>(M);
  
  // compute symplex part
  matrix_type Ms = rdm_.symplex(M);
  
  // diagonalize and compute transformation matrices
  rdm_.diagonalize(Ms,R,invR);
  
  /*
   * formula (38) in paper of Dr. Christian Baumgarten:
   * Geometrical method of decoupling
   * 
   * 		[0, 	alpha, 	0, 	0;
   * F_{d} =	-beta, 	0, 	0, 	0;
   * 		0, 	0, 	0, 	gamma;
   * 		0, 	0, 	-delta,	0]
   * 
   * 
   */
  vector_type eigen(4);
  eigen(0) = Ms(0,1);		// alpha
  eigen(1) = -Ms(1,0);		// beta
  eigen(2) = Ms(2,3);		// gamma
  eigen(3) = -Ms(3,2);		// delta
  return eigen;
}

template<typename Value_type, typename Size_type>
typename SigmaGenerator<Value_type, Size_type>::value_type SigmaGenerator<Value_type, Size_type>::isEigenEllipse(const matrix_type& Mturn, const matrix_type& sigma) {
  // compute sigma matrix after one turn
  matrix_type newSigma = matt_boost::gemmm<matrix_type>(Mturn,sigma,boost::numeric::ublas::trans(Mturn));
  
  // return L2 error
  return L2ErrorNorm(sigma,newSigma);
}

template<typename Value_type, typename Size_type>
inline typename SigmaGenerator<Value_type, Size_type>::matrix_type& SigmaGenerator<Value_type, Size_type>::getSigma() {
  return sigma_;
}

template<typename Value_type, typename Size_type>
inline typename SigmaGenerator<Value_type, Size_type>::size_type SigmaGenerator<Value_type, Size_type>::getIterations() {
  return (converged_) ? niterations_ : size_type(0);
}

// -----------------------------------------------------------------------------------------------------------------------
// PRIVATE MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type>
void SigmaGenerator<Value_type, Size_type>::chris_initialize(value_type ravg, value_type nuz, value_type nur) {
  
#ifdef PRINT
  std::ofstream out("data/maps/InitialSigmaPerAngle.dat");
#endif
  
  // formula (57)
  value_type nurf = wo_*nh_/(2.0*M_PI);		// nurf*1e-6 MHz
  value_type K3=3.0*physics::q0*I_/(20.0*std::sqrt(5.0)*M_PI*physics::eps0*1.0e6*physics::E0*nurf*beta_*beta_*gamma_*gamma2_);  
  
#ifdef DEBUG
  std::cout << "gamma = " << gamma_ << std::endl;
  std::cout << "beta = " << beta_ << std::endl;
  std::cout << "current = " << I_ << std::endl;
  std::cout << "K3 = " << K3*1e8 << std::endl;
#endif
  
  // helper constants
  value_type invbg = 1.0/(beta_*gamma_);
  value_type micro = 1.0e-6;
  value_type inv3gam = 1.0/(3.0*gamma_);
  
  
  /// WHY DON'T WE UPDATE THE EMITTANCE FOR OTHER ENERGIES ???
  emittance_[0] = emittance_[1] = emittance_[2] = 1.254129758715586;
  
  // mm mrad --> meter radian
  value_type ex = emittance_[0]*invbg*micro;	// [ex] = Pi mm mrad
  value_type ey = emittance_[1]*invbg*micro;	// [ey] = Pi mm mrad
  value_type ez = emittance_[2]*invbg*micro;	// [ez] = Pi mm mrad
  
#ifdef DEBUG
  std::cout << "beta*gamma = " << beta_*gamma_ << std::endl;
  std::cout << emittance_[0] << " " << emittance_[1] << " " << emittance_[2] << std::endl;
  std::cout << "chris_initialize: ex = " << ex*1e6 << " ey = " << ey*1e6 << " ez = " << ez*1e6 << std::endl;
#endif
  value_type exz = ex + ez;
  
  value_type f;
  value_type KXY, Kx, Ky, Kz;
  value_type a, a2, b;
  value_type W, w, W2, w2;
  value_type A, B;
  value_type dx, dy, dz, fac;
  
  value_type h = 1.0/ravg;
  value_type kx = h*h*nur*nur;
  value_type ky = h*h*nuz*nuz;
  
  value_type sx = std::sqrt(exz/(h*nur));
  value_type sy = std::sqrt(ey/(h*nuz));
  value_type sz = std::sqrt(ex/(h*nur*gamma2_));
  
  
#ifdef DEBUG
  std::cout << "sx = " << sx << " sy = " << sy << " sz = " << sz << std::endl;
#endif
      
  do {
    // formula (21)
    f = inv3gam*std::sqrt(sx*sy)/sz;
    KXY = K3*std::fabs(1.0-f)/((sx+sy)*sz);
    Kx = KXY/sx;
    Ky = KXY/sy;
    Kz = K3*f*gamma2_/(sx*sy*sz);
    
    // formula (22)
    a = 0.5*(kx-Kx-Kz);
    a2 = a*a;
    
    // formula (23) with h*depsilon/dr = 0
    b = Kz*Kx;
    
    if(a2<b) {
      sx *= 1.41;
      sy *= 1.1;
      sz *= 1.41;
    }
  } while(a2<b);
  
  if(b<a2) {
    value_type invAB, AB;
    value_type err;
    value_type olderr = std::sqrt(3.0);
    value_type f0 = 5.0;
    value_type erx, ery, erz;
    
    do {
      // formula (22)
      value_type tmp = std::sqrt(a2-b);
      W2 = a+tmp;
      w2 = a-tmp;
      
      if(w2>0) {
	W = std::sqrt(W2);
	w = std::sqrt(w2);
	
	// formula (26)
	A = h/(W2+Kz);
	B = h/(w2+Kz);
	
	invAB = 1.0/(B-A);
	AB = A*B;
	
	// formula (30) (std::sqrt ...)
	dx = sx - std::sqrt((B*ex/W+A*ez/w)*invAB);
	dy = sy - std::sqrt(ey/std::sqrt(ky-Ky));
	dz = sz - std::sqrt(AB*(A*ex*W+B*ez*w)*invAB);
	
	erx = dx*dx/(sx*sx);
	ery = dy*dy/(sy*sy);
	erz = dz*dz/(sz*sz);
	
// 	std::cerr << "dx = " << dx << " " << dy << " " << dz << " " << erx << " " << ery << " " << erz << std::endl; std::cin.get();
	
	err = std::sqrt(erx+ery+erz);
	
	if(err>olderr) {
	  f0 *= 1.25;
	}
	if(err<5e-3 && err/olderr>0.9) {
	  f0 /= 1.05;
	}
	olderr = err;
	
	fac = std::exp(-w2/W2)/f0;
	
	sx -= dx*fac;
	sy -= dy*fac;
	sz -= dz*fac;
	
      } else {
	PhysicalError::message("SigmaGenerator::chris_initialize()", PhysicalError::negative);
      }
      
      // formula (21)
      f = inv3gam*std::sqrt(sx*sy)/sz;
      KXY = K3*std::fabs(1.0-f)/((sx+sy)*sz);
      Kx = KXY/sx;
      Ky = KXY/sy;
      Kz = K3*f*gamma2_/(sx*sy*sz);    
      
      // formula (22)
      a = 0.5*(kx-Kx-Kz);
      a2 = a*a;
      // formula (23) with h*depsilon/dr = 0
      b = Kz*Kx;
    } while(err>1e-5);
    
    // meter --> millimeter
    sx *= 1000.0;
    sy *= 1000.0;
    sz *= 1000.0;
    
#ifdef DEBUG
  std::cout << "sx = " << sx << " sy = " << sy << " sz = " << sz << std::endl; std::cin.get();
#endif
    
//     std::cerr << "sx = " << sx << " " << sy << " " << sz << " " << AB << " " << invAB << " " << W << " " << w << std::endl; std::cin.get();
    
    // construct initial sigma-matrix (formula 29)
    matrix_type sigma = boost::numeric::ublas::zero_matrix<value_type>(6);
    sigma(0,0) = sx*sx;
    sigma(0,5) = sigma(5,0) = invAB*(ex/W+ex/w);
    sigma(1,1) = invAB*(B*ex*W+A*ez*w);
    sigma(1,4) = sigma(4,1) = AB*invAB*(ex*W+ez*w);
    sigma(2,2) = sy*sy;
    sigma(3,3) = ey*ey/sigma(2,2);
    sigma(4,4) = sz*sz;
    sigma(5,5) = std::fabs((ex/(B*W)+ez/(A*w))*invAB);
    
#ifdef PRINT
    out << sigma << std::endl;
#endif
    
    for(typename std::vector<matrix_type>::iterator it=sigmas_.begin(); it!=sigmas_.end(); ++it) {
      *it = sigma;
    }
  }
  
#ifdef PRINT
  out.close();
#endif
}

template<typename Value_type, typename Size_type>
void SigmaGenerator<Value_type, Size_type>::initialize(value_type nuz, value_type ravg) {
  /*
   * The initialization is based on the analytical solution of using a spherical symmetric beam in the paper:
   * Transverse-longitudinal coupling by space charge in cyclotrons
   * by Dr. Christian Baumgarten
   * (formulas: (46), (56), (57)
   */
  
  // helper constants
  value_type invbg = 1.0/(beta_*gamma_);
  value_type micro = 1.0e-6;
  value_type mega = 1.0e6;
  
  /// WHY DON'T WE UPDATE THE EMITTANCE FOR OTHER ENERGIES ???
  emittance_[0] = emittance_[1] = emittance_[2] = 1.254129758715586;
  
  // mm mrad --> meter radian
  value_type ex = emittance_[0]*invbg*micro;	// [ex] = Pi mm mrad
  value_type ey = emittance_[1]*invbg*micro;	// [ey] = Pi mm mrad
  value_type ez = emittance_[2]*invbg*micro;	// [ez] = Pi mm mrad
  
#ifdef PRINT
  std::ofstream out("data/maps/InitialSigmaPerAngle.dat");
#endif
  
  // initial guess of emittance
  value_type e = std::cbrt(ex*ey*ez);	// cbrt computes cubic root (C++11) <cmath>
  
  // normalized emittance
  e *= beta_*gamma_;
  
  // cyclotron radius
  value_type rcyc = ravg/beta_;
  
  // "average" inverse bending radius
  value_type h = 1.0/ravg;
  
  // formula (57)
  value_type lam = 2.0*M_PI*physics::c/(wo_*nh_*micro); // wavelength, wo*nh*1e-6 real orbital frequency in MHz
  value_type K3 = 3.0*physics::q0*I_*lam/(20.0*std::sqrt(5.0)*M_PI*physics::eps0*m_*physics::c*physics::c*physics::c*beta_*beta_*gamma2_*gamma_);
  value_type alpha = physics::q0*physics::mu0*I_/(5.0*std::sqrt(10.0)*m_*physics::c*gamma_*nh_)*std::sqrt(rcyc*rcyc*rcyc/(e*e*e));
  value_type sig0 = std::sqrt(2.0*rcyc*e)/gamma_;
  
  // formula (56)
  value_type sig;
  if(alpha>=2.5) {
    sig = sig0*std::cbrt(1.0+alpha);		// cbrt computes cubic root (C++11) <cmath>
  } else if(alpha>=0) {
    sig = sig0*(1+alpha*(0.25-0.03125*alpha));
  } else {
    Error::message("SigmaGenerator::initialize()",Error::range);
  }
  
  // K = Kx = Ky = Kz
  value_type K = K3*gamma_/(3.0*sig*sig*sig);		// formula (46)
  value_type kx = h*h*gamma2_;			// formula (46) (assumption of an isochronous cyclotron)
  
  
  value_type a = 0.5*kx-K;	// formula (22) (with K = Kx = Kz)
  value_type b = K*K;		// formula (22) (with K = Kx = Kz and kx = h^2*gamma^2)
    
#ifdef DEBUG
  std::cout << "e = " << e << std::endl;
  std::cout << "gamma = " << gamma_ << std::endl;
  std::cout << "beta = " << beta_ << std::endl;
  std::cout << "current = " << I_ << std::endl;
  std::cout << "K3 = " << K3*1e8 << std::endl;
  std::cout << "K = " << K << std::endl;
  std::cout << "kx = " << kx << std::endl;
  std::cout << "a = " << a << std::endl;
  std::cout << "b = " << b << std::endl;
  std::cout << "tmp = " << a*a-b << std::endl;
#endif
  
  // b must be positive, otherwise no real-valued frequency
  if(b<0) {
    PhysicalError::message("SigmaGenerator::initialize()",PhysicalError::negative);
  }
  
  value_type tmp = a*a - b;
  if(tmp<0) {
    PhysicalError::message("SigmaGenerator::initialize()", PhysicalError::negative);
  }
  
  tmp = std::sqrt(tmp);
  
  if(a<tmp) {
    PhysicalError::message("SigmaGenerator::initialize()", PhysicalError::negative);
  }
  
  value_type Omega = std::sqrt(a+tmp);		// formula (22)
  value_type omega = std::sqrt(a-tmp);		// formula (22)
  
  value_type A = h/(Omega*Omega+K);		// formula (26)
  value_type B = h/(omega*omega+K);		// formula (26)
  value_type invAB = 1/(B-A);
    
#ifdef DEBUG
  std::cout << "K = " << K << std::endl
	    << "kx = " << kx << std::endl
	    << "a = " << a << std::endl
	    << "b = " << b << std::endl
	    << "Omega = " << Omega << std::endl
	    << "omega = " << omega << std::endl
	    << "A = " << A << std::endl
	    << "B = " << B << std::endl
	    << "invAB = " << invAB << std::endl; std::cin.get();
#endif
  
  // construct initial sigma-matrix (formula (29, 30, 31)
  matrix_type sigma = boost::numeric::ublas::zero_matrix<value_type>(6);
  sigma(0,0) = invAB*(B*ex/Omega + A*ez/omega)*mega;				// formula (30)
  sigma(0,5) = sigma(5,0) = invAB*(ex/Omega + ez/omega)*mega;
  sigma(1,1) = invAB*(B*ex*Omega+A*ez*omega)*mega;
  sigma(1,4) = sigma(4,1) = invAB*(ex*Omega+ez*omega)/(K*gamma2_)*mega;
  sigma(2,2) = sigma(3,3) = invAB*ey/(std::sqrt(h*h*nuz*nuz-K))*mega;		// formula (31)
  sigma(4,4) = invAB*(A*ex*Omega+B*ez*omega)/(K*gamma2_)*mega;
  sigma(5,5) = invAB*(ex/(B*Omega)+ez/(A*omega))*mega;			// formula (30)
  
#ifdef PRINT
  out << sigma << std::endl << std::endl;
#endif
  
  // fill in initial guess of the sigma matrix (for each angle the same guess)
  sigmas_.resize(/*N_*/nSector_);
  for(typename std::vector<matrix_type>::iterator it=sigmas_.begin(); it!=sigmas_.end(); ++it) {
    *it = sigma;
  }
  
#ifdef PRINT
  out.close();
#endif
}

template<typename Value_type, typename Size_type>
template<class matrix>
void SigmaGenerator<Value_type, Size_type>::reduce(matrix& M) {
  /* The 6x6 matrix gets reduced to a 4x4 matrix in the following way:
   * 
   * a11 a12 a13 a14 a15 a16
   * a21 a22 a23 a24 a25 a26		a11 a12 a15 a16
   * a31 a32 a33 a34 a35 a36	-->	a21 a22 a25 a26
   * a41 a42 a43 a44 a45 a46		a51 a52 a55 a56
   * a51 a52 a53 a54 a55 a56		a61 a62 a65 a66
   * a61 a62 a63 a64 a65 a66
   */
  
  // copy x- and z-direction to a 4x4 matrix_type
  matrix_type M4x4(4,4);
  for(size_type i=0; i<2; ++i) {
      // upper left 2x2 [a11,a12;a21,a22]
      M4x4(i,0) = M(i,0);
      M4x4(i,1) = M(i,1);
      // lower left 2x2 [a51,a52;a61,a62]
      M4x4(i+2,0) = M(i+4,0);
      M4x4(i+2,1) = M(i+4,1);
      // upper right 2x2 [a15,a16;a25,a26]
      M4x4(i,2) = M(i,4);
      M4x4(i,3) = M(i,5);
      // lower right 2x2 [a55,a56;a65,a66]
      M4x4(i+2,2) = M(i+4,4);
      M4x4(i+2,3) = M(i+4,5);
  }
  
  M.resize(4,4,false);
  M.swap(M4x4);
}

template<typename Value_type, typename Size_type>
template<class matrix>
void SigmaGenerator<Value_type, Size_type>::expand(matrix/*_type*/& M) {
  /* The 4x4 matrix gets expanded to a 6x6 matrix in the following way:
   * 
   * 				a11 a12 0 0 a13 a14
   * a11 a12 a13 a14		a21 a22 0 0 a23 a24
   * a21 a22 a23 a24	-->	0   0   1 0 0   0 
   * a31 a32 a33 a34		0   0   0 1 0   0
   * a41 a42 a43 a44		a31 a32 0 0 a33 a34
   * 				a41 a42 0 0 a43 a44
   */
  
  matrix/*_type*/ M6x6 = boost::numeric::ublas::identity_matrix<value_type>(6,6);
  
    for(size_type i=0; i<2; ++i) {
      // upper left 2x2 [a11,a12;a21,a22]
      M6x6(i,0) = M(i,0);
      M6x6(i,1) = M(i,1);
      // lower left 2x2 [a31,a32;a41,a42]
      M6x6(i+4,0) = M(i+2,0);
      M6x6(i+4,1) = M(i+2,1);
      // upper right 2x2 [a13,a14;a23,a24]
      M6x6(i,4) = M(i,2);
      M6x6(i,5) = M(i,3);
      // lower right 2x2 [a22,a34;a43,a44]
      M6x6(i+4,4) = M(i+2,2);
      M6x6(i+4,5) = M(i+2,3);
  }
  
  // exchange
  M.swap(M6x6);
}

template<typename Value_type, typename Size_type>
typename SigmaGenerator<Value_type, Size_type>::matrix_type SigmaGenerator<Value_type, Size_type>::updateInitialSigma(const matrix_type& M, const vector_type& eigen,
							       sparse_matrix_type& R, sparse_matrix_type& invR) {
  
  /* This function is taken from the program of Dr. Christian Baumgarten.
   * File:	cyc.c
   * Function:	MatchedEllipse
   */
  
  // normalize emittances
  value_type invbg = 1.0/(beta_*gamma_);
  value_type ex = emittance_[0]*invbg;
  value_type ey = emittance_[1]*invbg;
  value_type ez = emittance_[2]*invbg;
  
  // new sigma matrix
  matrix_type sigma(6,6);
  
  // compute sinus
  value_type invsinx = 1.0/(matt::sign(eigen(0)+eigen(1))*std::sqrt(std::fabs(-eigen(1)*eigen(0))));
  value_type sinz = matt::sign(eigen(2)+eigen(3))*std::sqrt(std::fabs(-eigen(3)*eigen(2)));
  
#ifdef DEBUG
  std::cout << "eigen = " << eigen << std::endl;
  std::cout << "1/smx = " << invsinx << std::endl;
  std::cout << "smz = " << sinz << std::endl;
#endif
  
  // compute Twiss parameters in x-direction (alphax = 0)
  value_type betax = eigen(0)*invsinx;
  value_type gammax = eigen(1)*invsinx;
  
  // compute Twiss parameters in z-direction
  value_type betaz = 0, gammaz = 0;
  
  if(std::fabs(sinz)>1e-4) {
    betaz = eigen(2)/sinz;
    gammaz = eigen(3)/sinz;
  }
  
  // compute Twiss parameters in y-direction
  value_type alphay = 0, betay = 0, gammay = 0;
  
  value_type cosy = 0.5*(M(2,2)+M(3,3));		// cosinus in y
  value_type siny = std::sqrt(std::fabs(1-cosy*cosy));	// sinus in y
  
  if(M(2,3)<0)
    siny *= -1;
  
  if(std::fabs(siny)>=1e-6) {
    siny = 1./siny;
    alphay = 0.5*siny*(M(2,2)-M(3,3));
    betay = M(2,3)*siny;
    gammay = -M(3,2)*siny;
    
    if(betay<0) {
      alphay *= -1.0;
      betay *= -1.0;
      gammay *= -1.0;
    }
  }
  
#ifdef DEBUG
  std::cout << "ay = " << alphay << std::endl;
  std::cout << "bx = " << betax << " by = " << betay << " bz = " << betaz << std::endl;
  std::cout << "gx = " << gammax << " gy = " << gammay << " gz = " << gammaz << std::endl; std::cin.get();
#endif
  
  // how is this matrix called
  matrix_type D = boost::numeric::ublas::zero_matrix<value_type>(6,6);
  // x-direction
  D(0,1) = betax*ex/*emittance_[0]*/;
  D(1,0) = -gammax*ex/*emittance_[0]*/;
  // y-direction
  D(2,2) = alphay*ey/*emittance_[1]*/;
  D(3,3) = -alphay*ey/*emittance_[1]*/;
  D(2,3) = betay*ey/*emittance_[1]*/;
  D(3,2) = -gammay*ey/*emittance_[1]*/;
  // z-direction
  D(4,5) = betaz*ez/*emittance_[2]*/;
  D(5,4) = -gammaz*ez/*emittance_[2]*/;
  
  // expand 4x4 transformation matrices to 6x6
  expand<sparse_matrix_type>(R);
  expand<sparse_matrix_type>(invR);
  
  // symplectic matrix
  sparse_matrix_type S(6,6,6);
  S(0,1) = S(2,3) = S(4,5) = 1;
  S(1,0) = S(3,2) = S(5,4) = -1;
  
  sigma = matt_boost::gemmm<matrix_type>(invR,D,R);
  sigma = boost::numeric::ublas::prod(sigma,S);
  
  if(sigma(0,0)<0)
    sigma *= -1.0;
  
#ifdef DEBUG
  std::cout << "sigma:" << std::endl << std::setprecision(16) << std::fixed << sigma << std::endl; std::cin.get();
#endif
  return sigma;
}

template<typename Value_type, typename Size_type>
void SigmaGenerator<Value_type, Size_type>::updateSigma(const std::vector<matrix_type>& Mscs, const std::vector<matrix_type>& Mcycs) {
  matrix_type M = boost::numeric::ublas::matrix<value_type>(6,6);
  
  // initial sigma is already computed
  for(size_type i=1; i</*N_*/nSector_; ++i) {
    // transfer matrix for one angle
    M = boost::numeric::ublas::prod(Mscs[i-1],Mcycs[i-1]);
    // transfer the matrix sigma
    sigmas_[i] = matt_boost::gemmm<matrix_type>(M,sigmas_[i-1],boost::numeric::ublas::trans(M));
  }
}

template<typename Value_type, typename Size_type>
typename SigmaGenerator<Value_type, Size_type>::value_type SigmaGenerator<Value_type, Size_type>::L2ErrorNorm(const matrix_type& oldS, const matrix_type& newS) {  
  // compute difference
  matrix_type diff = newS-oldS;
  
  // sum squared error up and take square root
  return std::sqrt(std::inner_product(diff.data().begin(),diff.data().end(),diff.data().begin(),0.0));
}


// template<typename Value_type, typename Size_type>
// typename SigmaGenerator<Value_type, Size_type>::matrix_type SigmaGenerator<Value_type, Size_type>::test(const matrix_type& M, const vector_type& eigen,
// 							       sparse_matrix_type& R, sparse_matrix_type& invR) {
//   
//   /*
//    * This function is based on the paper of Dr. Christian Baumgarten:
//    * Transverse-Longitudinal Coupling by Space Charge in Cyclotrons (2012)
//    */
//   
//   /*
//    * Function input:
//    * - M: one turn transfer matrix
//    * - eigen = {alpha, beta, gamma, delta}
//    * - R: transformation matrix (in paper: E)
//    * - invR: inverse transformation matrix (in paper: E^{-1}
//    */
//   
//   /* formula (18):
//    * sigma = -E*D*E^{-1}*S
//    * with diagonal matrix D (stores eigenvalues of sigma*S (emittances apart from +- i),
//    * skew-symmetric matrix (formula (13)), and tranformation matrices E, E^{-1}
//    */
//   
//   // normalize emittances
//   value_type invbg = 1.0/(beta_*gamma_);
//   value_type ex = emittance_[0]*invbg;
//   value_type ey = emittance_[1]*invbg;
//   value_type ez = emittance_[2]*invbg;
//   
//   // determinants of each block of the transfer matrix
//   value_type invdetx = 1.0/(M(0,0)*M(1,1) - M(0,1)*M(1,0));
//   value_type invdety = 1.0/(M(2,2)*M(3,3) - M(2,3)*M(3,2));
//   value_type invdetz = 1.0/(M(4,4)*M(5,5) - M(4,5)*M(5,4));
//   
//   // normalize blocks by determinants
//   M(0,0) *= invdetx;
//   M(1,1) *= invdetx;
//   M(1,0) *= invdetx;
//   M(0,1) *= invdetx;
//   
//   M(2,2) *= invdety;
//   M(3,3) *= invdety;
//   M(2,3) *= invdety;
//   M(3,2) *= invdety;
//   
//   M(4,4) *= invdetz;
//   M(5,5) *= invdetz;
//   M(4,5) *= invdetz;
//   M(5,4) *= invdetz;
//   
//   
//   
//   
// }


#endif