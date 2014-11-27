#ifndef CLOSEDORBITFINDER_H
#define CLOSEDORBITFINDER_H

#include <array>
#include <cmath>
#include <functional>
#include <numeric>
#include <tuple>
#include <vector>

#include "error.h"
#include "physics.h"
#include "physical_error.h"

#include "MagneticField.h" // ONLY FOR STAND-ALONE PROGRAM


#include <fstream>

// include headers for integration
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>

/// Finds a closed orbit of a cyclotron for a given energy
/*! The algorithm is based on the paper of M. M. Gordon: "Computation of closed orbits and basic focusing properties for
 * sector-focused cyclotrons and the design of 'cyclops'" (1983)
 * As template arguments one chooses the type of the variables and the integrator for the ODEs. The supported steppers can be found on
 * http://www.boost.org/ where it is part of the library Odeint.
 */
template<typename Value_type, typename Size_type, class Stepper>
class ClosedOrbitFinder
{
public:
  /// Type of variables
  typedef Value_type value_type;
  /// Type for specifying sizes
  typedef Size_type size_type;
  /// Type of container for storing quantities (path length, orbit, etc.)
  typedef std::vector<value_type> container_type;
  /// Type for holding state of ODE values
  typedef std::vector<value_type> state_type;
  
  /// Sets the initial values for the integration and calls findOrbit().
  /*!
   * @param E is the (kinetic) energy to which the closed orbit should be found
   * @param wo is the nominal orbital frequency (see paper of Dr. C. Baumgarten: "Transverse-Longitudinal Coupling by Space Charge in Cyclotrons" (2012), formula (1))
   * @param N specifies the number of splits (2pi/N) (default: 360, for each degree a step)
   * @param accuracy specifies the accuracy of the closed orbit (default: 10^-8)
   * @param maxit is the maximal number of iterations done. Program stops if closed orbit not found within this time.
   * @param Emin is the minimum energy [MeV] needed in cyclotron
   * @param Emax is the maximum energy [MeV] reached in cyclotron
   */
  ClosedOrbitFinder(value_type, value_type, size_type, value_type, size_type, value_type, value_type, bool);	// this flag is only temporary
  
  /// Returns the inverse bending radius
  container_type& getInverseBendingRadius();
  
  /// Returns the step lengths of the path
  container_type& getPathLength();
  
  /// Returns the field index
  container_type& getFieldIndex();
  
  /// Returns the radial and vertical tunes (in that order)
  std::tuple<value_type,value_type> getTunes(bool);
  
  /// Returns the closed orbit
  container_type& getOrbit();
  
  /// Returns the relativistic factor gamma
  value_type getGamma();
  
  /// Returns the average orbit radius
  value_type getAverageRadius();
  
  /// Returns the phase
  value_type getPhase();
  
  /// Returns if a closed orbit could be found
  bool isConverged();
  
private:
  /// Computes the closed orbit
  /*!
   * @param accuracy specifies the accuracy of the closed orbit
   * @param maxit is the maximal number of iterations done for finding the closed orbit
   */
  bool findOrbit(value_type, size_type);
  
//   bool findOrbit_new(value_type, size_type);
  
  /// Fills in the values of h_, ds_, fidx_. It gets called by in by constructor.
  void computeOrbitProperties();
  
  /// This function is called by the function getTunes().
  /*! Transfer matrix Y = [y11, y12; y21, y22] (see Gordon paper for more details).
   * @param y are the positions (elements y11 and y12 of Y)
   * @param py2 is the momentum of the second solution (element y22 of Y)
   * @param ncross is the number of sign changes (\#crossings of zero-line)
   */
  value_type computeTune(const std::array<value_type,2>&, value_type, size_type);
  
  value_type christian_computeTune(const std::array<value_type,2>&, value_type, size_type);
  
  /// This function computes nzcross_ which is used to compute the tune in z-direction.
  void computeVerticalOscillations();
  
  /// Stores current position in horizontal direction for the solutions of the ODE with different initial values
  std::array<value_type,2> x_; // x_ = [x1, x2]
  /// Stores current momenta in horizontal direction for the solutions of the ODE with different initial values
  std::array<value_type,2> px_; // px_ = [px1, px2]
  /// Stores current position in longitudinal direction for the solutions of the ODE with different initial values
  std::array<value_type,2> z_; // z_ = [z1, z2]
  /// Stores current momenta in longitudinal direction for the solutions of the ODE with different initial values
  std::array<value_type,2> pz_; // pz_ = [pz1, pz2]
  
  /// Stores the inverse bending radius
  container_type h_;
  /// Stores the step length
  container_type ds_;
  /// Stores the radial orbit
  container_type r_;
  /// Stores the radial momentum
  container_type pr_;
  /// Stores the field index
  container_type fidx_;
  
  /// Counts the number of zero-line crossings in horizontal direction (used for computing horizontal tune)
  size_type nxcross_;
  /// Counts the number of zero-line crossings in vertical direction (used for computing vertical tune)
  size_type nzcross_; //#crossings of zero-line in x- and z-direction
  
  /// Is the energy for which the closed orbit should be found
  value_type E_;
  /// Is the nominal orbital frequency
  value_type wo_;
  /// Is the number of angle splits
  size_type N_;
  /// Is the angle step size
  value_type dtheta_;
  
  /// Is the relativistic factor
  value_type gamma_;
  
  /// Is the average radius
  value_type ravg_;
  
  /// Is the phase
  value_type phase_;
  
  /// Boolean which tells if a closed orbit for this configuration could be found (get set by the function findOrbit)
  bool converged_;
  
  /// Minimum energy needed in cyclotron
  value_type Emin_;
  
  /// Maximum energy reached in cyclotron
  value_type Emax_;
  
  /// Lambda-function
  std::function<value_type(value_type, value_type)> ptheta_;
  
  /// Defines the stepper for integration of the ODE's
  Stepper stepper_;
  
  /// ONLY FOR STAND-ALONE PROGRAM
  float** bmag_;
  
  /// Christian's function for finding equilibrium orbit
  void christian_findOrbit(value_type, size_type);
//   void christian_test(value_type, size_type);
};

// -----------------------------------------------------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type, class Stepper>
ClosedOrbitFinder<Value_type, Size_type, Stepper>::ClosedOrbitFinder(value_type E, value_type wo, size_type N, value_type accuracy,
								     size_type maxit,value_type Emin, value_type Emax, bool flag=true)
: h_(0), ds_(0), r_(N), pr_(N), fidx_(0),
  nxcross_(0), nzcross_(0), E_(E), wo_(wo), N_(N), Emin_(Emin), Emax_(Emax),
  dtheta_(2.0*M_PI/double(N)), gamma_(E/physics::E0+1.0), stepper_()
{
  ptheta_ = [&](value_type p2, value_type x) {
    value_type pts = p2-x*x;
    if(pts<=0) {
      Error::message("ClosedOrbitFinder",Error::invalid); //,"SQUARE ROOT OF NEGATIVE NUMBER");
    }
    return std::sqrt(p2-x*x);
  };
    
  if(Emin_>Emax_ || E_<Emin_ || E>Emax_) {
    Error::message("ClosedOrbitFinder",Error::invalid);
  }
  
  // my functions
  if(flag) {
    // compute closed orbit
    converged_ = findOrbit(accuracy, maxit);
    
    // compute h, ds, fidx, rav (average radius)
    computeOrbitProperties();
  } else {
    // call Christian's functions instead
//     christian_test(accuracy, maxit);
    christian_findOrbit(accuracy, maxit);   
  }
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::container_type& ClosedOrbitFinder<Value_type, Size_type, Stepper>::getInverseBendingRadius() {
  return h_;
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::container_type& ClosedOrbitFinder<Value_type, Size_type, Stepper>::getPathLength() {
  return ds_;
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::container_type& ClosedOrbitFinder<Value_type, Size_type, Stepper>::getFieldIndex() {
  return fidx_;
}

template<typename Value_type, typename Size_type, class Stepper>
std::tuple<Value_type,Value_type> ClosedOrbitFinder<Value_type, Size_type, Stepper>::getTunes(bool flag) {
  // compute radial tune
  value_type nur = computeTune(x_,px_[1],nxcross_);
  // compute nzcross_
  if(flag) {
    computeVerticalOscillations();
  }
  // compute vertical tune
  value_type nuz = computeTune(z_,pz_[1],nzcross_);
    
  return std::make_tuple(nur,nuz);
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::container_type& ClosedOrbitFinder<Value_type, Size_type, Stepper>::getOrbit() {
  return r_;
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::value_type ClosedOrbitFinder<Value_type, Size_type, Stepper>::getGamma() {
  return gamma_;
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::value_type ClosedOrbitFinder<Value_type, Size_type, Stepper>::getAverageRadius() {
  return ravg_;
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::value_type ClosedOrbitFinder<Value_type, Size_type, Stepper>::getPhase() {
  return phase_;
}

template<typename Value_type, typename Size_type, class Stepper>
inline bool ClosedOrbitFinder<Value_type, Size_type, Stepper>::isConverged() {
  return converged_;
} 

// -----------------------------------------------------------------------------------------------------------------------
// PRIVATE MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type, class Stepper>
bool ClosedOrbitFinder<Value_type, Size_type, Stepper>::findOrbit(value_type accuracy, size_type maxit) {
  /* REMARK TO GORDON
   * q' = 1/b = 1/bcon
   * a' = a = acon
   */
  
  // READ IN MAGNETIC FIELD: ONLY FOR STAND-ALONE PROGRAM
  int nsc = 8, nr = 141, Nth = 1440, nth = 1440/8; value_type r0 = 1.8, dr = 0.02;
  bmag_ = MagneticField::malloc2df(Nth,nr);
  MagneticField::ReadSectorMap(bmag_,nr,Nth,1,"data/ring590_bfld.dat",0.0);
  MagneticField::MakeNFoldSymmetric(bmag_,Nth,nr,nth,nsc);
  value_type bint, brint, btint;
  
  // velocity: beta = v/c = sqrt(1-1/(gamma*gamma))
  if(gamma_ == 0) {
    PhysicalError::message("ClosedOrbitFinder::findOrbit",PhysicalError::undefined);
  }
  
  // store acon and bcon locally
  value_type acon = physics::acon(wo_);		// [acon] = m
  value_type invbcon = 1.0/physics::bcon(wo_);	// [bcon] = MeV*s/(C*m^2) = 10^6 T = 10^7 kG (kilo Gauss)
  
  // helper constants
  value_type p2;				// p^2 = p*p
  value_type pr2;				// squared radial momentum (pr^2 = pr*pr)
  value_type ptheta, invptheta;			// Gordon, formula (5c)
  value_type invdenom;				// denominator for computing dr,dpr
  value_type invgamma4/* = 1.0/(gamma2*gamma2)*/;
  value_type xold = 0.0;			// for counting nxcross
  
  // index for reaching next element of the arrays r and pr (no nicer way found yet)
  size_type idx = 0;
  // observer for storing the current value after each ODE step (e.g. Runge-Kutta step) into the containers of r and pr
  auto store = [&](state_type& y, const value_type t)
  {
    r_[idx] = y[0];
    pr_[idx] = y[1];
    
    // count number of crossings (excluding starting point --> idx>0)
    nxcross_ += (idx>0)*(y[4]*xold<0);
    xold = y[4];
    
    ++idx;
  };
  
  // define the six ODEs (using lambda function)
  std::function<void(const state_type&, state_type&, const double)> orbit_integration = [&](const state_type &y, state_type &dydt, const double theta){
    pr2 = y[1]*y[1];
    if(p2 < pr2) {
      PhysicalError::message("ClosedOrbitFinder::findOrbit",PhysicalError::negative);
    }
    // Gordon, formula (5c)
    ptheta = std::sqrt(p2-pr2);
    invptheta = 1.0/ptheta;
    
    // intepolate values of magnetic field
    MagneticField::interpolate(&bint,&brint,&btint,theta*180/M_PI,nr,Nth,y[0],r0,dr,bmag_);
    bint *= invbcon;
    brint *= invbcon;
    
    // Gordon, formula (5a)
    dydt[0] = y[0]*y[1]*invptheta;
    // Gordon, formula (5b)
    dydt[1] = ptheta - y[0]*bint;
    // Gordon, formulas (9a) and (9b)
    for(size_type i=2; i<5; i+=2) {
      dydt[i] = (y[1]*y[i]+y[0]*p2*y[i+1]*invptheta*invptheta)*invptheta;
      dydt[i+1] = -y[1]*y[i+1]*invptheta - (bint+y[0]*brint)*y[i];
    }
    
//     std::cout << std::setprecision(16) << y[2] << " " << y[3] << " " << y[4] << " " << y[5] << std::endl; std::cin.get();
    
  };
  
  // define initial state container for integration: y = {r, pr, x1, px1, x2, px2}
  state_type y(6);
  
  container_type err(2);			// difference of last and first value of r (1. element) and pr (2. element)
  container_type delta = {0.0, 0.0};		// correction term for initial values: r = r + dr, pr = pr + dpr; Gordon, formula (17)
  value_type error;				// amplitude of error; Gordon, formula (18) (a = a')
  size_type niterations = 0;			// if niterations > maxit --> stop iteration
  
  /*
   * Christian:
   * N = 1440 ---> N = 720 ---> dtheta = 2PI/720 --> N = 721
   * 
   * 0, 2, 4, ... ---> jeden zweiten berechnene: 1, 3, 5, ... interpolieren --> 1440 Werte
   * 
   * Matthias:
   * N = 1440 --> dtheta = 2PI/1440 --> N = 1441
   * 
   * 0, 1, 2, 3, 4, 5, ... --> 1440 Werte
   * 
   */
  N_ += 1;				/// ---> that was the error: instead of N, do N+1 steps
  
  r_.resize(N_);
  pr_.resize(N_);
  
  // iterate until suggested energy (start with minimum energy)
  value_type E = Emin_;
  
  // energy dependent values
  value_type en = E/physics::E0;		// en = E/E0 = E/(mc^2) (E0 is kinetic energy)
  value_type p = acon*std::sqrt(en*(2.0+en));	// momentum [p] = m; Gordon, formula (3)
  value_type gamma2 = (1.0+en)*(1.0+en);
  value_type beta = std::sqrt(1.0-1.0/gamma2);
  p2 = p*p;					// p^2 = p*p
  invgamma4 = 1.0/(gamma2*gamma2);
  
  // set initial values for radius and radial momentum for lowest energy Emin
  // orbit, [r] = m; Gordon, formula (20)
  // radial momentum; Gordon, formula (20)
  container_type init = {beta*acon, 0.0};
  
  // store initial values for updating values for higher energies
  container_type previous_init = {0.0, 0.0};
  
  // difference in initial values (current and previous)
//   container_type diff = {0.0, 0.0};
  
  while(E <= E_) {
    
    // (re-)set inital values for r and pr
    r_[0] = init[0]; 
    pr_[0] = init[1];
    
    // integrate until error smaller than user-define accuracy
    do {
      // (re-)set inital values
      x_[0] = 1.0;		// x1; Gordon, formula (10)
      px_[0] = 0.0;		// px1; Gordon, formula (10)
      x_[1] = 0.0;		// x2; Gordon, formula (10)
      px_[1] = 1.0;		// px2; Gordon, formula (10)
      nxcross_ = 0;		// counts the number of crossings of x-axis (excluding first step)
      idx = 0;			// index for looping over r and pr arrays
      
      // fill container with initial states
      y = {/*r_[0]+delta[0]*/init[0],init[1]/*pr_[0]+delta[1]*/, x_[0], px_[0], x_[1], px_[1]};
      
//       std::cerr << std::setprecision(16) << "init: " << y[0] << " " << y[1] << std::endl;
      
      // integrate: assume no imperfections --> only integrate over a single sector (dtheta_ = 2pi/N)
      boost::numeric::odeint::integrate_n_steps(stepper_,orbit_integration,y,0.0,dtheta_,N_,store);
      
      // write new state
      x_[0] = y[2];
      px_[0] = y[3];
      x_[1] = y[4];
      px_[1] = y[5];
      
      // compute error
      err[0] = r_[N_-1] - r_[0];	// Gordon, formula (14)
      err[1] = pr_[N_-1] - pr_[0];	// Gordon, formula (14)
      
//       std::cerr << std::setprecision(16) << "last: " << r_[N_-1] << " " << pr_[N_-1] << std::endl; std::cin.get();
      
      // correct inital values of r and pr
      invdenom = 1.0/(x_[0]+px_[1]-2.0);
      delta[0] = ((px_[1]-1.0)*err[0] - x_[1]*err[1])*invdenom;	// dr; Gordon, formula (16a)
      delta[1] = ((x_[0]-1.0)*err[1]-px_[0]*err[0])*invdenom;	// dpr; Gordon, formula (16b)
      
      // improved initial values; Gordon, formula (17) (here it's used for higher energies)
      init[0] += delta[0];
      init[1] += delta[1];
      
      // compute amplitude of the error
      error = std::sqrt(delta[0]*delta[0]+delta[1]*delta[1]*invgamma4)/r_[0];
      
    } while(error>accuracy && niterations++<maxit);
    
    // reset iteration counter
    niterations = 0;
    
    // reset correction term
    delta[0] = delta[1] = 0.0;
    
    // increase energy by one
    E += 1.0;
    
    // set constants for new energy E
    en = E/physics::E0;			// en = E/E0 = E/(mc^2) (E0 is kinetic energy)
    p = acon*std::sqrt(en*(2.0+en));	// momentum [p] = m; Gordon, formula (3)
    p2 = p*p;				// p^2 = p*p
    gamma2 = (1.0+en)*(1.0+en);
    invgamma4 = 1.0/(gamma2*gamma2);
  }
  
  // remove last entry
  r_.pop_back();
  pr_.pop_back();
  N_ -= 1;
  
  // returns true if converged, otherwise false
  return error<accuracy;
}

// template<typename Value_type, typename Size_type, class Stepper>
// bool ClosedOrbitFinder<Value_type, Size_type, Stepper>::findOrbit_new(value_type accuracy, size_type maxit) {
//   /* REMARK TO GORDON
//    * q' = 1/b = 1/bcon
//    * a' = a = acon
//    */
//   
//   // READ IN MAGNETIC FIELD: ONLY FOR STAND-ALONE PROGRAM
//   int nsc = 8, nr = 141, Nth = 1440, nth = 1440/8; value_type r0 = 1.8, dr = 0.02;
//   bmag_ = MagneticField::malloc2df(Nth,nr);
//   MagneticField::ReadSectorMap(bmag_,nr,Nth,1,"data/ring590_bfld.dat",0.0);
//   MagneticField::MakeNFoldSymmetric(bmag_,Nth,nr,nth,nsc);
//   value_type bint, brint, btint;
//   
//   // velocity: beta = v/c = sqrt(1-1/(gamma*gamma))
//   if(gamma_ == 0) {
//     PhysicalError::message("ClosedOrbitFinder::findOrbit",PhysicalError::undefined);
//   }
//   value_type gamma2 = gamma_*gamma_;
//   value_type beta = std::sqrt(1.0-1.0/gamma2);
//   
//   // acon, bcon
//   value_type acon = physics::c/wo_;				// [acon] = m
//   value_type bcon = physics::E0*1.0e7/(physics::q0*physics::c*acon);	// [bcon] = MeV*s/(C*m^2) = 10^6 T = 10^7 kG (kilo Gauss)
//   /*
//    * christian's bcon = 5.53726858635484
//    */
//   
//   // acon = 5.65213418581506 bcon = 5.53726858635484
// //   std::cout << std::setprecision(15) << "acon = " << acon << " bcon = " << bcon << std::endl; std::cin.get();
//   
//   // helper constants
//   value_type p = acon*std::sqrt(E_/physics::E0*(2.0+E_/physics::E0));	// momentum [p] = m; Gordon, formula (3)
//   value_type p2 = p*p;
//   value_type pr2;							// squared radial momentum
//   value_type ptheta, invptheta;
//   value_type invdenom;							// denominator for computing dr,dpr
//   value_type invgamma4 = 1.0/(gamma2*gamma2);
//   value_type xold;				// for counting nxcross
//   
//   // set initial values for radius and radial momentum
//   r_[0] = /*beta*acon*/ p;	// orbit, [r] = m; Gordon, formula (20)
//   pr_[0] = 0.0;		// radial momentum; Gordon, formula (20)
//   
//   // define initial state container for integration: y = {r, pr, x1, px1, x2, px2}
//   state_type y(6);
//   
//   container_type err(2);
//   container_type delta = {0.0, 0.0};		// correction term for initial values: r = r + dr, pr = pr + dpr; Gordon, formula (17)
//   value_type error;				// amplitude of error; Gordon, formula (18) (a = a')
//   
//   size_type niterations = 0;			// if niterations > maxit --> stop iteration
//   
//   
//   double ark[4]={0.5, 0.5, 1.0};
//   double brk[4]={1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
//   double crk[4]={0, 0.5, 0.5, 1};
//   double k[4*6];
//   
//   N_ *= 0.5;
//   N_ += 1;				/// ---> that was the error: instead of 720, do 721 steps
//   std::cerr << N_ << std::endl;
//   
//   r_.resize(N_);
//   pr_.resize(N_);
//   
//   std::cerr << 2.0*dtheta_ << std::endl;
//   // integrate until error smaller than user-define accuracy
//   do {
//     // (re-)set inital values
//     x_[0] = 1.0;		// x1; Gordon, formula (10)
//     px_[0] = 0.0;		// px1; Gordon, formula (10)
//     x_[1] = 0.0;		// x2; Gordon, formula (10)
//     px_[1] = 1.0;		// px2; Gordon, formula (10)
//     nxcross_ = 0;
//     
//     
//     
//     // fill container with initial states
//     y = {r_[0]+delta[0],pr_[0]+delta[1], x_[0], px_[0], x_[1], px_[1]};
//     
//     std::cerr << std::setprecision(16) << "init: " << y[0] << " " << y[1] << std::endl; //std::cin.get();
//     
//     // integrate: assume no imperfections --> only integrate over a single sector (dtheta_ = 2pi/N)
//     state_type yold(6);
//     state_type dydt(6);
//     
//     
//     for(int i=0; i<N_; ++i) {
//       
//       r_[i] = y[0];
//       pr_[i] = y[1];
//       
//       for(int n=0; n<6 ; ++n) {
// 	yold[n] = y[n];
//       }
//       
//       for(int j=0; j<4; ++j) {
// 	
// 	pr2 = y[1]*y[1];
// 	if(p2 < pr2) {
// 	  PhysicalError::message("ClosedOrbitFinder::findOrbit",PhysicalError::negative);
// 	}
// 	
// 	// Gordon, formula (5c)
// 	ptheta = std::sqrt(p2-pr2);
// 	invptheta = 1.0/ptheta;
// 	
// 	int ith = (j+1)/2 + 2*i;
// 	MagneticField::interpolate1(&bint,&brint,&btint,ith,nr,Nth,y[0],r0,dr,bmag_,false);
// 	
// 	bint /= bcon;
// 	brint /= bcon;
// 	
// // 	std::cout << std::setprecision(15) << "bint = " << bint << " brint = " << brint << std::endl; std::cin.get(); 
// 	
// 	// Gordon, formula (5a)
// 	dydt[0] = y[0]*y[1]*invptheta;
// 	// Gordon, formula (5b)
// 	dydt[1] = ptheta - y[0]*bint;
// 	// Gordon, formulas (9a) and (9b)
// 	for(size_type l=2; l<5; l+=2) {
// 	  dydt[l] = (y[1]*y[l]+y[0]*p2*y[l+1]*invptheta*invptheta)*invptheta;
// 	  dydt[l+1] = -y[1]*y[l+1]*invptheta - (bint+y[0]*brint)*y[l];
// 	}
// 	
// 	for (int n=0;n<6;n++) {
// 	  k[j+4*n] = dydt[n];
// 	  y[n] = yold[n] + 2*dtheta_*ark[j]*k[j+4*n];
// 	}
//       }
//       
//       for(int n=0; n<6; ++n) {
// 	y[n] = yold[n] + 2*dtheta_*(brk[0]*k[0+4*n] + brk[1]*k[1+4*n] + brk[2]*k[2+4*n] + brk[3]*k[3+4*n]);
//       }
//     }
//     
//     
//     // write new state
//     x_[0] = y[2];
//     px_[0] = y[3];
//     x_[1] = y[4];
//     px_[1] = y[5];
//     
//     // compute error
//     err[0] = r_[N_-1] - r_[0];		// Gordon, formula (14)
//     err[1] = pr_[N_-1] - pr_[0];	// Gordon, formula (14)
//     
//     std::cerr << "r[0] = " << r_[0] << " pr[0] = " << pr_[0] << std::endl; std::cin.get();
//     std::cerr << "r[N-1] = " << r_[N_-1] << " pr[N-1] = " << pr_[N_-1] << std::endl;
//     
//     // correct inital values of r and pr
//     invdenom = 1.0/(x_[0]+px_[1]-2.0);
//     delta[0] = ((px_[1]-1.0)*err[0] - x_[1]*err[1])*invdenom;	// dr; Gordon, formula (16a)
//     delta[1] = ((x_[0]-1.0)*err[1]-px_[0]*err[0])*invdenom;	// dpr; Gordon, formula (16b)
//     
//     // compute amplitude of the error
//     error = std::sqrt(delta[0]*delta[0]+delta[1]*delta[1]*invgamma4)/r_[0];
//     std::cerr << "error = " << error << std::endl;
//     
//   } while(error>accuracy && niterations<maxit);
//   
//   // returns true if converged, otherwise false
//   return error<accuracy;
// }

template<typename Value_type, typename Size_type, class Stepper>
Value_type ClosedOrbitFinder<Value_type, Size_type, Stepper>::computeTune(const std::array<value_type,2>& y, value_type py2, size_type ncross) {
  // Y = [y1, y2; py1, py2]
  
  // cos(sig)
  value_type cos = 0.5*(y[0]+py2);
  // sig' (see Gordon)
  value_type sig_prime;
  value_type twopi = 2.0*M_PI;
  // store the number of crossings
  value_type n = ncross;
  
  if(std::fabs(cos)>1) {
    // sig complex --> sig' purely imaginary
    
    // if sig complex -->n = n' ord n = n' + 1
    // if n (= ncross) uneven --> n' = n - 1
    n -= (ncross%2)*1.0;
    // Gordon, formula (36b)
    sig_prime = -std::acosh(std::fabs(cos));
  } else {
    // sig real --> 0 <= sig' <= pi
    if(ncross%2) {
      // if negative --> acos(-cos), else -acos(-cos); Gordon, formula (36a)
      // sign of sinus has to match sign of y12
      sig_prime = (std::signbit(y[1])) ? std::acos(-cos) : -std::acos(-cos);
    } else {
      // if negative --> -acos(cos), else acos(cos); Gordon, formula (36b)
      // sign of sinus has to match sign of y12
      sig_prime = (std::signbit(y[1])) ? -std::acos(cos) : std::acos(cos); 
      /* shift sig' into [0,2*pi]
       * (didn't match the result of Christian, when using the proposed interval of [0,pi] by Gordon)
       */
      sig_prime += twopi*std::abs(std::floor(sig_prime/twopi));
    }
  }
  // Gordon, formula (36)
  /*
   * sig complex --> sig = n*pi + sig'
   * sig real	 --> sig = n*pi + |sig'|
   */
  double sig = n*M_PI + sig_prime;
  
  // \nu = sig/theta, where theta = integration domain
  return sig/(2.0*M_PI);
}

template<typename Value_type, typename Size_type, class Stepper>
void ClosedOrbitFinder<Value_type, Size_type, Stepper>::computeOrbitProperties() {
  /* 
   * The formulas for h, fidx and ds are from the paper:
   * "Tranverse-Longitudinal Coupling by Space Charge in Cyclotrons"
   * written by Dr. Christian Baumgarten (2012, PSI)
   * p. 6
   */
  
  
  // READ IN MAGNETIC FIELD: ONLY FOR STAND-ALONE PROGRAM
  int nsc = 8, nr = 141, Nth = 1440, nth = 1440/8; value_type r0 = 1.8, dr = 0.02;
  value_type bint, brint, btint; // B, dB/dr, dB/dtheta
  
  value_type invbcon = 1.0/physics::bcon(wo_);
  value_type en = E_/physics::E0;				// en = E/E0 = E/(mc^2) (E0 is kinetic energy)
  value_type p = physics::acon(wo_)*std::sqrt(en*(2.0+en));	// momentum [p] = m; Gordon, formula (3)
  value_type p2 = p*p;
  value_type theta = 0.0;					// angle for interpolating
  value_type ptheta, drdtheta;
  
  // resize of container
  h_.resize(N_);
  fidx_.resize(N_);
  ds_.resize(N_);
  
  for(size_type i=0; i<N_; ++i) {
    // interpolate magnetic field
    MagneticField::interpolate(&bint,&brint,&btint,theta*180.0/M_PI,nr,nth*nsc,r_[i],r0,dr,bmag_);
    bint *= invbcon;
    brint *= invbcon;
    btint *= invbcon;
    
    // inverse bending radius (B0/bcon = 1 ?)
    h_[i] = bint/p;
    
    // local field index
    ptheta = std::sqrt(p2-pr_[i]*pr_[i]);
    fidx_[i] = (brint*ptheta - btint*pr_[i]/r_[i])/p2; //(bint*bint);
    
    // path length element
    ds_[i] = std::hypot(r_[i]*pr_[i]/ptheta,r_[i])*dtheta_;	// C++11 function
    
    // increase angle
    theta += dtheta_;
  }
  
  // compute average radius
  ravg_ = std::accumulate(r_.begin(),r_.end(),0.0)/value_type(r_.size());
  
}

template<typename Value_type, typename Size_type, class Stepper>
void ClosedOrbitFinder<Value_type, Size_type, Stepper>::computeVerticalOscillations() {
  
  // READ IN MAGNETIC FIELD: ONLY FOR STAND-ALONE PROGRAM
  int nsc = 8, nr = 141, Nth = 1440, nth = 1440/8; value_type r0 = 1.8, dr = 0.02;
  value_type bint, brint, btint; // B, dB/dr, dB/dtheta
  
  value_type en = E_/physics::E0;				// en = E/E0 = E/(mc^2) with kinetic energy E0
  value_type p = physics::acon(wo_)*std::sqrt(en*(en+2.0));	// Gordon, formula (3)
  value_type  p2 = p*p;						// p^2 = p*p
  size_type idx = 0;						// index for going through container
  value_type pr2;						// pr^2 = pr*pr
  value_type ptheta, invptheta;					// Gordon, formula (5c)
  value_type zold = 0.0;					// for counting nzcross
  
  // store bcon locally
  value_type invbcon = 1.0/physics::bcon(wo_);		// [bcon] = MeV*s/(C*m^2) = 10^6 T = 10^7 kG (kilo Gauss)
  
  // define the ODEs (using lambda function)
  std::function<void(const state_type&, state_type&, const double)> vertical = [&](const state_type &y, state_type &dydt, const double theta){
    pr2 = y[1]*y[1];
    if(p2 < pr2) {
      PhysicalError::message("ClosedOrbitFinder::findOrbit",PhysicalError::negative);
    }
    // Gordon, formula (5c)
    ptheta = std::sqrt(p2-pr2);
    invptheta = 1.0/ptheta;
    
    // We have to integrate r and pr again, otherwise we don't have the Runge-Kutta of the B-field
    // Gordon, formula (5a)
    dydt[0] = y[0]*y[1]*invptheta;
    // Gordon, formula (5b)
    dydt[1] = ptheta - y[0]*bint;
    
    // interpolate magnetic field
    MagneticField::interpolate(&bint,&brint,&btint,theta*180.0/M_PI,nr,nth*nsc,y[0],r0,dr,bmag_);
    bint *= invbcon;	
    brint *= invbcon;
    btint *= invbcon;
    
    // Gordon, formulas (22a) and (22b)
    for(size_type i=2; i<5; i+=2) {
      dydt[i] = y[0]*y[i+1]*invptheta;
      dydt[i+1] = (y[0]*brint - y[1]*invptheta*btint)*y[i];
    }
  };
  
  // to get next index for r and pr (to iterate over container)
  auto next = [&](state_type& y, const value_type t) {
    // number of times z2 changes sign
    nzcross_ += (idx>1)*(y[4]*zold<0);
    zold = y[4];
    ++idx;
  };
  
  // set initial state container for integration: y = {r, pr, z1, pz1, z2, pz2}
  state_type y = {r_[0], pr_[0], 1.0, 0.0, 0.0, 1.0};
  
  
  // integrate: assume no imperfections --> only integrate over a single sector (dtheta_ = 2pi/N)
  boost::numeric::odeint::integrate_n_steps(stepper_,vertical,y,0.0,dtheta_,N_,next);
    
  // write new state
  z_[0] = y[2];
  pz_[0] = y[3];
  z_[1] = y[4];
  pz_[1] = y[5];
}

// ------------------------------------------------------------------------------------------------------------------------
//					CHRISTIAN'S FUNCTIONS
// ------------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type, class Stepper>
Value_type ClosedOrbitFinder<Value_type, Size_type, Stepper>::christian_computeTune(const std::array<value_type,2>& y, value_type py2, size_type ncross) { // Y = [y1, y2; py1, py2]
  // compute cosinus and sinus
  value_type cos = 0.5*(y[0]+py2);	// formula (33c)
  value_type sin2 = 1.0-cos*cos;
  
  value_type two_pi = 2.0*M_PI;
  
  value_type tmp;
  
  if(sin2<0) {	// if abs(cos) > 1 --> motion is unstable
    tmp = -std::log(std::fabs(cos)+std::sqrt(-sin2));
  } else {
    sin2 = std::sqrt(sin2);
    if(y[1]<0) {
      sin2 *= -1;
    }
    tmp = std::atan2(sin2,cos);
    if(tmp<0) {
      tmp += two_pi;
    }
  }
  tmp /= two_pi;
  
  return tmp + std::floor(0.5*ncross);
}

template<typename Value_type, typename Size_type, class Stepper>
void ClosedOrbitFinder<Value_type, Size_type, Stepper>::christian_findOrbit(value_type accuracy, size_type maxit) {
  
  // READ IN MAGNETIC FIELD: ONLY FOR STAND-ALONE PROGRAM
  int nsc = 8, nr = 141, Nth = 1440, nth = 1440/8; value_type r0 = 1.8, dr = 0.02;
  bmag_ = MagneticField::malloc2df(Nth,nr);
  MagneticField::ReadSectorMap(bmag_,nr,Nth,1,"data/ring590_bfld.dat",0.0);
  MagneticField::MakeNFoldSymmetric(bmag_,Nth,nr,nth,nsc);
  
  // READ IN MAGNETIC FIELD: ONLY FOR STAND-ALONE PROGRAM
  value_type bint, brint, btint; // B, dB/dr, dB/dtheta
  
  std::cout << "Christian's findOrbit" << std::endl;
  
  double ark[4]={0.5,0.292893219,1.707106781,0.1666666667};
  double brk[4]={2.0,1.0,1.0,2.0};
  double crk[4]={-0.5,-0.292893219,-1.707106781,-0.5};
  
  h_.resize(N_);
  ds_.resize(N_);
  fidx_.resize(N_);
  
  double stp,psq,pp,gam,pt,pts,tpi,rg,prg;
    double xc0,xc1,xc2,xc3,xc4,ep1,ep2,den,yp5,yp10;
    double rp,prp,dth;
    // Start values of R,Pr:
    double rrk,y[12],ck[12],q[12];
    double sect,de;
    double phs;
    double e,rav,rmax,rmin;
    double cn,sn;
    double rznt,rznu[2];
//     double bint,brint,btint;
    double tmp;
    FILE *f;
    char fname[200];
    int nstp,ntry,neq,nxcros,nzcros,eqdst;
    int i,j,n0;
    int ir0,ith,istp;
    int done=0;
    /***/
    sprintf(fname,"eo_params.dat");
    f=fopen(fname,"w");
    /***/
    de=1.0;
    eqdst=(de>0.0);
    tpi=2.0*M_PI;
    nth=nsc*nth;			// nsc = #sectors
    dth=tpi/(double)nth;			// dth = angle step (d\theta)
    sect=tpi;//(double)(cyc->nsc);
    nstp=nth/2;
//     n0=(int)((cyc->emax-cyc->emin)/de+1.5);
    n0 = 0;
    // Allocate Memory:
    stp=tpi/(double)(nstp);
    
    rp=0.0;
    prp=0.0;
    double nuz=0.0;
    double nur=0.0;
    phs=0.0;
//     ir0=(int)(1.5-cyc->r0/cyc->dr);
    //dpeo[0]=phs;
    //==============================
    // PSI RING:
    eqdst=1;
    double emin = 72.0;
    e=emin;
    int neo = 0;
    while(e < E_+de)
    {
    //==============================
    // Loop over all energies:
//     do {
	// Note Energy in "de":
// 	if (eqdst) {
// 	  //if (fabs(e-0.87)<0.1) e=1.0; else
// 	  e+=de;
// 	} else {
// 	    de=e;
// 	    // Calculate new energy:
// 	    tmp=SQR((double)(cyc->neo+1)/(double)(n0-1));
// 	    e=cyc->emax * 3.0*tmp/(1.+2.0*tmp);
// 	    // Calculate new "de":
// 	    de=e-de;
// 	}
	tmp=e/physics::E0;
	// Relativistic gamma-factor = (E+E0)/E0 
        // because E0 is rest mass energy and E is kinetic energy 
	gam=1.0+tmp;
	// psq = squared total momentum 
	// 
	// P^2 = A^2 * (2 E/E0+(E/E0)^2)
	//
	// A="acon" is "cyclotron unit". 
	psq=tmp*(2.0+tmp)*physics::acon(wo_)*physics::acon(wo_);			// Gordon: formula (3) squared
	pp=std::sqrt(psq);						// Gordon: formula (3)
	// "ntry" is the iteration counter.
	//  ntry > 20 : No success.
	//  ntry < 0  : Last iteration (use ALL equations of motion)
	ntry=0;
	done=0;
// 	cyc->E[cyc->neo]=e;
// 	cyc->PC[cyc->neo]=cyc->e0*std::sqrt(tmp*(2.0+tmp));
// 	cyc->gam[cyc->neo]=gam;
// 	cyc->gam2[cyc->neo]=SQR(gam);
	if (neo<=1) {
	    rg=std::sqrt(psq);
	    prg=prp;
	}
	// Iterate until conversion:
	do {
	    // Set initial values for the integration variables:
	    // Radius:
	  y[0]=rg;
	  // Radial Momentum:
	  y[1]=prg;
//     The Equations 3..6 describe the x-motion.
//     There are two identical pairs of equations for
//     two (orthogonal) pairs of starting conditions
	    y[2]=1.0;
	    y[3]=0.0;
	    y[4]=0.0;
	    y[5]=1.0;
//     Equation 7 is the "Phase-Equation"...
	    y[6]=0.0;
//     The Equations 8..11 describe the z-motion.
//     There are two identical pairs of equations for
//     two (orthogonal) pairs of starting conditions
	    y[7]=1.0;
	    y[8]=0.0;
	    y[9]=0.0;
	    y[10]=1.0;
//     Equation 12 the Equation for Rav:
	    y[11]=0.0;
	    for (i=0;i<12;i++) q[i]=0.0;
	    istp=0;
	    nxcros=nzcros=0;
	    // Runge-Kutta Integration over Theta:
	    do {
		r_[2*istp]=y[0];
		// Remember last values y[4] and y[9]:
		yp5=y[4];
		yp10=y[9];
		for (j=0;j<4;j++) { 
		    ith=((j+1)/2+2*istp);
		    MagneticField::interpolate1(&bint,&brint,&btint,ith,nr,nth,y[0],r0,dr,bmag_,ntry<0);
		    // "bint" = interpolated magn. field value
		    bint/=physics::bcon(wo_);
		    // "brint" = radial derivative of magn. field.
		    brint/=physics::bcon(wo_);
                    // "btint" = theta derivative (last time only);
		    if (ntry<0) btint/=physics::bcon(wo_);
		    // "pts" = Ptheta^2
		    pts=psq-y[1]*y[1];
		    if (pts<=0.0) { 
		      std::cerr << "Error in ClosedOrbitFinder." << std::endl;
			return;
		    }
		    pt=std::sqrt(pts);
		    xc0=y[0]/pt;
	            // R - Equation:
		    // dR = dth * r/pt*pr;
		    ck[0]=xc0*y[1];
	            // Pr - Equation:
                    // dPR = dth*(PT-Q*R*Bz)
		    // For the correct units, Pr is given in cyclotron units
		    // and Q*Bz in units of Bcon.
		    ck[1]=pt-y[0]*bint;
		    if ((ntry<0)&&(j==0)) {
		      // dR/dth:
// 		      r_[2*istp]=ck[0];
		      // d2R/dth^2:
		      // Inverse Bending-Radius:
		      h_[2*istp]=bint/pp;
		      // Feldindex:
		      fidx_[2*istp]=(brint*pt-btint*y[1]/y[0])/psq;
		      ds_[2*istp]=std::sqrt(y[0]*y[0]+ck[0]*ck[0])*dth;
		    }
		    //     The Equations 3..6 describe the x-motion.
		    //     There are two identical pairs of equations for
		    //     two (orthogonal) pairs of starting conditions
		    //     Refer to M.M. Gordon in:
		    //     Particle Accelerators Vol. 16, pp 39-62 (1984)
		    xc1=y[1]/pt;
		    xc2=psq*xc0/pts;
		    xc3=y[0]*brint+bint;
		    // dx = dth ( (pr/pt) * x + R * p^2/pt^3 * px)
		    ck[2]= xc1*y[2] + xc2*y[3];
		    // dpx = - dth ( (pr/pt) * px + (r*dB/dr+B) * x)
		    ck[3]=-xc1*y[3] - xc3*y[2];
		    ck[4]= xc1*y[4] + xc2*y[5];
		    ck[5]=-xc1*y[5] - xc3*y[4];
// 		    std::cout << std::setprecision(16) << ck[1] << std::endl; std::cin.get();
		    neq=6;
		    if (ntry<0) {
                        //     phase and z equations (last time only);
			neq=12;
                        //     Equation 7 is the "Phase-Equation" 
			// dphi=dth*(r/pt*gamma-1)
			ck[6]=xc0*gam-1.0;
                        //     The Equations 8..11 describe the z-motion.
                        //     There are two identical pairs of equations for
                        //     two (orthogonal) pairs of starting conditions
			//  
			// dz=dth*(r/pt*pz)
			xc4=y[0]*brint-xc1*btint;
			ck[7]=xc0*y[8];
			// dpz=dth*(r*dB/dr-pr/pt*dB/dth)*z;
			ck[8]=xc4*y[7];
			// dz=dth*(r/pt*pz)
			ck[9]=xc0*y[10];
			// dpz=dth*(r*dB/dr-pr/pt*dB/dth)*z;
			ck[10]=xc4*y[9];
                        //     Equation 12 is the Equation for Rav
			// drav=dth * r
			ck[11]=y[0];
		    }
		    // Normierung auf Schrittweite:
		    for (i=0;i<neq;i++) ck[i]*=stp;
		    // Runge-Kutta-Integration:
		    for (i=0;i<neq;i++) { 
			rrk=ark[j]*(ck[i]-brk[j]*q[i]);
			y[i]+=rrk;
			q[i]+=3.0 * rrk+crk[j]*ck[i];
		    }
		}// for (j=0;j<4;j++)
		istp++;
		if (istp>1) {
		    if (yp5*y[4]<0.0) nxcros++;
		    if (yp10*y[9]<0.0) nzcros++;
		}
	    } while (istp<nstp);
	    if (ntry>=0) {
		// Compare to M.M. Gordon: "Computation of closed orbits and basic
		// focusing properties for sector-focused cyclotrons and the design
		// of cyclops", Part. Acc. Vol 16 (1984), pp. 39-62. Here page 45:
		ep1=y[0]-rg;  // "eps1"
		ep2=y[1]-prg; // "eps2"
                // Denominator (determinant von X(thf,thi)-1):
		den=(1.0-y[2])*(1.0-y[5])-y[3]*y[4];
		rg+=(ep1*(1.0-y[5])+ep2*y[4])/den;
		prg+=(ep1*y[3]+ep2*(1.0-y[2]))/den;
		if (std::fabs(ep1)+std::fabs(ep2)>=1.0E-8*rg) {
		    ntry++;
		} else ntry=-1;
	    } else done = 1;
	} while (!done);
	phs=y[6]/sect;
	rav=y[11]/sect;
	ravg_=y[11]/sect;
	phase_=y[6]/sect;
// 	rmax=0.0;
// 	rmin=1.0E20;
	for (i=0;i<nstp;i++) {
	  ith=i*2+1;
	  r_[ith % nth]=0.5*(r_[(ith-1) % nth]+r_[(ith+1) % nth]);
// 	  cyc->r1eo[cyc->neo][ith % nth]=0.5*(cyc->r1eo[cyc->neo][(ith-1) % nth]+cyc->r1eo[cyc->neo][(ith+1) % nth]);
	  h_[ith % nth]=0.5*(h_[(ith-1) % nth]+h_[(ith+1) % nth]);
	  fidx_[ith % nth]=0.5*(fidx_[(ith-1) % nth]+fidx_[(ith+1) % nth]);
	  ds_[ith % nth]=0.5*(ds_[(ith-1) % nth]+ds_[(ith+1) % nth]);
// 	  if (cyc->reo[cyc->neo][2*i]>rmax) rmax=cyc->reo[cyc->neo][2*i];
// 	  if (cyc->reo[cyc->neo][2*i]<rmin) rmin=cyc->reo[cyc->neo][2*i];
	}
	//tmp=0.0;
	//for (i=0;i<nth;i++) tmp+=cyc->dss[cyc->neo][i];
// 	cyc->Rmax[cyc->neo]=rmax;
// 	cyc->Rmin[cyc->neo]=rmin;
	//fprintf(stderr,"2 pi Rmax=%g, S=%g, 2 pi Rmin=%g",tpi*rmin,tmp,tpi*rmax);
// 	for (i=0;i<2;i++) { 
// 	    cn=0.5*(y[2+5*i]+y[5+5*i]);
// 	    sn=1.-cn*cn;
// 	    if (sn<0) { 
// 		rznt=-std::log(std::fabs(cn)+std::sqrt(-sn));
// 	    } else {
// 		sn=std::sqrt(sn);
// 		if (y[4+5*i]<0) sn=-sn;
// 		rznt=std::atan2(sn,cn);
// 		if (rznt<0) rznt+=tpi;
// 	    }
// 	    rznu[i]=rznt/sect;
// 	}
// 	// Eintragen der Radius-Daten:
//     
// 	//dpeo[cyc->neo]=phs;
// 	while (nxcros>1) { rznu[0]+=1.0; nxcros-=2; }	// == rznu[0] += std::floor(nxcros/2)
// 	while (nzcros>1) { rznu[1]+=1.0; nzcros-=2; }	// == rznu[1] += std::floor(nrcros/2)
// 	nuz=rznu[1];//+(double)(nzcros/2);
// 	nur=rznu[0];//+(double)(nxcros/2);
	neo++;
// 	if (cyc->neo<n0) {
	    ep1=rg-rp;
	    ep2=prg-prp;
	    rp=rg;
	    rg+=ep1;
	    prp=prg;
	    prg+=ep2;
// 	}
//     } while (cyc->neo<n0);
    x_[0] = y[2];
    px_[0] = y[3];
    
    x_[1] = y[4];
    px_[1] = y[5];
    
    z_[0] = y[7];
    pz_[0] = y[8];
    z_[1] = y[9];
    pz_[1] = y[10];
    
    nxcross_ = nxcros;
    nzcross_ = nzcros;
    
    e += 1.0;
    }
    
//     std::cerr << std::setprecision(16) << y[2] << ", " << y[4] << ", 0.0, " << y[5] << " " << nxcross_ << std::endl;
//     std::cerr << std::setprecision(16) << y[7] << ", " << y[9] << ", 0.0, " << y[10] << " " << nzcross_ << std::endl;
    
}

//-----------------------------------------------------------------------------------------------------------------------------
// template<typename Value_type, typename Size_type, class Stepper>
// void ClosedOrbitFinder<Value_type, Size_type, Stepper>::christian_test(value_type accuracy, size_type maxit) {
//   
//   // READ IN MAGNETIC FIELD: ONLY FOR STAND-ALONE PROGRAM
//   int nsc = 8, nr = 141, Nth = 1440, nth = 1440/8; value_type r0 = 1.8, dr = 0.02;
//   bmag_ = MagneticField::malloc2df(Nth,nr);
//   MagneticField::ReadSectorMap(bmag_,nr,Nth,1,"data/ring590_bfld.dat",0.0);
//   MagneticField::MakeNFoldSymmetric(bmag_,Nth,nr,nth,nsc);
//   
//   // READ IN MAGNETIC FIELD: ONLY FOR STAND-ALONE PROGRAM
//   value_type bint, brint, btint; // B, dB/dr, dB/dtheta
//   
//   double ark[4]={0.5, 0.5, 1.0};
//   double brk[4]={1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
//   double crk[4]={0, 0.5, 0.5, 1};
//   double k[4*12];
//   double yold[12];
//   
//   h_.resize(N_);
//   ds_.resize(N_);
//   fidx_.resize(N_);
//   
//   double stp,psq,pp,gam,pt,pts,tpi,rg,prg;
//     double xc0,xc1,xc2,xc3,xc4,ep1,ep2,den,yp5,yp10;
//     double rp,prp,dth;
//     // Start values of R,Pr:
//     double rrk,y[12],ck[12],q[12];
//     double sect,de;
//     double phs;
//     double e,rav,rmax,rmin;
//     double cn,sn;
//     double rznt,rznu[2];
// //     double bint,brint,btint;
//     double tmp;
//     FILE *f;
//     char fname[200];
//     int nstp,ntry,neq,nxcros,nzcros,eqdst;
//     int i,j,n0;
//     int ir0,ith,istp;
//     int done=0;
//     /***/
//     sprintf(fname,"eo_params.dat");
//     f=fopen(fname,"w");
//     /***/
//     de=1.0;
//     eqdst=(de>0.0);
//     tpi=2.0*M_PI;
//     nth=nsc*nth;			// nsc = #sectors
//     dth=tpi/(double)nth;			// dth = angle step (d\theta)
//     sect=tpi;//(double)(cyc->nsc);
//     nstp=nth/2;
// //     n0=(int)((cyc->emax-cyc->emin)/de+1.5);
//     n0 = 0;
//     // Allocate Memory:
//     stp=tpi/(double)(nstp);
//     
//     rp=0.0;
//     prp=0.0;
//     double nuz=0.0;
//     double nur=0.0;
//     phs=0.0;
// //     ir0=(int)(1.5-cyc->r0/cyc->dr);
//     //dpeo[0]=phs;
//     //==============================
//     // PSI RING:
//     eqdst=1;
//     double emin = 72.0;
//     e=emin;
//     int neo = 0;
//     while(e < E_+de)
//     {
//     //==============================
// 	tmp=e/physics::E0;
// 	// Relativistic gamma-factor = (E+E0)/E0 
//         // because E0 is rest mass energy and E is kinetic energy 
// 	gam=1.0+tmp;
// 	// psq = squared total momentum 
// 	// 
// 	// P^2 = A^2 * (2 E/E0+(E/E0)^2)
// 	//
// 	// A="acon" is "cyclotron unit". 
// 	psq=tmp*(2.0+tmp)*physics::acon(wo_)*physics::acon(wo_);			// Gordon: formula (3) squared
// 	pp=std::sqrt(psq);						// Gordon: formula (3)
// 	// "ntry" is the iteration counter.
// 	//  ntry > 20 : No success.
// 	//  ntry < 0  : Last iteration (use ALL equations of motion)
// 	ntry=0;
// 	done=0;
// // 	cyc->E[cyc->neo]=e;
// // 	cyc->PC[cyc->neo]=cyc->e0*std::sqrt(tmp*(2.0+tmp));
// // 	cyc->gam[cyc->neo]=gam;
// // 	cyc->gam2[cyc->neo]=SQR(gam);
// 	if (neo<=1) {
// 	    rg=std::sqrt(psq);
// 	    prg=prp;
// 	}
// 	// Iterate until conversion:
// 	do {
// 	    // Set initial values for the integration variables:
// 	    // Radius:
// 	  y[0]=rg;
// 	  // Radial Momentum:
// 	  y[1]=prg;
// //     The Equations 3..6 describe the x-motion.
// //     There are two identical pairs of equations for
// //     two (orthogonal) pairs of starting conditions
// 	    y[2]=1.0;
// 	    y[3]=0.0;
// 	    y[4]=0.0;
// 	    y[5]=1.0;
// //     Equation 7 is the "Phase-Equation"...
// 	    y[6]=0.0;
// //     The Equations 8..11 describe the z-motion.
// //     There are two identical pairs of equations for
// //     two (orthogonal) pairs of starting conditions
// 	    y[7]=1.0;
// 	    y[8]=0.0;
// 	    y[9]=0.0;
// 	    y[10]=1.0;
// //     Equation 12 the Equation for Rav:
// 	    y[11]=0.0;
// 	    for (i=0;i<12;i++) q[i]=0.0;
// 	    istp=0;
// 	    nxcros=nzcros=0;
// 	    // Runge-Kutta Integration over Theta:
// 	    do {
// 	      
// // 	          std::cout << std::setprecision(15) << "r = " << y[0] << " pr = " << y[1] << std::endl;
// //     std::cout << std::setprecision(15) << "x1 = " << y[2] << " x2 = " << y[3] << " px1 = " << y[4] << " px2 = " << y[5] << std::endl; std::cin.get();
// 		r_[2*istp]=y[0];
// 		// Remember last values y[4] and y[9]:
// 		yp5=y[4];
// 		yp10=y[9];
// 		
// 		for(int i=0; i<12 ; ++i) {
// 		  yold[i] = y[i];
// 		}
// 		
// 		for (j=0;j<4;j++) { 
// 		    ith=((j+1)/2+2*istp);
// 		    MagneticField::interpolate1(&bint,&brint,&btint,ith,nr,nth,y[0],r0,dr,bmag_,ntry<0);
// 		    // "bint" = interpolated magn. field value
// 		    bint/=physics::bcon(wo_);
// 		    // "brint" = radial derivative of magn. field.
// 		    brint/=physics::bcon(wo_);
//                     // "btint" = theta derivative (last time only);
// 		    if (ntry<0) btint/=physics::bcon(wo_);
// 		    // "pts" = Ptheta^2
// 		    pts=psq-y[1]*y[1];
// 		    if (pts<=0.0) { 
// 		      std::cerr << "Error in ClosedOrbitFinder." << std::endl;
// 			return;
// 		    }
// 		    pt=std::sqrt(pts);
// 		    xc0=y[0]/pt;
// 	            // R - Equation:
// 		    // dR = dth * r/pt*pr;
// 		    ck[0]=xc0*y[1];
// 	            // Pr - Equation:
//                     // dPR = dth*(PT-Q*R*Bz)
// 		    // For the correct units, Pr is given in cyclotron units
// 		    // and Q*Bz in units of Bcon.
// 		    ck[1]=pt-y[0]*bint;
// 		    if ((ntry<0)&&(j==0)) {
// 		      // dR/dth:
// // 		      r_[2*istp]=ck[0];
// 		      // d2R/dth^2:
// 		      // Inverse Bending-Radius:
// 		      h_[2*istp]=bint/pp;
// 		      // Feldindex:
// 		      fidx_[2*istp]=(brint*pt-btint*y[1]/y[0])/psq;
// 		      ds_[2*istp]=std::sqrt(y[0]*y[0]+ck[0]*ck[0])*dth;
// 		    }
// 		    //     The Equations 3..6 describe the x-motion.
// 		    //     There are two identical pairs of equations for
// 		    //     two (orthogonal) pairs of starting conditions
// 		    //     Refer to M.M. Gordon in:
// 		    //     Particle Accelerators Vol. 16, pp 39-62 (1984)
// 		    xc1=y[1]/pt;
// 		    xc2=psq*xc0/pts;
// 		    xc3=y[0]*brint+bint;
// 		    // dx = dth ( (pr/pt) * x + R * p^2/pt^3 * px)
// 		    ck[2]= xc1*y[2] + xc2*y[3];
// 		    // dpx = - dth ( (pr/pt) * px + (r*dB/dr+B) * x)
// 		    ck[3]=-xc1*y[3] - xc3*y[2];
// 		    ck[4]= xc1*y[4] + xc2*y[5];
// 		    ck[5]=-xc1*y[5] - xc3*y[4];
// // 		    std::cout << std::setprecision(16) << ck[1] << std::endl; std::cin.get();
// 		    neq=6;
// 		    if (ntry<0) {
//                         //     phase and z equations (last time only);
// 			neq=12;
//                         //     Equation 7 is the "Phase-Equation" 
// 			// dphi=dth*(r/pt*gamma-1)
// 			ck[6]=xc0*gam-1.0;
//                         //     The Equations 8..11 describe the z-motion.
//                         //     There are two identical pairs of equations for
//                         //     two (orthogonal) pairs of starting conditions
// 			//  
// 			// dz=dth*(r/pt*pz)
// 			xc4=y[0]*brint-xc1*btint;
// 			ck[7]=xc0*y[8];
// 			// dpz=dth*(r*dB/dr-pr/pt*dB/dth)*z;
// 			ck[8]=xc4*y[7];
// 			// dz=dth*(r/pt*pz)
// 			ck[9]=xc0*y[10];
// 			// dpz=dth*(r*dB/dr-pr/pt*dB/dth)*z;
// 			ck[10]=xc4*y[9];
//                         //     Equation 12 is the Equation for Rav
// 			// drav=dth * r
// 			ck[11]=y[0];
// 		    }
// 		    // Normierung auf Schrittweite:
// // 		    for (i=0;i<neq;i++) ck[i]*=stp;
// 		    // Runge-Kutta-Integration:
// 		    for (int i=0;i<neq;i++) {
// 			k[j+4*i] = ck[i];
// 			y[i] = yold[i] + stp*ark[j]*k[j+4*i];
// 		    }
// 		}// for (j=0;j<4;j++)
// 		
// 		for(int i=0; i<neq; ++i) {
// 		  y[i] = yold[i] + stp*(brk[0]*k[0+4*i] + brk[1]*k[1+4*i] + brk[2]*k[2+4*i] + brk[3]*k[3+4*i]);
// 		}
// 		
// 		istp++;
// 		if (istp>1) {
// 		    if (yp5*y[4]<0.0) nxcros++;
// 		    if (yp10*y[9]<0.0) nzcros++;
// 		}
// 	    } while (istp<nstp);
// 	    if (ntry>=0) {
// 		// Compare to M.M. Gordon: "Computation of closed orbits and basic
// 		// focusing properties for sector-focused cyclotrons and the design
// 		// of cyclops", Part. Acc. Vol 16 (1984), pp. 39-62. Here page 45:
// 		ep1=y[0]-rg;  // "eps1"
// 		ep2=y[1]-prg; // "eps2"
// 		
//                 // Denominator (determinant von X(thf,thi)-1):
// 		den=(1.0-y[2])*(1.0-y[5])-y[3]*y[4];
// 		rg+=(ep1*(1.0-y[5])+ep2*y[4])/den;
// 		prg+=(ep1*y[3]+ep2*(1.0-y[2]))/den;
// 		if (std::fabs(ep1)+std::fabs(ep2)>=1.0E-8*rg) {
// 		    ntry++;
// 		} else ntry=-1;
// 	    } else done = 1;
// 	} while (!done);
// 	phs=y[6]/sect;
// 	rav=y[11]/sect;
// 	ravg_=y[11]/sect;
// 	phase_=y[6]/sect;
// // 	rmax=0.0;
// // 	rmin=1.0E20;
// 	for (i=0;i<nstp;i++) {
// 	  ith=i*2+1;
// 	  r_[ith % nth]=0.5*(r_[(ith-1) % nth]+r_[(ith+1) % nth]);
// // 	  cyc->r1eo[cyc->neo][ith % nth]=0.5*(cyc->r1eo[cyc->neo][(ith-1) % nth]+cyc->r1eo[cyc->neo][(ith+1) % nth]);
// 	  h_[ith % nth]=0.5*(h_[(ith-1) % nth]+h_[(ith+1) % nth]);
// 	  fidx_[ith % nth]=0.5*(fidx_[(ith-1) % nth]+fidx_[(ith+1) % nth]);
// 	  ds_[ith % nth]=0.5*(ds_[(ith-1) % nth]+ds_[(ith+1) % nth]);
// 	}
// 	neo++;
// // 	if (cyc->neo<n0) {
// 	    ep1=rg-rp;
// 	    ep2=prg-prp;
// 	    rp=rg;
// 	    rg+=ep1;
// 	    prp=prg;
// 	    prg+=ep2;
// // 	}
// //     } while (cyc->neo<n0);
//     x_[0] = y[2];
//     px_[0] = y[3];
//     
//     x_[1] = y[4];
//     px_[1] = y[5];
//     
//     z_[0] = y[7];
//     pz_[0] = y[8];
//     z_[1] = y[9];
//     pz_[1] = y[10];
//     
//     nxcross_ = nxcros;
//     nzcross_ = nzcros;
//     
//     e += 1.0;
//     }
// }

#endif