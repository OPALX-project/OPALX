#ifndef MAPGENERATOR_H
#define MAPGENERATOR_H

#include <cmath>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>

#include "matrix_vector_operation.h"
#include "physics.h"

#ifdef PRINT
#include <fstream>
#endif

/// This class generates the matrices for the one turn matrix of a cyclotron.
/*!
 * These computations are based on the paper "Transverse-Longitudinal Coupling by Space Charge in Cyclotrons" (2012)
 * of Dr. Christian Baumgarten.
 * It has one template parameter that specifies the type of the variables and containers.
 */
template<typename Value_type, typename Size_type>
class MapGenerator
{
public:
  /// Type of variables
  typedef Value_type value_type;
  /// Type for specifying sizes
  typedef Size_type size_type;
  /// Type for specifying matrices
  typedef boost::numeric::ublas::matrix<value_type> matrix_type;
  /// Type for specifying vectors
  typedef std::vector<value_type> vector_type;
  
  /// Constructor
  /*!
   * @param N is the number of angle splits
   * @param I is the current
   * @param gamma is the relativistic factor
   * @param wo is the orbital frequency
   * @param nh is the harmonic number
   */
  MapGenerator(size_type, value_type, value_type, value_type, value_type);
  
  /// Computes the space charge map for a given angle
  /*!
   * @param sigma is the sigma matrix
   * @param ds is the step size (angle dependent)
   */
  matrix_type generateMsc(const matrix_type&, value_type);
  
  /// Computes the cyclotron map for a given angle
  /*!
   * @param h is the inverse bending radius (angle dependent)
   * @param fidx is the field index (angle dependent)
   * @param ds is the step size (angle dependent)
   */
  matrix_type generateMcyc(value_type,value_type,value_type,size_type);
  
  /// Combines the space charge maps (for each angle one) and the cyclotron maps (for each angle one) to the one turn map
  /*!
   * @param h is vector of inverse bending radiuses for each angle
   * @param fidx is a vector of field indexes for each angle
   * @param ds is vector of path lengths for each angle
   * @param sigma is a list of sigma matrices for each angle
   */
  void combine(const vector_type&, const vector_type&, const vector_type&, /*const*/ std::vector<matrix_type>&, size_type);
  
  /// Combines the space charge maps (for each angle one) and the cyclotron maps (for each angle one) to the ont turn map, taking lists of maps
  /*!
   * @param Mscs is a list of space charge maps (the higher the index, the higher the angle)
   * @param Mcycs is a list of cyclotron maps (the higher the index, the higher the angle)
   */
  void combine(std::vector<matrix_type>&, std::vector<matrix_type>&);
  
  /// Returns the one turn map
  matrix_type getMap();
  
private:
  /// Number of angle splits
  size_type N_;
  /// Relativistic factor
  value_type gamma_;
  /// Relativistic factor squared
  value_type gamma2_;
  /// One-turn matrix
  matrix_type Mturn_;
  /// Constant
  value_type K3_;
};

// -----------------------------------------------------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type>
MapGenerator<Value_type, Size_type>::MapGenerator(size_type N, value_type I, value_type gamma, value_type wo, value_type nh)
: N_(N), gamma_(gamma), gamma2_(gamma*gamma), Mturn_(6,6)
{
  K3_ = 3.0*physics::q0*I*1e-6/(10.0*std::sqrt(5.0)*physics::eps0*physics::E0*wo*nh*gamma_*(gamma2_-1.0));
#ifdef DEBUG
  // K3 = 6.5062e-10
  std::cout << "K3 = " << K3_ << std::endl;
#endif
}

template<typename Value_type, typename Size_type>
typename MapGenerator<Value_type, Size_type>::matrix_type MapGenerator<Value_type, Size_type>::generateMsc(const matrix_type& sigma, value_type ds) {
  /*
   * Computes space charge matrix according to 
   * formula (20) and (21) from paper:
   * Transverse-Longitudinal Coupling by Space Charge in Cyclotrons
   * 
   * The functions "__M6SC" and "__M6sc" of mtx.c of Dr. Christian Baumgarten are used as implementation guide.
   */
  // space charge map for one angle
  matrix_type Msc = boost::numeric::ublas::identity_matrix<value_type>(6);
  
  value_type sx = std::sqrt(std::fabs(sigma(0,0)))*0.001;
  value_type sy = std::sqrt(std::fabs(sigma(2,2)))*0.001;
  value_type sz = std::sqrt(std::fabs(sigma(4,4)))*0.001;
  
  value_type tmp = sx*sy;
  
  value_type f = std::sqrt(tmp)/(3.0*gamma_*sz);
  value_type kxy = K3_*std::fabs(1.0-f)*ds/((sx+sy)*sz);
  
  Msc(1,0) = kxy/sx;
  Msc(3,2) = kxy/sy;
  Msc(5,4) = K3_*f*ds*gamma2_/(tmp*sz);
  
#ifdef PRINT
  std::ofstream out("data/maps/SpaceChargeMapPerAngle.dat",std::ios::app);
  out << Msc << std::endl;
#endif
  
  return Msc;
}

template<typename Value_type, typename Size_type>
typename MapGenerator<Value_type, Size_type>::matrix_type MapGenerator<Value_type, Size_type>::generateMcyc(value_type h,value_type fidx,
													    value_type ds, size_type order) {
  
  /*
   * Computes the cyclotron matrix according to formula
   * (20) and (44) from the paper:
   * Transverse-Longitudinal Coupling by Space Charge in Cyclotrons
   * 
   * The implementation in the function "scMakeCycFromFieldData" of ring_ana5.c
   * of Dr. Christian Baumgarten is used as implementation guide
   */
  
  // force matrix for one angle
  matrix_type F = boost::numeric::ublas::zero_matrix<value_type>(6);
  F(0,1) = 1.0;
  F(1,0) = -h*h-fidx;	// p. 6 of paper: "Transverse-Longitudinal Coupling by Space Charge in Cyclotrons"
  F(1,5) = h;
  F(2,3) = 1.0;
  F(3,2) = fidx;
  F(4,0) = -h;
  F(4,5) = 1.0/gamma2_;
  
  // cyclotron map for one angle;
  matrix_type Mcyc = matt_boost::taylor_exp<value_type>(F,ds,order);
  
#ifdef PRINT
  std::ofstream out("data/maps/CyclotronMapPerAngle.dat",std::ios::app);
  out << Mcyc << std::endl;
#endif
  
  return Mcyc;
}

template<typename Value_type, typename Size_type>
void MapGenerator<Value_type, Size_type>::combine(const vector_type& h, const vector_type& fidx, const vector_type& ds,
			      /*const*/ std::vector<matrix_type>& sigma, size_type order) {
    
  Mturn_ = boost::numeric::ublas::identity_matrix<value_type>(6);
  
//   matrix_type tmp(6,6);
  
#ifdef PRINT
  std::ofstream out("data/maps/OneTurnMapSummedUp.dat");
#endif
  
  for(size_type i=0; i<N_; ++i) {
    Mturn_ = matt_boost::gemmm<matrix_type>(generateMsc(sigma[i],ds[i]),generateMcyc(h[i],fidx[i],ds[i],order),Mturn_);
#ifdef PRINT
    out << Mturn_ << std::endl;
#endif
//     tmp = boost::numeric::ublas::prod(generateMsc(sig,ds[i]),generateMcyc(h[i],fidx[i],ds[i]));
//     Mturn_ = boost::numeric::ublas::prod(tmp,Mturn_);
  }
#ifdef PRINT
  out.close();
#endif
}

template<typename Value_type, typename Size_type>
void MapGenerator<Value_type, Size_type>::combine(std::vector<matrix_type>& Mscs, std::vector<matrix_type>& Mcycs) {
  if(N_ != Mscs.size() || N_ != Mcycs.size()) {
    Error::message("MapGenerator<Value_type, Size_type>:combine(std::vector<matrix_type>&, std::vector<matrix_type>&", Error::size);
  }
  
  Mturn_ = boost::numeric::ublas::identity_matrix<value_type>(6);
  for(size_type i=0; i<N_; ++i) {
    Mturn_ = matt_boost::gemmm<matrix_type>(Mscs[i],Mcycs[i],Mturn_);
  }
}

template<typename Value_type, typename Size_type>
typename MapGenerator<Value_type, Size_type>::matrix_type MapGenerator<Value_type, Size_type>::getMap() {
  return Mturn_;
}

#endif