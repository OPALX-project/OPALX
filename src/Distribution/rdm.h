#ifndef RDM_H
#define RDM_H

#include <cmath>
#include <iostream>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "error.h"
#include "physical_error.h"
#include "matrix_vector_operation.h"

// Remark: Comments starting with "///" are doxygen comments, comments with "//" are normal comments.

#ifdef DEBUG
  #include <boost/numeric/ublas/io.hpp>
#endif

/// Real Dirac Matrix class
/*!
 * They're ordered after the paper of Dr. C. Baumgarten: "Use of real Dirac matrices in two-dimensional coupled linear optics".
 * The diagonalizing method is based on the paper "Geometrical method of decoupling" (2012) of Dr. C. Baumgarten.
 * The template parameter Value_type is the datatype of the matrices and variables.
 */
template<typename Value_type, typename Size_type>
class RDM
{
public:
  /// Type of variables
  typedef Value_type value_type;
  /// Type for specifying sizes
  typedef Size_type size_type;
  /// Sparse matrix type definition
  typedef boost::numeric::ublas::compressed_matrix<value_type,boost::numeric::ublas::row_major> sparse_matrix_type;
  /// Dense matrix type definition
  typedef boost::numeric::ublas::matrix<value_type> matrix_type;
  /// Dense vector type definition
  typedef boost::numeric::ublas::vector<value_type> vector_type;
  
  /// Returns the i-th Real Dirac matrix
  /*!
   * @param i specifying the matrix (has to be in the range from 0 to 15)
   */
  sparse_matrix_type getRDM(short);
  
  /// Decomposes a real-valued 4x4 matrix into a linear combination and returns a vector containing the coefficients
  /*!
   * @param M an arbitrary real-valued 4x4 matrix
   */
  vector_type decompose(const matrix_type&);
  
  /// Takes a vector of coefficients, evaluates the linear combination of RDMs with these coefficients and returns a 4x4 matrix
  /*!
   * @param coeffs is a vector of coefficients (at most length NumOfRDMs)
   */
  matrix_type combine(const vector_type&);
  
  /// Brings a 4x4 symplex matrix into Hamilton form and computes the transformation matrix and its inverse
  /*!
   * @param Ms is a 4x4 symplex matrix
   * @param R is the 4x4 transformation matrix (gets computed)
   * @param invR is the 4x4 inverse transformation matrix (gets computed)
   */
  void diagonalize(matrix_type&, sparse_matrix_type&, sparse_matrix_type&);
  
  /// Returns the symplex part of a 4x4 real-valued matrix
  /*!
   * @param M 4x4 real-valued matrix
   */
  matrix_type symplex(const matrix_type&);
  
  /// Returns the cosymplex part of a 4x4 real-valued matrix
  /*!
   * @param M 4x4 real-valued matrix
   */
  matrix_type cosymplex(const matrix_type&);
  
  /// The number of real Dirac matrices
  static short NumOfRDMs;
  /// The matrix dimension (4x4)
  static short DimOfRDMs;
  
private:
  /// Applies a rotation to the matrix M by a given angle
  /*!
   * @param M is the matrix to be transformed
   * @param i is the i-th RDM used for transformation
   * @param phi is the angle of rotation
   * @param Rtot is a reference to the current transformation matrix
   * @param invRtot is a reference to the inverse of the current transformation matrix
   */
  void transform(matrix_type&, short, value_type, sparse_matrix_type&, sparse_matrix_type&);
};

// -----------------------------------------------------------------------------------------------------------------------
// PUBLIC STATIC MEMBER VARIABLES
// -----------------------------------------------------------------------------------------------------------------------

short NumOfRDMs = 16;
short DimOfRDMs = 4;

// -----------------------------------------------------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type>
typename RDM<Value_type, Size_type>::sparse_matrix_type RDM<Value_type, Size_type>::getRDM(short i) {
  sparse_matrix_type rdm(4,4,4);	// #nrows, #ncols, #non-zeros
  switch(i) {
    case 0: rdm(0,1) = rdm(2,3) = 1; rdm(1,0) = rdm(3,2) = -1; 	break;
    case 1: rdm(0,1) = rdm(1,0) = -1; rdm(2,3) = rdm(3,2) = 1; 	break;
    case 2: rdm(0,3) = rdm(1,2) = rdm(2,1) = rdm(3,0) = 1; 	break;
    case 3: rdm(0,0) = rdm(2,2) = -1; rdm(1,1) = rdm(3,3) = 1; 	break;
    case 4: rdm(0,0) = rdm(3,3) = -1; rdm(1,1) = rdm(2,2) = 1;	break;
    case 5: rdm(0,2) = rdm(2,0) = 1; rdm(1,3) = rdm(3,1) = -1; 	break;
    case 6: rdm(0,1) = rdm(1,0) = rdm(2,3) = rdm(3,2) = 1;	break;
    case 7: rdm(0,3) = rdm(2,1) = 1; rdm(1,2) = rdm(3,0) = -1;	break;
    case 8: rdm(0,1) = rdm(3,2) = 1; rdm(1,0) = rdm(2,3) = -1;	break;
    case 9: rdm(0,2) = rdm(1,3) = -1; rdm(2,0) = rdm(3,1) = 1;	break;
    case 10: rdm(0,2) = rdm(3,1) = 1; rdm(1,3) = rdm(2,0) = -1;	break;
    case 11: rdm(0,2) = rdm(1,3) = rdm(2,0) = rdm(3,1) = -1;	break;
    case 12: rdm(0,0) = rdm(1,1) = -1; rdm(2,2) = rdm(3,3) = 1;	break;
    case 13: rdm(0,3) = rdm(3,0) = -1; rdm(1,2) = rdm(2,1) = 1;	break;
    case 14: rdm(0,3) = rdm(1,2) = -1; rdm(2,1) = rdm(3,0) = 1;	break;
    case 15: rdm(0,0) = rdm(1,1) = rdm(2,2) = rdm(3,3) = 1;	break;
    default: Error::message("RDM<Value_type, Size_type>::generateRDM(short& i)",Error::notdefined); break;
  }
  return rdm;
}

template<typename Value_type, typename Size_type>
typename RDM<Value_type, Size_type>::vector_type RDM<Value_type, Size_type>::decompose(const matrix_type& M) {
  /*
   * Formula (11) from paper:
   * Geometrical method of decoupling
   */
  if(M.size1() != 4 && M.size1() != M.size2())
    Error::message("RDM<Value_type, Size_type>::decompose(const matrix_type&)", Error::dim);
  
  vector_type coeffs(16);
  
  matrix_type left(4,4), right(4,4), rdm2(4,4);
  
  value_type inv32 = 1.0/32.0;
  
  // do it with iterators
  for(short i=0; i<16; ++i) {
    sparse_matrix_type rdm = getRDM(i);
    left = boost::numeric::ublas::prod(M,rdm);
    right = boost::numeric::ublas::prod(rdm,M);
    rdm2 = boost::numeric::ublas::prod(rdm,rdm);
    left *= inv32;
    right *= inv32;
    left += right;
    coeffs(i) = matt_boost::trace(rdm2)*matt_boost::trace(left);
  }
  return coeffs;
}

template<typename Value_type, typename Size_type>
typename RDM<Value_type, Size_type>::matrix_type RDM<Value_type, Size_type>::combine(const vector_type& coeffs) {
  if(coeffs.size() > 16)
    Error::message("RDM<Value_type, Size_type>::combine(const vector_type&)",Error::size);
  
  // initialize a 4x4 zero matrix
  matrix_type M = boost::numeric::ublas::zero_matrix<value_type>(4,4);
  
  // evaluate linear combination
  for(short i=0; i<16; ++i)
    M += coeffs(i)*getRDM(i);
  
  return M;
}

template<typename Value_type, typename Size_type>
void RDM<Value_type, Size_type>::diagonalize(matrix_type& Ms, sparse_matrix_type& R, sparse_matrix_type& invR) {
  
  // R and invR store the total transformation
  R = boost::numeric::ublas::identity_matrix<value_type>(4);
  invR = R;
  
  vector_type P(3), E(3), B(3), b;
  value_type mr, mg, mb, eps;
  
  // Lambda function to compute vectors E, P, B and scalar eps (it takes the current Ms as reference argument (--> [&])
  auto mult = [&](short i) {
    /*
     * For computing E, P, B, eps according to formula (C4) from paper:
     * Geometrical method of decoupling
     */
    matrix_type tmp = boost::numeric::ublas::prod(Ms,getRDM(i))+boost::numeric::ublas::prod(getRDM(i),Ms);
    return 0.125*matt_boost::trace(tmp);
  };
  
  // 1. Transformation with \gamma_{0}
  P(0) = mult(1); P(1) = mult(2); P(2) = mult(3);
  E(0) = mult(4); E(1) = mult(5); E(2) = mult(6);
  B(0) = -mult(7); B(1) = -mult(8); B(2) = -mult(9);
  mr = boost::numeric::ublas::inner_prod(E,B);		// formula (31), paper: Geometrical method of decoupling
  mg = boost::numeric::ublas::inner_prod(B,P);		// formula (31), paper: Geometrical method of decoupling
  
#ifdef DEBUG
  eps = -mult(0);
  P(0) = mult(1); P(1) = mult(2); P(2) = mult(3);
  E(0) = mult(4); E(1) = mult(5); E(2) = mult(6);
  B(0) = -mult(7); B(1) = -mult(8); B(2) = -mult(9);
  mr = boost::numeric::ublas::inner_prod(E,B);
  mg = boost::numeric::ublas::inner_prod(B,P);
  mb = boost::numeric::ublas::inner_prod(E,P);
  b = eps*B+matt_boost::cross_prod(E,P);
  std::cout << "Mr = " << mr << " Mb = " << mb << " Mg = " << mg << std::endl;
  std::cout << "b = " << b << std::endl << "P = " << P << std::endl;
  std::cout << "E = " << E << std::endl << "B = " << B << std::endl;
#endif
  
  transform(Ms,short(0),0.5*std::atan2(mg,mr),R,invR);
  
#ifdef DEBUG
    std::cout << "rdm0: mg, mr" << std::endl << Ms << std::endl; std::cin.get();
#endif
  
  // 2. Transformation with \gamma_{7}
  eps = -mult(0);
  P(0) = mult(1); P(1) = mult(2); P(2) = mult(3);
  E(0) = mult(4); E(1) = mult(5); E(2) = mult(6);
  B(0) = -mult(7); B(1) = -mult(8); B(2) = -mult(9);
  b = eps*B+matt_boost::cross_prod(E,P);		// formula (32), paper: Geometrical method of decoupling
  
#ifdef DEBUG
//   std::cout << decompose(Ms) << std::endl;
  eps = -mult(0);
  P(0) = mult(1); P(1) = mult(2); P(2) = mult(3);
  E(0) = mult(4); E(1) = mult(5); E(2) = mult(6);
  B(0) = -mult(7); B(1) = -mult(8); B(2) = -mult(9);
  mr = boost::numeric::ublas::inner_prod(E,B);
  mg = boost::numeric::ublas::inner_prod(B,P);
  mb = boost::numeric::ublas::inner_prod(E,P);
  b = eps*B+matt_boost::cross_prod(E,P);
  std::cout << "Mr = " << mr << " Mb = " << mb << " Mg = " << mg << std::endl; 
  std::cout << "b = " << b << std::endl << "P = " << P << std::endl;
  std::cout << "E = " << E << std::endl << "B = " << B << std::endl;
#endif
  
  transform(Ms,short(7),0.5*std::atan2(b(2),b(1)),R,invR);
  
#ifdef DEBUG
  std::cout << "rdm7: b(2), b(1)" << std::endl << Ms << std::endl; std::cin.get();
#endif
  
  // 3. Transformation with \gamma_{9}
  eps = -mult(0);
  P(0) = mult(1); P(1) = mult(2); P(2) = mult(3);
  E(0) = mult(4); E(1) = mult(5); E(2) = mult(6);
  B(0) = -mult(7); B(1) = -mult(8); B(2) = -mult(9);
  b = eps*B+matt_boost::cross_prod(E,P);  
  
#ifdef DEBUG
  eps = -mult(0);
  P(0) = mult(1); P(1) = mult(2); P(2) = mult(3);
  E(0) = mult(4); E(1) = mult(5); E(2) = mult(6);
  B(0) = -mult(7); B(1) = -mult(8); B(2) = -mult(9);
  mr = boost::numeric::ublas::inner_prod(E,B);
  mg = boost::numeric::ublas::inner_prod(B,P);
  mb = boost::numeric::ublas::inner_prod(E,P);
  b = eps*B+matt_boost::cross_prod(E,P);
//   std::cout << b(0) << " " << b(1) << std::endl;
  std::cout << "Mr = " << mr << " Mb = " << mb << " Mg = " << mg << std::endl; 
  std::cout << "b = " << b << std::endl << "P = " << P << std::endl;
  std::cout << "E = " << E << std::endl << "B = " << B << std::endl;
#endif
  
  transform(Ms,short(9),-0.5*std::atan2(b(0),b(1)),R,invR);
  
#ifdef DEBUG
  std::cout << "rdm9: b(0), b(1)" << std::endl << Ms << std::endl; std::cin.get();
#endif
  
  // 4. Transformation with \gamma_{2}
  eps = -mult(0);
  P(0) = mult(1); P(1) = mult(2); P(2) = mult(3);
  E(0) = mult(4); E(1) = mult(5); E(2) = mult(6);
  B(0) = -mult(7); B(1) = -mult(8); B(2) = -mult(9);
  mr = boost::numeric::ublas::inner_prod(E,B);
  b = eps*B+matt_boost::cross_prod(E,P);  
  
#ifdef DEBUG
  eps = -mult(0);
  P(0) = mult(1); P(1) = mult(2); P(2) = mult(3);
  E(0) = mult(4); E(1) = mult(5); E(2) = mult(6);
  B(0) = -mult(7); B(1) = -mult(8); B(2) = -mult(9);
  mr = boost::numeric::ublas::inner_prod(E,B);
  mg = boost::numeric::ublas::inner_prod(B,P);
  mb = boost::numeric::ublas::inner_prod(E,P);
  b = eps*B+matt_boost::cross_prod(E,P);
  std::cout << "min = " << std::fmin(std::fabs(mr),std::fabs(b(1))) << " max = " << std::fmax(std::fabs(mr),std::fabs(b(1))) << std::endl;
  std::cout << "Mr = " << mr << " Mb = " << mb << " Mg = " << mg << std::endl; 
  std::cout << "b = " << b << std::endl << "P = " << P << std::endl;
  std::cout << "E = " << E << std::endl << "B = " << B << std::endl;
#endif
  
  if(std::fabs(mr) < std::fabs(b(1))) {
    transform(Ms,short(2),0.5*std::atanh(mr/b(1)),R,invR);
  } else {
    transform(Ms,short(2),0.5*std::atanh(b(1)/mr),R,invR);
  }
  
#ifdef DEBUG
  std::cout << "rdm2: b(1), mr" << std::endl << Ms << std::endl; std::cin.get();
#endif  
  
//   mr = boost::numeric::ublas::inner_prod(E,B);
  
  eps = -mult(0);
  P(0) = mult(1); P(1) = mult(2); P(2) = mult(3);
  E(0) = mult(4); E(1) = mult(5); E(2) = mult(6);
  B(0) = -mult(7); B(1) = -mult(8); B(2) = -mult(9);
  
  // formula (31), paper: Geometrical method of decoupling
  mr = boost::numeric::ublas::inner_prod(E,B);
  mg = boost::numeric::ublas::inner_prod(B,P);
  mb = boost::numeric::ublas::inner_prod(E,P);
  
  // formula (29), paper: Geometrical method of decoupling
//   value_type EP2 = boost::numeric::ublas::inner_prod(E,P); EP2 *= EP2;
//   value_type EB2 = boost::numeric::ublas::inner_prod(E,B); EB2 *= EB2;
//   value_type PB2 = boost::numeric::ublas::inner_prod(P,B); PB2 *= PB2;
  
//   vector_type eb = matt_boost::cross_prod(E,B);
  
  value_type P2 = boost::numeric::ublas::inner_prod(P,P);
  value_type E2 = boost::numeric::ublas::inner_prod(E,E);
//   value_type B2 = boost::numeric::ublas::inner_prod(B,B);
  
//   value_type k2 = eps*eps*B2+E2*P2-EP2-EB2-PB2-2.0*eps*boost::numeric::ublas::inner_prod(P,eb);
//   value_type k1 = eps*eps+B2-E2-P2;
  
  // check if stable system (see p. 5, left column at the top, paper: Geometrical method of decoupling)
/*   if(k2 < 0) {
#ifdef DEBUG
    std::cout << "k2 = " << k2 << std::endl;
#endif
    PhysicalError::message("void RDM<Value_type, Size_type>::diagonalize(matrix_type&)",PhysicalError::match);
  }*/
  // Transformation with \gamma_{0}
//   mb = boost::numeric::ublas::inner_prod(E,P);
  
#ifdef DEBUG
  eps = -mult(0);
  P(0) = mult(1); P(1) = mult(2); P(2) = mult(3);
  E(0) = mult(4); E(1) = mult(5); E(2) = mult(6);
  B(0) = -mult(7); B(1) = -mult(8); B(2) = -mult(9);
  mr = boost::numeric::ublas::inner_prod(E,B);
  mg = boost::numeric::ublas::inner_prod(B,P);
  mb = boost::numeric::ublas::inner_prod(E,P);
  b = eps*B+matt_boost::cross_prod(E,P);
  std::cout << "Mr = " << mr << " Mb = " << mb << " Mg = " << mg << std::endl; 
  std::cout << "b = " << b << std::endl << "P = " << P << std::endl;
  std::cout << "E = " << E << std::endl << "B = " << B << std::endl;
#endif
  
  transform(Ms,short(0),0.25*std::atan2(mb,0.5*(E2-P2)),R,invR);
  
#ifdef DEBUG
  std::cout << "rdm0: mb, 0.5*(E2-P2)" << std::endl << Ms << std::endl; std::cin.get();
#endif
  
  // Transformation with \gamma_{8}
  P(0) = mult(1); P(2) = mult(3);
  
#ifdef DEBUG
  eps = -mult(0);
  P(0) = mult(1); P(1) = mult(2); P(2) = mult(3);
  E(0) = mult(4); E(1) = mult(5); E(2) = mult(6);
  B(0) = -mult(7); B(1) = -mult(8); B(2) = -mult(9);
  mr = boost::numeric::ublas::inner_prod(E,B);
  mg = boost::numeric::ublas::inner_prod(B,P);
  mb = boost::numeric::ublas::inner_prod(E,P);
  b = eps*B+matt_boost::cross_prod(E,P);
  std::cout << "Mr = " << mr << " Mb = " << mb << " Mg = " << mg << std::endl; 
  std::cout << "b = " << b << std::endl << "P = " << P << std::endl;
  std::cout << "E = " << E << std::endl << "B = " << B << std::endl;
#endif
  
  transform(Ms,short(8),-0.5*std::atan2(P(2),P(0)),R,invR);
  
#ifdef DEBUG
  std::cout << "rdm8: P(2), P(0)" << std::endl << std::fixed << Ms << std::endl;
  std::cout << "R: " << R << std::endl;
  std::cout << "invR: " << invR << std::endl; std::cin.get();
#endif
  
}

template<typename Value_type, typename Size_type>
typename RDM<Value_type, Size_type>::matrix_type RDM<Value_type, Size_type>::symplex(const matrix_type& M) {
  sparse_matrix_type rdm0 = getRDM(0);
//   matrix_type Mt = boost::numeric::ublas::trans(M);
  return 0.5*(M+matt_boost::gemmm<matrix_type>(rdm0,boost::numeric::ublas::trans(M)/*Mt*/,rdm0));
}


template<typename Value_type, typename Size_type>
typename RDM<Value_type, Size_type>::matrix_type RDM<Value_type, Size_type>::cosymplex(const matrix_type& M) {
  sparse_matrix_type rdm0 = getRDM(0);
//   matrix_type Mt = boost::numeric::ublas::trans(M);
  return 0.5*(M-matt_boost::gemmm<matrix_type>(rdm0,boost::numeric::ublas::trans(M)/*Mt*/,rdm0));
}

// -----------------------------------------------------------------------------------------------------------------------
// PRIVATE MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type>
void RDM<Value_type, Size_type>::transform(matrix_type& M, short i, value_type phi, sparse_matrix_type& Rtot, sparse_matrix_type& invRtot) {
  
/*   if(i == 7 || i == 8 || i ==  9) {
    if (phi > 0.25*M_PI) phi -= 0.5*M_PI;
    if (phi < -0.25*M_PI) phi += 0.5*M_PI;
  }*/
#ifdef DEBUG
  std::cout << "phi = " << phi << std::endl;
#endif
  if(phi) {	// if phi == 0 --> nothing happens, since R and invR would be identity_matrix matrix
    sparse_matrix_type R(4,4), invR(4,4);
    sparse_matrix_type I = boost::numeric::ublas::identity_matrix<value_type>(4);
//     sparse_matrix_type I = matt_boost::prod(getRDM(i),getRDM(i));
//     I(0,0) = I(1,1) = I(2,2) = I(3,3) = 1;
    
    if(i < 7 && i != 0 && i < 11 && i != 14) {
      R = I*std::cosh(phi) + getRDM(i)*std::sinh(phi);
      invR = I*std::cosh(phi) - getRDM(i)*std::sinh(phi);
    } else {
      R = I*std::cos(phi) + getRDM(i)*std::sin(phi);
      invR = I*std::cos(phi) - getRDM(i)*std::sin(phi);
    }
    // update matrices
    M = matt_boost::gemmm<matrix_type>(R,M,invR);
    Rtot = boost::numeric::ublas::prod(R,Rtot);
    invRtot = boost::numeric::ublas::prod(invRtot,invR);
  }
}

#endif