#ifndef MATRIX_OPERATION_H
#define MATRIX_OPERATION_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include <boost/numeric/ublas/vector.hpp>

#include "error.h"

/*!
 *  \addtogroup matt_boost
 *  @{
 */

/// Expands the existing functions of the boost library uBLAS (http://www.boost.org/).
namespace matt_boost {
  using namespace boost::numeric::ublas;
  
  /// Computes the trace of a square matrix
  template<class V>
  BOOST_UBLAS_INLINE
  V trace(matrix<V> &e) {
    V tr = 0;
    if(e.size1()==e.size2()) {
      matrix_vector_range<matrix<V> > diag(e,range(0,e.size1()),range(0,e.size2()));
      tr = sum(diag);
    }
    else
      Error::message("matt_boost::trace",Error::dim);
    return tr;
  }
  
  /// Computes the cross product \f$ v_{1}\times v_{2}\f$ of two vectors in \f$ \mathbb{R}^{3} \f$
  template<class V>
  BOOST_UBLAS_INLINE
  vector<V> cross_prod(vector<V> &v1, vector<V> &v2) {
    vector<V> v(v1.size());
    if(v1.size() == v2.size() && v1.size() == 3) {
      v(0) = v1(1)*v2(2)-v1(2)*v2(1);
      v(1) = v1(2)*v2(0)-v1(0)*v2(2);
      v(2) = v1(0)*v2(1)-v1(1)*v2(0);
    }
    else
      Error::message("matt_boost::cross_prod",Error::dim);
    return v;
  }
  
  /// Computes Taylor-Series of M(s) = exp(F*s)
  template<class V>
  BOOST_UBLAS_INLINE
  matrix<V> taylor_exp(const matrix<V>& F, const V ds, const unsigned int order) {
    double fac = 1.0;
    matrix<V> Fn = identity_matrix<V>(6);
    matrix<V> M = Fn;
    
    for(unsigned int k=1; k<order; ++k) {
      fac *= ds/V(k);
      Fn = prod(Fn,F);
      M = M + fac*Fn;
    }
    return M;
  }
  
  /// Generalized matrix-matrix-matrix multiplication \f$ e_{1}\cdot e_{2}\cdot e_{3} \f$
  template<class M, class E1, class E2, class E3>
  BOOST_UBLAS_INLINE
  M gemmm(const matrix_expression<E1>& e1, const matrix_expression<E2>& e2, const matrix_expression<E3>& e3) {
    M tmp = prod(e2,e3);
    return prod(e1,tmp);
  }
  
//   template<class V>
//   BOOST_UBLAS_INLINE
//   matrix<V> gemmm(const matrix<V>& A, const matrix<V>& B, const matrix<V>& C) {
//     matrix<V> tmp(A.size1(),C.size2());
//     if(A.size2() == B.size1() && B.size2() == C.size1()) {
//       matrix<V> tmp1 = prod(A,B);
//       tmp = prod(tmp1,C);
//     } else
//       Error::message("matt_boos::gemmm",Error::dim);
//     return tmp;
//   }
//   
//     /// Compute a matrix-matrix-matrix multiplication A*B*C
//   template<class V>
//   BOOST_UBLAS_INLINE
//   matrix<V> gemmm(const compressed_matrix<V,row_major>& A, const matrix<V>& B, const compressed_matrix<V,row_major>& C) {
//     matrix<V> tmp(A.size1(),C.size2());
//     if(A.size2() == B.size1() && B.size2() == C.size1()) {
//       matrix<V> tmp1 = prod(A,B);
//       tmp = prod(tmp1,C);
//     } else
//       Error::message("matt_boost::gemmm",Error::dim);
//     return tmp;
//   }
  
}

/*! @} End of Doxygen Groups*/

#endif