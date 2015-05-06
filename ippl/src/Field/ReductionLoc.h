// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef REDUCTION_LOC_H
#define REDUCTION_LOC_H

// include files
#include "PETE/IpplExpressions.h"

// forward declarations
template<unsigned D> class NDIndex;
template<class T, unsigned int D> class BareFieldIterator;
template<class T, unsigned int D> class IndexedBareFieldIterator;

//////////////////////////////////////////////////////////////////////

// Fundamental reduction function for minloc and maxloc

template<class T, class CompOp, class AccOp, unsigned D>
typename T::PETE_Return_t
Reduction(const PETE_Expr<T>& const_expr,
	  CompOp comp_op,
	  AccOp  acc_op,
	  NDIndex<D>& loc);

// Struct for saving current location

template <unsigned int D>
struct FindLocation
{
  typedef int PETE_Return_t;
  NDIndex<D>* Loc;
  FindLocation(NDIndex<D>& loc) : Loc(&loc) {}
};

// Behavior of this struct

template <class T, class C, unsigned int D>
inline int
for_each(const IndexedBareFieldIterator<T,D>& expr,
         FindLocation<D>& find_loc, C)
{
  int loc[D];
  expr.GetCurrentLocation(loc);
  for (unsigned d=0; d<D; d++) (*(find_loc.Loc))[d] = Index(loc[d],loc[d]);
  return 0;
}

template <class T, class C, unsigned int D>
inline int
for_each(const BareFieldIterator<T,D>& expr, FindLocation<D>& find_loc, C)
{
  int loc[D];
  expr.GetCurrentLocation(loc);
  for (unsigned d=0; d<D; d++) (*(find_loc.Loc))[d] = Index(loc[d],loc[d]);
  return 0;
}

//////////////////////////////////////////////////////////////////////
//
// Minloc and Maxloc require slightly different functionality.
// They need to keep track of the location in addition to finding
// the min and max.
//
//////////////////////////////////////////////////////////////////////

#ifdef __MWERKS__
// Workaround for CodeWarrior 4 bug
// Funny "MMMin" name is to avoid conflict with Min defined elsewhere in r1,
// and "MMin defined iin IpplExpressions.h.
template<class T, unsigned D> 
struct MMMin {
  typedef typename T::PETE_Expr_t::PETE_Return_t type_t;
  static inline type_t apply(const PETE_Expr<T>& expr, NDIndex<D>& loc) {
    return Reduction(Expressionize<typename T::PETE_Expr_t>::apply( expr.PETE_unwrap().MakeExpression() ), 
	      OpLT(), OpMinAssign(),loc);
  }
};

template<class T, unsigned D> 
inline typename MMMin<T,D>::type_t 
min(const PETE_Expr<T>& expr, NDIndex<D>& loc) {
  return MMMin<T,D>::apply(expr, loc);
}
#else
template<class T, unsigned D> 
inline typename T::PETE_Expr_t::PETE_Return_t
min(const PETE_Expr<T>& expr, NDIndex<D>& loc)
{
  return Reduction(Expressionize<typename T::PETE_Expr_t>::apply(expr.PETE_unwrap().MakeExpression()),
		   OpLT(),OpMinAssign(),loc);
}
#endif // __MWERKS__

#ifdef __MWERKS__
// Workaround for CodeWarrior 4 bug
// Funny "MMMax" name is to avoid conflict with Max defined elsewhere in r1,
// and "MMax defined iin IpplExpressions.h.
template<class T, unsigned D> 
struct MMMax {
  typedef typename T::PETE_Expr_t::PETE_Return_t type_t;
  static inline type_t apply(const PETE_Expr<T>& expr, NDIndex<D>& loc) {
    return Reduction(Expressionize<typename T::PETE_Expr_t>::apply( expr.PETE_unwrap().MakeExpression() ), 
	      OpGT(), OpMaxAssign(),loc);
  }
};

template<class T, unsigned D> 
inline typename MMMax<T,D>::type_t 
max(const PETE_Expr<T>& expr, NDIndex<D>& loc) {
  return MMMax<T,D>::apply(expr, loc);
}
#else
template<class T, unsigned D>
inline typename T::PETE_Expr_t::PETE_Return_t
max(const PETE_Expr<T>& expr, NDIndex<D>& loc)
{
  return Reduction(Expressionize<typename T::PETE_Expr_t>::apply(expr.PETE_unwrap().MakeExpression()),
		   OpGT(),OpMaxAssign(),loc);
}
#endif // __MWERKS__

//////////////////////////////////////////////////////////////////////

#include "Field/ReductionLoc.hpp"

#endif // REDUCTION_LOC_H

/***************************************************************************
 * $RCSfile: ReductionLoc.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:26 $
 * IPPL_VERSION_ID: $Id: ReductionLoc.h,v 1.1.1.1 2003/01/23 07:40:26 adelmann Exp $ 
 ***************************************************************************/




