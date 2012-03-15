// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by PSI. 
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "FieldLayout/ConejoBalancer.h"
#include "FieldLayout/MultiBalancer.h"
#include "Utility/IplInfo.h"

//////////////////////////////////////////////////////////////////////

//
// reduceLocalWeights
//
// Input: 
//   BareField<T,D> with weights
//   bool dropCompressed, which is true if compressed vnodes in weights
//        should have the weight of a single element instead of multiplied
//        by the number of elements.
// Outout: a container of doubles.
//
// Given a BareField and a container,
// loop over each vnode, finding the maximum weight
// on each vnode.
// Put those reduced weights in the container.
//
// This is a bare function, since it does not depend on anything
// in the class ConejoBalancer.
// This is an inline function so it does not generate an additional
// prelink step.
//

template<class T, unsigned int D>
inline void
reduceLocalWeights(BareField<T,D>& weights, 
                   vector<double>& vnodeWeights,
		   bool dropCompressed)
{
  // Get an iterator and loop over the LFields of the BareField.
  typename BareField<T,D>::iterator_if lf = weights.begin_if(); 
  for ( ; lf != weights.end_if() ; ++lf )
    {
      // Get an iterator and loop over the contents of this LField.
      typename LField<T,D>::iterator lp = (*lf).second->begin(); 

      // A place to record the total weight.
      // Start with just the first element.
      double x = *lp;

      // If the LField is compressed and and dropCompressed is true,
      // then the weight from this vnode comes from just the one element.
      if ( !(dropCompressed && (*lf).second->IsCompressed()) )
	{
	  // Add up all the values in the vnode.
	  for ( ++lp ; lp != (*lf).second->end() ; ++lp )
	    x += *lp;
	}

      // Append the largest value to the container.
      vnodeWeights.push_back( x );
    }
}

//////////////////////////////////////////////////////////////////////
//
// PRIVATE member functions
//
//////////////////////////////////////////////////////////////////////
//mwerks Moved into class definition (.h file).


/***************************************************************************
 * $RCSfile: ConejoBalancer.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:27 $
 * IPPL_VERSION_ID: $Id: ConejoBalancer.cpp,v 1.1.1.1 2003/01/23 07:40:27 adelmann Exp $ 
 ***************************************************************************/
