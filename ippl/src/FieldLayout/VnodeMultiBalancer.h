// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef VNODE_MULTI_BALANCER_H
#define VNODE_MULTI_BALANCER_H

//-----------------------------------------------------------------------------
// Description:
// Vnode-granularity load balancer, based on input container of bool weight
// BareFields
//-----------------------------------------------------------------------------

// include files
#include <vector>

// Forward declarations:
template<unsigned Dim> class FieldLayout;
template<class T, unsigned Dim> class BareField;

//-----------------------------------------------------------------------------
// Full Description:

// A fairly simple load balancer inspired by Dan Quinlan's MLB, which balances
// at the level of granularity of the vnode, rather than the individual Field
// element (BinaryBalancer does the latter).

// It requires a FieldLayout expressed as a product of vnodes along each
// direction, meaning a power-of-two vnode count or else use of the special
// [Centered]FieldLayout constructors specifying vnodes per direction.

// Internally, it constructs a Field<...,double,....> having one element for
// each vnode along each direction in the input FieldLayout, called vf. It
// assigns a weight of 1.0 to each element equal to the total number of unique
// Field's having any true elements in each vnode for all the Field's in the
// input container (STL vector) of boolean weight Fields. If any vnode is
// compressed to false, the corresponding vf element value is left unchanged;
// if a vnode is uncompressed, this means that there is at least one true
// element in it, and the vf element value is incremented by one; if a vnode is
// compressed to true, this means all elements are true and the vf element is
// also incremented by 1.  Then, it invokes BinaryBalancer on this small Field
// of weights, vf, to partition it elementwise among the processors, one
// element per PE. Finally, it maps this back to the input FieldLayout (used by
// the input Fields of weights), and repartitions the input FieldLayout's
// vnodes in a corresponding way among the PE's.

// There is one function defined here:

// Calculate and apply a local domain for a binary repartition.
template<unsigned Dim>
void VnodeMultiRepartition(FieldLayout<Dim>& layout, 
			   std::vector<BareField<bool,Dim>* >& weights);

#include "FieldLayout/VnodeMultiBalancer.hpp"

#endif // VNODE_MULTI_BALANCER_H

/***************************************************************************
 * $RCSfile: VnodeMultiBalancer.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:27 $
 * IPPL_VERSION_ID: $Id: VnodeMultiBalancer.h,v 1.1.1.1 2003/01/23 07:40:27 adelmann Exp $ 
 ***************************************************************************/

