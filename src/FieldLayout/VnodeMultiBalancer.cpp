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
 * Visit www.amas.web.psi for more details
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

//-----------------------------------------------------------------------------
// Description:
// Vnode-granularity load balancer, based on input container of bool weight
// BareFields
//-----------------------------------------------------------------------------

// include files
#include "FieldLayout/VnodeMultiBalancer.h"
#include "FieldLayout/BinaryBalancer.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/BareField.h"
#include "Profile/Profiler.h"

///////////////////////////////////////////////////////////////////////////////
// Implementation of VnodeMultiRepartition().
// 
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
///////////////////////////////////////////////////////////////////////////////


template<unsigned Dim>
void VnodeMultiRepartition(FieldLayout<Dim>& layout, 
   			   std::vector<BareField<bool,Dim>* >& weights) {
  
  TAU_TYPE_STRING(taustr, " (" + CT(layout) + ", "+ CT(weights) + " )");
  TAU_PROFILE("VnodeMultiRepartition()", taustr, TAU_LAYOUT);

  int npe = Ippl::Comm->getNodes(); // Total number of PE's (pnodes)
  if (npe == 1) return; // Not much hope of balancing on 1 node, eh?

  // Get numbers of vnodes per direction (vnode "array" extents) and total
  // number of vnodes:
  unsigned vpd[Dim];
  int vnodes = 1;
  for (int d=0; d<Dim; d++) {
    vpd[d] = layout.getVnodesPerDirection(d);
    vnodes *= vpd[d];
  }

  // Construct the Field sized as vpd[d] elements along each dimension d:
  NDIndex<Dim> vfndi;
  for (int d=0; d<Dim; d++) vfndi[d] = Index(vpd[d]);
  // SERIAL/PARALLEL specification must match input FieldLayout layout:
  e_dim_tag edt[Dim];
  for (int d=0; d<Dim; d++) edt[d] = layout.getDistribution(d);
  // Because "true" for recurse parameter selected here, this algorithm will
  // run faster the first time through if the input "layout" was also
  // constructed with recurse=true. Should be correct even if not, though, and
  // subsequent re-calls of VnodeMultiRepartition() should go faster in either
  // case:
  FieldLayout<Dim> l(vfndi, edt, vpd, true, vnodes);
  BareField<double,Dim> vf(l);
  vf = 0.0;

  
  // Loop through the PE's owned LField's (<==> Vnode's) in each boolean Field
  // of weights. Check each for compression, and increment the value the value
  // in an array of double's with one element for each vnode; increment by 1 if
  // compressed to true or uncompressed; if compressed to false don't
  // increment:
  int* anyTrueElements = new int[vnodes];
  for (int v=0; v < vnodes; v++) anyTrueElements[v] = 0;
  
  // Iterate through the BareFields in the container:
  typename std::vector<BareField<bool,Dim>* >::iterator bfi;
  for (bfi = weights.begin(); bfi != weights.end(); ++bfi) {
    BareField<bool,Dim>& weight = *(*bfi);
    typename BareField<bool,Dim>::iterator_if weightItr;
    // Go through the LFields in the BareField:
    for (weightItr = weight.begin_if();
	 weightItr != weight.end_if(); 
	 ++weightItr) {
      if ((*weightItr).second->IsCompressed()) {
	if ((*weightItr).second->getCompressedData()) {
	  anyTrueElements[(*weightItr).second->getVnode()] += 1;
	}
      } else {
	anyTrueElements[(*weightItr).second->getVnode()] += 1;
      }
    }
  }


  // Now broadcast so that every PE has complete copy of "anyTrueElements"
  // array: Send to PE 0, which combines all, then broadcasts back:
  int pe = Ippl::Comm->myNode();   // My processor ID
  Message *msg;
  int partialTag  = Ippl::Comm->next_tag(VNMB_PARTIAL_TAG, VNMB_TAG_CYCLE);
  int completeTag = Ippl::Comm->next_tag(VNMB_COMPLETE_TAG, VNMB_TAG_CYCLE);
  if (pe == 0) {
    // Receive partially-filled arrays from other PE's:
    int* anyTrueElementsPartial = new int[vnodes];
    int notReceived = npe - 1;
    while (notReceived > 0) {
      int otherPE = COMM_ANY_NODE;
      int tag = partialTag;
      Message* msg2 = Ippl::Comm->receive_block(otherPE, tag);
      msg2->get(anyTrueElementsPartial);
      delete msg2;
      // Put values into anyTrueElements:
      for (int v=0; v < vnodes; v++) {
	// Any nonzero element in anyTrueElementsPartial *must* correspond to a
	// zero element in current status of anyTrueElements on PE 0; check for
	// otherwise and give error message if otherwise:
	if ((anyTrueElements[v] != 0) && (anyTrueElementsPartial[v] != 0)) {
	  ERRORMSG("VnodeMultiRepartition(): anyTrueElements[" << v << "] = "
		   << anyTrueElements[v] << " and anyTrueElementsPartial[" 
		   << v << "] = " << anyTrueElementsPartial[v] 
		   << " ; inconsistent!" << endl);
	}
	anyTrueElements[v] += anyTrueElementsPartial[v];
      }
      notReceived--;
    }
    msg = new Message();
    msg->put(anyTrueElements, anyTrueElements + vnodes);
    // Broadcast fully-filled array to other PE's:
    Ippl::Comm->broadcast_others(msg, completeTag);
    delete [] anyTrueElementsPartial;
  } else {
    // Send my partially-filled array to PE 0:
    msg = new Message();
    msg->put(anyTrueElements, anyTrueElements + vnodes);
    Ippl::Comm->send(msg, 0, partialTag);
    // Receive fully-filled array from PE 0:
    int pe0 = 0;
    msg = Ippl::Comm->receive_block(pe0, completeTag);
    msg->get(anyTrueElements);
  }

  // Loop through the PE's owned LField's (<==> Vnode's) in vf:
  typename BareField<double,Dim>::iterator_if vfItr;
  for (vfItr = vf.begin_if(); vfItr != vf.end_if(); ++vfItr) {
    
    // Global integer index of this vnode:
    int vnode = (*vfItr).second->getVnode();

    // Assign the integer value from anyTrueElements to the corresponding
    // single element in vf's corresponding vnode: The use of the special
    // LField::Compress(double val) function is a trick; we know that the vf
    // Field can be compressed (and probably always is); using this function
    // assigns the CompressedData value directly, or else compresses the LField
    // and sets CompressedData to the requested value if for some reason it
    // wasn't already compressed:
    (*vfItr).second->Compress(anyTrueElements[vnode]);

  }


  // Now invoke BinaryBalancer on the small Field vf:
  BinaryRepartition(l, vf);

  // Try this more rational alternative to handling zero-length vnodes:

  // Find and record all zero-length vnodes in l; each PE can check for this
  // independently, since all have domain size info for both local and remote
  // vnodes:
  int nZeroSizeVnodes = 0;
  // Locals:
  typename FieldLayout<Dim>::iterator_iv lLocals;
  for (lLocals = l.begin_iv(); lLocals != l.end_iv(); ++lLocals) {
    if ((*lLocals).second->getDomain().size() == 0) {nZeroSizeVnodes += 1;}
  }
  // Remotes:
  typename FieldLayout<Dim>::iterator_dv lRemotes;
  for (lRemotes = l.begin_rdv(); lRemotes != l.end_rdv(); ++lRemotes) {
    if ((*lRemotes).second->getDomain().size() == 0) {nZeroSizeVnodes += 1;}
  }

  // For now, punt on handling this and just return without changing the input
  // FieldLayout at all; report this as an error. The stuff ifdef'd by "GROSS"
  // below is partially-completed (and partially-compiling) kludge to try and
  // make BinaryRepartition do something in this case by tweaking the vf Field
  // values with random numbers.
  if (nZeroSizeVnodes != 0) {
    WARNMSG("VnodeMultiRepartition() was not able to get a successful "
	    << "rebalance. So, it is leaving the FieldLayout vnode "
	    << "PE-partioning the way it was whenyou called it. Sorry about "
	    << "that; there just aren't enough noncompressed vnodes to go "
	    << "around for this many PEs---at least not enough located in a "
	    << "way that the underlying BinaryBalancer algorithm can deal "
	    << "with." << endl);
    // Cleanup:
    delete [] anyTrueElements;
    return;
  }

#ifdef GROSS
  // Moved this gross stuff to the bottom of the file; hoist it back in here
  // later if needed. Of course, this module will not compile when GROSS is
  // defined until you do this hoisting, and fix remaining compile errors in
  // the gross stuff.
#endif // GROSS


  // This has changed the FieldLayout l so that it no longer has the same
  // number of vnodes as layout; it has one vnode per pnode.

  // Now must go through each *element* of vf, find which PE owns it now, map
  // the element index back to a vnode-array (and then global vnode ID) index
  // in layout, and reassign ownership of layout's corresponding vnode to that
  // same PE.

  int* peOwner = new int[vnodes];
  // Mask value -1 needed for inter-PE-exchange of array values below:
  for (int v=0; v < vnodes; v++) peOwner[v] = -1;

  // Loop through the elements of the Field vf. We know there is one per vnode
  // in the Field vf:

  // Outer loop over PE's owned Lfields in vf; vfItr constructed already:
  for (vfItr = vf.begin_if(); vfItr != vf.end_if(); ++vfItr) {

    // Inner loop over elements in LField:
    typename LField<double,Dim>::iterator lfi;
    for (lfi = (*vfItr).second->begin(); 
	 lfi != (*vfItr).second->end(); ++lfi) {
      // Global integer index values of the Field element this refers to:
      int vfIndex[Dim];
      // To compute this from result of Lfield::iterator::GetOffset(), must
      // have the global-index-space base index values of this LField's
      // subdomain so you can add it on:
      int lfBase[Dim];
      for (int d=0; d<Dim; d++) {
	lfBase[d] = (*vfItr).second->getOwned()[d].first();
      }
      for (int d=0; d<Dim; d++) vfIndex[d] = lfi.GetOffset(d) + lfBase[d];
      // Global integer index of this vnode:
      int vnode = vfIndex[0];
      int multipplier = 1;
      for (int d=1; d<Dim; d++) {
	multipplier *= vpd[d-1];
	vnode += vfIndex[d]*multipplier;
      }
      if (vnode >= vnodes) {
	ERRORMSG("VnodeMultiRepartition(): vnode = " << vnode
		 << " but vnodes is only " << vnodes << " ; inconsistent!"
		 << endl);
	PInsist(vnode < vnodes, 
		"VnodeMultiRepartition: exit because of vnode value error.");
      }

      // Record PE ownership. This is SPMD code here. The PE that calculated
      // the value of vnode is the one that owns it, so assign ownership to
      // it's PE ID number:
      peOwner[vnode] = pe;
    }
  }


  // Now broadcast so that every PE has complete copy of "peOwner" array:
  partialTag  = Ippl::Comm->next_tag(VNMB_PARTIAL_TAG, VNMB_TAG_CYCLE);
  completeTag = Ippl::Comm->next_tag(VNMB_COMPLETE_TAG, VNMB_TAG_CYCLE);
  Message* msg4;
  if (pe == 0) {
    // Receive partially-filled arrays from other PE's:
    int* peOwnerPartial = new int[vnodes];
    int notReceived = npe - 1;
    while (notReceived > 0) {
      int otherPE = COMM_ANY_NODE;
      int tag = partialTag;
      Message* msg2 = Ippl::Comm->receive_block(otherPE, tag);
      msg2->getmsg((void *)peOwnerPartial);
      delete msg2;
      // Put values into peOwner:
      for (int v=0; v < vnodes; v++) {
	if (peOwnerPartial[v] != -1) {
	  // Any non-minus-one element in peOwnerPartial *must* correspond to a
	  // minus-one element in current status of peOwner on PE 0; check for
	  // otherwise and give error message if otherwise:
	  if (peOwner[v] != -1) {
	    ERRORMSG("VnodeMultiRepartition(): peOwner[" << v << "] = "
		     << peOwner[v] << " and peOwnerPartial[" 
		     << v << "] = " << peOwnerPartial[v] 
		     << " ; inconsistent!" << endl);
	  }
	  peOwner[v] = peOwnerPartial[v];
	}
      }
      notReceived--;
    }
    msg4 = new Message();
    msg4->put(peOwner, peOwner + vnodes);
    // Broadcast fully-filled array to other PE's:
    Ippl::Comm->broadcast_others(msg4, completeTag);
    delete [] peOwnerPartial;
  } else {
    // Send my partially-filled array to PE 0:
    msg4 = new Message();
    msg4->put(peOwner, peOwner + vnodes);
    Ippl::Comm->send(msg4, 0, partialTag);
    // Receive fully-filled array from PE 0:
    int pe0 = 0;
    msg4 = Ippl::Comm->receive_block(pe0, completeTag);
    msg4->get(peOwner);
  }
  delete msg4;


  // Now repartition layout to have the same PE ownership as l; must construct
  // an array of NDIndex's before FieldLayout::Repartition can be invoked:
  
  // Find out how many vnodes I (PE) own:
  int nVnodesIOwn = 0;
  for (int v=0; v < vnodes; v++) if (peOwner[v] == pe) ++nVnodesIOwn;

  // Array of Vnodes that I own:
  Vnode<Dim>* domains = new Vnode<Dim>[nVnodesIOwn];

  // Get the values of the domains from the original layout:

  // The ones I owned in the original layout (locals):
  typename FieldLayout<Dim>::iterator_iv oldLocals;
  int domain = 0; // counter
  for (oldLocals = layout.begin_iv(); oldLocals != layout.end_iv(); 
       ++oldLocals) {
    // Global integer index of this vnode:
    int vnode = (*oldLocals).second->getVnode();
    if (peOwner[vnode] == pe) {
      domains[domain] = 
	Vnode<Dim>((*oldLocals).second->getDomain(), pe, vnode);
      ++domain;
    }
  }

  // The ones I didn't own in the original layout (remotes):
  typename FieldLayout<Dim>::iterator_dv oldRemotes;
  for (oldRemotes = layout.begin_rdv(); oldRemotes != layout.end_rdv(); 
       ++oldRemotes) {
    // Global integer index of this vnode:
    int vnode = (*oldRemotes).second->getVnode();
    if (peOwner[vnode] == pe) {
      domains[domain] = 
	Vnode<Dim>((*oldRemotes).second->getDomain(), pe, vnode);
      ++domain;
    }
  }


  // Finally, call FieldLayout::Repartition() to repartition layout:
  layout.Repartition(domains, domains + nVnodesIOwn);


  // Cleanup:
  delete [] anyTrueElements;
  delete [] peOwner;
  delete [] domains;

  return;
}

#ifdef GROSS
  // See comments above in first "GROSS" block about putting this stuff back up there.

  // If none are zero, go on from here. If any are zero, have to go back and
  // find them, and find suitable nonzero ones to divide to get one nonzero one
  // for each PE that owns a zero-size one:
  FieldLayout<Dim> l2(vfndi, edt, vpd, true, vnodes); // Outside if block for scoping.
  if (nZeroSizeVnodes != 0) {
    // Reconstruct a replacement vf Field, add a random epsilon fudge factor to
    // it, and see if the repartitioner works; if not, *don't* repartition and
    // just return (load balance attempt fails, so leave load as it was):

    // The replacement vf Field:
    //outside if block; see above    FieldLayout<Dim> l2(vfndi, edt, vpd, true, vnodes);
    BareField<double,Dim> vf2(l2);
    vf2 = 0.0;

    // -----------------------------------------------------------------------------------
    // Duplicate the original assignment of the vf Field:

    // Loop through the PE's owned LField's (<==> Vnode's) in each boolean Field
    // of weights. Check each for compression, and increment the value the value
    // in an array of double's with one element for each vnode; increment by 1 if
    // compressed to true or uncompressed; if compressed to false don't
    // increment:
    for (int v=0; v < vnodes; v++) anyTrueElements[v] = 0;
  
    // Iterate through the BareFields in the container:
    std::vector<BareField<bool,Dim>* >::iterator bfi2;
    for (bfi2 = weights.begin(); bfi2 != weights.end(); ++bfi2) {
      BareField<bool,Dim>& weight2 = *(*bfi2);
      BareField<bool,Dim>::iterator_if weight2Itr;
      // Go through the LFields in the BareField:
      for (weight2Itr = weight2.begin_if();
	   weight2Itr != weight2.end_if(); 
	   ++weight2Itr) {
	if ((*weight2Itr).second->IsCompressed()) {
	  if ((*weight2Itr).second->getCompressedData()) {
	    anyTrueElements[(*weight2Itr).second->getVnode()] += 1;
	  }
	} else {
	  anyTrueElements[(*weight2Itr).second->getVnode()] += 1;
	}
      }
    }
    // -----------------------------------------------------------------------------------

    // Now add the epsilon tweak:
    double normfact = sum(Abs(vf2));
    normfact = normfact/vfndi.size();
    double epsilon = 1.0/normfact;
    vf2 += epsilon*IpplRandom;

    // ***********************************************************************
    // Once again go through the code up to the attempted BinaryRepartition():
    // This should maybe a function....

    // Now broadcast so that every PE has complete copy of "anyTrueElements"
    // array: Send to PE 0, which combines all, then broadcasts back:
    partialTag  = Ippl::Comm->next_tag(VNMB_PARTIAL_TAG, VNMB_TAG_CYCLE);
    completeTag = Ippl::Comm->next_tag(VNMB_COMPLETE_TAG, VNMB_TAG_CYCLE);
    Message *msg3;
    if (pe == 0) {
      // Receive partially-filled arrays from other PE's:
      int* anyTrueElementsPartial = new int[vnodes];
      int notReceived = npe - 1;
      while (notReceived > 0) {
	int otherPE = COMM_ANY_NODE;
	int tag = partialTag;
	Message* msg2 = Ippl::Comm->receive_block(otherPE, tag);
	msg2->get(anyTrueElementsPartial);
	delete msg2;
	// Put values into anyTrueElements:
	for (int v=0; v < vnodes; v++) {
	  // Any nonzero element in anyTrueElementsPartial *must* correspond to a
	  // zero element in current status of anyTrueElements on PE 0; check for
	  // otherwise and give error message if otherwise:
	  if ((anyTrueElements[v] != 0) && (anyTrueElementsPartial[v] != 0)) {
	    ERRORMSG("VnodeMultiRepartition(): anyTrueElements[" << v << "] = "
		     << anyTrueElements[v] << " and anyTrueElementsPartial[" 
		     << v << "] = " << anyTrueElementsPartial[v] 
		     << " ; inconsistent!" << endl);
	  }
	  anyTrueElements[v] += anyTrueElementsPartial[v];
	}
	notReceived--;
      }
      msg3 = new Message();
      msg3->put(anyTrueElements, anyTrueElements + vnodes);
      // Broadcast fully-filled array to other PE's:
      Ippl::Comm->broadcast_others(msg3, completeTag);
      delete [] anyTrueElementsPartial;
    } else {
      // Send my partially-filled array to PE 0:
      msg3 = new Message();
      msg3->put(anyTrueElements, anyTrueElements + vnodes);
      Ippl::Comm->send(msg3, 0, partialTag);
      // Receive fully-filled array from PE 0:
      int pe0 = 0;
      msg3 = Ippl::Comm->receive_block(pe0, completeTag);
      msg3->get(anyTrueElements);
    }
    delete msg3;

    // Loop through the PE's owned LField's (<==> Vnode's) in vf2:
    BareField<double,Dim>::iterator_if vf2Itr;
    for (vf2Itr = vf2.begin_if(); vf2Itr != vf2.end_if(); ++vf2Itr) {
    
      // Global integer index of this vnode:
      int vnode = (*vf2Itr).second->getVnode();

      // Assign the integer value from anyTrueElements to the corresponding
      // single element in vf2's corresponding vnode: The use of the special
      // LField::Compress(double val) function is a trick; we know that the vf2
      // Field can be compressed (and probably always is); using this function
      // assigns the CompressedData value directly, or else compresses the LField
      // and sets CompressedData to the requested value if for some reason it
      // wasn't already compressed:
      (*vf2Itr).second->Compress(anyTrueElements[vnode]);

    }

    // Now invoke BinaryBalancer on the small Field vf2:
    BinaryRepartition(l2, vf2);

    // ***********************************************************************

    // Check once again for boo-boos (zero-length vnodes); if found, abandon balance:
    // Find and record all zero-length vnodes in l; each PE can check for this
    // independently, since all have domain size info for both local and remote
    // vnodes:
    int nZeroSizeVnodes2 = 0;
    // Locals:
    FieldLayout<Dim>::iterator_iv l2Locals;
    for (l2Locals = l2.begin_iv(); l2Locals != l2.end_iv(); ++l2Locals) {
      if ((*l2Locals).second->getDomain().size() == 0) {nZeroSizeVnodes2 += 1;}
    }
    // Remotes:
    FieldLayout<Dim>::iterator_dv l2Remotes;
    for (l2Remotes = l2.begin_rdv(); l2Remotes != l2.end_rdv(); ++l2Remotes) {
      if ((*l2Remotes).second->getDomain().size() == 0) {nZeroSizeVnodes2 += 1;}
    }
    // If none are zero, go on from here. If any are zero, have to go back and
    // find them, and find suitable nonzero ones to divide to get one nonzero one
    // for each PE that owns a zero-size one:
    if (nZeroSizeVnodes2 != 0) {
      WARNMSG("VnodeMultiRepartition(): even on a desperate 2nd attempt by adding in some"
	      << "random nonzero vnodes, was not able to get a successful rebalance. So, "
	      << "leaving the FieldLayout partioning the way it was when you called"
	      << " VnodeMultiRepartition(). Sorry about that." << endl);
      return;
    } else {
      // Success! (Of some sort....); repartition vf and assign it to vf2. Do
      // this by setting l equal to l2 and doing l.Repartition():
      l = l2;
      l2.Repartition
      vf = vf2;
    }
  }
#endif // GROSS

/***************************************************************************
 * $RCSfile: VnodeMultiBalancer.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:27 $
 * IPPL_VERSION_ID: $Id: VnodeMultiBalancer.cpp,v 1.1.1.1 2003/01/23 07:40:27 adelmann Exp $ 
 ***************************************************************************/
