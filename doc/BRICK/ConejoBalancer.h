// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

//////////////////////////////////////////////////////////////////////
// Class ConejoBalancer
//////////////////////////////////////////////////////////////////////

#ifndef CONEJO_BALANCER_H
#define CONEJO_BALANCER_H

//////////////////////////////////////////////////////////////////////
//
// Description:
// A load balancer designed for the Conejo code.
// It distributes the vnodes so that each of the materials is 
// simultaneously load balanced.
//
//////////////////////////////////////////////////////////////////////

// include files
#include "Field/BareField.h"

#ifdef IPPL_STDSTL
#include <vector>
using std::vector;
#else
#include <vector.h>
#endif

// forward declarations
class MultiBalancer;
template<unsigned int D> class FieldLayout;
template<unsigned int D> class NDIndex;

/*

  ConejoBalancer is an interface class for MultiBalancer.

  **************************************************
  GENERAL DESCRIPTION
  **************************************************

  ConejoBalancer does the following things:

  1. Inputs a series of BareFields with weights for each
     grid location for each material.

  2. Sum the weights to a single value per vnode.

  3. Uses MultiBalancer to find the new distribution of vnodes.

  4. Rebalances one or more FieldLayouts for that new distribution.

 */

class ConejoBalancer
{
public:

  // Construct the ConejoBalancer.
  // This needs no arguments since it gets its information
  // from member functions below.
  ConejoBalancer();

  // Destroy the balancer.
  ~ConejoBalancer();

  // Add a new set of weights for a new material.
  // if dropCompressed is true, complressed vnodes in weights
  // are considered to have weight of just one element.
  //mwerks  template<class T, unsigned int D>
  //mwerks  void addMaterial(BareField<T,D>& weights, bool dropCompressed=false);
  //////////////////////////////////////////////////////////////////////

  //
  // ConejoBalancer::addMaterial
  //
  // Input: 
  //   BareField with weights for this material.
  //   bool dropCompressed, which is true if compressed vnodes in weights
  //        should have the weight of a single element instead of multiplied
  //        by the number of elements.
  // Output: The state of the ConejoBalancer includes the new weights.
  //
  // Extracting the weights goes through several phases:
  // 1. If this is the first time this has been called, 
  //    initialize the MultiBalancer.
  // 2. Calculate the weights for the local part of the BareField.
  // 3a. If this is not processor zero, send those weights to processor zero.
  // 3b. If it is processor zero, collect the weights.
  //

  template<class T, unsigned int D>
  void addMaterial(BareField<T,D>& weights, bool dropCompressed)
  {
    // Initialize (or check consistency).
    setupVnodes(weights.size_if(), weights.getLayout().size_rdv());

    // A container to hold the reduced weights for the local vnodes.
    vector<double> vnodeWeights;
    vnodeWeights.reserve(m_localVnodes);

    // Get the local weights.
    reduceLocalWeights(weights,vnodeWeights,dropCompressed);

    // Get a message tag.
    int tag = Ipl::Comm->next_tag(F_CONEJO_BALANCER_TAG, F_TAG_CYCLE);

    // Everybody sends their data to processor zero.
    sendWeights(vnodeWeights,tag);

    // If we are processor zero, process messages.
    if ( m_myProc == 0 )
      receiveWeights(vnodeWeights,tag);
  }


  // Redistribute a FieldLayout using the stored weights.
  //mwerks  template<unsigned int D>
  //mwerks  void redistribute(FieldLayout<D>& layout);
  //////////////////////////////////////////////////////////////////////

  //
  // ConejoBalancer::redistribute
  //
  // Redistribute a FieldLayout using the stored weights.
  //

  template<unsigned int D>
  void redistribute(FieldLayout<D>& layout)
  {
    // Get message tags.
    int bcasttag = Ipl::Comm->next_tag(F_CB_BCAST_TAG, F_TAG_CYCLE);
    int domaintag = Ipl::Comm->next_tag(F_CB_DOMAIN_TAG, F_TAG_CYCLE);

    // On proc 0, figure things out and send them.
    if ( m_myProc == 0 )
      {
	// Tell the MultiBalancer to figure out the new distribution.
	m_balancer->distribute();

      // Broadcast the vnode ids that each processor will have to send.
	broadcastVnodesToSend(bcasttag);
      }

    // Everywhere receives the id's of the vnodes it will send.
    vector<int> vnodeDestinations;
    receiveVnodesToSend(vnodeDestinations,bcasttag);

    // Send the domains for the vnodes to their new homes.
    sendVnodeDomains(vnodeDestinations,layout,domaintag);

    // Receive the domains for the vnodes that will live here.
    vector< NDIndex<D> > vnodeDomains;
    receiveVnodeDomains(vnodeDomains,domaintag);

  // Redistribute the FieldLayout using the new local domains.
    NDIndex<D> *p = &*(vnodeDomains.begin());
    layout.Repartition( p , p + vnodeDomains.size() );
  }


private:

  // Keep a pointer to the object that encapsulates the algorithm.
  MultiBalancer *m_balancer;

  // Remember the number of local vnodes.
  int m_localVnodes;

  // Remember the total number of vnodes.
  int m_totalVnodes;

  // Remember the number of processors.
  int m_procs;

  // Remember my processor.
  int m_myProc;

  // Record the number of vnodes on each processor.
  vector<int> m_vnodeCounts;

  // Support functions for internal use.
  void sendWeights(vector<double>& vnodeWeights, int tag);
  void receiveWeights(vector<double>& vnodeWeights, int tag);
  void setupVnodes(int localVnodes, int remoteVnodes);
  void recordVnodeCount(int count, int proc);
  void broadcastVnodesToSend(int tag);
  void receiveVnodesToSend(vector<int>& vnodeDestinations,int tag);

  //mwerks  template<unsigned int D>
  //mwerks  void sendVnodeDomains(vector<int>& vnodeDests,
  //mwerks			FieldLayout<D>& layout,
  //mwerks			int tag);
  //
  // ConejoBalancer::sendVnodeDomains
  //
  // Send to all the other processors the domains for the vnodes 
  // that will be sent.
  // Here just the NDIndexes are being sent.
  // The contents of the Fields on those domains will be sent later.
  //

  template<unsigned int D>
  void sendVnodeDomains(vector<int>& vnodeDests, 
			FieldLayout<D>& layout,
			int tag)
  {
    // A buffer for writing down the domains to be sent.
    vector< NDIndex<D> > send;

  // Loop over processors, sending a message to each.
    for ( int proc = 0; proc < m_procs; ++proc )
      {
	// Loop over the local vnodes, figuring out where each should go.
	vector<int>::iterator vp = vnodeDests.begin();
	typename FieldLayout<D>::iterator_iv fp = layout.begin_iv();
	for ( ; vp!=vnodeDests.end(); ++vp, ++fp)
	  {
	    // If this vnode is going to processor proc
	    if ( *vp == proc ) 
	      // Record the domain for sending.
	      send.push_back( (*fp).second->getDomain() );
	  }

	// Build the message to be sent.
	Message *msg = new Message;

	// Add the domains to the message.
	NDIndex<D> *p = &*(send.begin());
	msg->put(send.size());
	putMessage(*msg, p , p + send.size() );

	// Send the message.
	Ipl::Comm->send(msg,proc,tag);

	// Clear the send container.
	send.clear();
      }
  }

  //mwerks  template<unsigned int D>
  //mwerks  void receiveVnodeDomains(vector< NDIndex<D> >& vnodeDomains, int tag);
  //////////////////////////////////////////////////////////////////////
  //
  // ConejoBalancer::receiveVnodeDomains
  //
  // Each processor receives from all the other processors
  // the domains it will have after the redistribution.
  //
  template<unsigned int D>
  void receiveVnodeDomains(vector< NDIndex<D> >& vnodeDomains, 
			   int tag)
  {
    // Loop over all the processors, receiving a message from each.
    for (int proc=0; proc < m_procs; ++proc)
      {
	// Receive a message from any processor.
	int any_proc = -1;
	Message *msg = Ipl::Comm->receive_block(any_proc,tag);

	// Get the number of NDIndexes in this message.
	long int s;
	msg->get(s);

      // Make sure the size isn't negative.
	PAssert(s>=0);

	// If there are any there, unpack them.
	if ( s != 0 ) 
	  {
	    // Add room to the container.
	    vnodeDomains.resize( vnodeDomains.size() + s );

          // Unpack onto the end of the container.
	    getMessage_iter(*msg, vnodeDomains.end()-s );
	  }

	// Delete the message.
	delete msg;
      }
  }
};

#include "FieldLayout/ConejoBalancer.cpp"

#endif // CONEJO_BALANCER_H

/***************************************************************************
 * $RCSfile: ConejoBalancer.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:27 $
 * IPPL_VERSION_ID: $Id: ConejoBalancer.h,v 1.1.1.1 2003/01/23 07:40:27 adelmann Exp $ 
 ***************************************************************************/
