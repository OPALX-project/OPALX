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

// include files
#include "FieldLayout/ConejoBalancer.h"
#include "FieldLayout/MultiBalancer.h"
#include "Utility/IpplInfo.h"

//////////////////////////////////////////////////////////////////////
//
// PUBLIC FUNCTIONS
//
//////////////////////////////////////////////////////////////////////

//
// The constructor for ConejoBalancer.
// Record the number of processors in the system.
//

ConejoBalancer::ConejoBalancer()
:
  // Initialize the pointer to the balancer to zero
  // because we can't construct it yet.
  m_balancer(0), 
  // Initialize the number of vnodes to -1 because we don't know that yet.
  m_localVnodes(-1),
  m_totalVnodes(-1),
  // Record the number of processors.
  m_procs( IpplInfo::getNodes() ),
  // Record which processor is mine.
  m_myProc( IpplInfo::myNode() )
{
  // Record the number of vnodes on each processor.
  // Start the value out with -1 so we know if it has been set or not.
  m_vnodeCounts.resize(m_procs);
  for (int i=0; i<m_procs; ++i)
    m_vnodeCounts[i] = -1;
}

//////////////////////////////////////////////////////////////////////

//
// The destructor for ConejoBalancer.
// Just delete the balancer.
//

ConejoBalancer::~ConejoBalancer()
{
  delete m_balancer;
}

//////////////////////////////////////////////////////////////////////
//
// PRIVATE FUNCTIONS
//
//////////////////////////////////////////////////////////////////////

//
// ConejoBalancer::sendWeights
//
// INPUTS:
// vnodeWeights: A vector<double> of weights for vnodes on this proc.
// tag: an integer tag for the message.
//
// OUTPUTS:
// Message is sent.
//

void
ConejoBalancer::sendWeights(vector<double>& vnodeWeights, int tag)
{
  // Build a message to be sent.
  Message *msg = new Message();

  // Add data to the message.
  // First the number of elements, then the elements themselves.
  msg->put(vnodeWeights.size());
  putMessage(*msg, vnodeWeights.begin(), vnodeWeights.end());

  // Send the message.
  Ippl::Comm->send(msg, 0, tag);
}

//////////////////////////////////////////////////////////////////////

//
// ConejoBalancer::receiveWeights
//
// INPUTS:
// vnodeWeights: A vector<double>. Initially has the data from this
//               processor. Used as a place to put incoming data.
// tag: an integer tag for the messages.
//
// OUTPUTS:
// Updated state of the ConejoBalancer with the data from 
// the other processors.
//

void
ConejoBalancer::receiveWeights(vector<double>& vnodeWeights, int tag)
{
  // Tell the MultiBalancer that it is about to get a new material.
  m_balancer->newMaterial();

  // Receive messages from all the other processors.
  for ( int proc=0; proc < m_procs ; ++proc )
    {
      // Receive a message from the next proc.
      Message *msg = Ippl::Comm->receive_block(proc,tag);

      // Get the number of vnodes.
      long int count;
      msg->get(count);
      
      // Record the number of vnodes on that proc.
      recordVnodeCount(count,proc);

      // Allocate space for the incoming data.
      vnodeWeights.resize( count );

      // Unpack it.
      getMessage_iter(*msg, vnodeWeights.begin() );

      // Record it by extracting a double* from the vnodeWeights iterator.
      double *p = &*(vnodeWeights.begin());
      m_balancer->appendWeights(p,p+vnodeWeights.size());

      // Delete the message.
      delete msg;
    }
}

//////////////////////////////////////////////////////////////////////

//
// Tell the ConejoBalancer the number of vnodes it is dealing with.
//
// INPUT:
//    localVnodes:  The number of vnodes on this processor.
//    totalVnodes:  The total number of vnodes on all the processors.
//
// OUTPUT:
//    Updated state of the ConejoBalancer.
//    The number of vnodes is recorded.
//    The MultiBalancer is initialized.
//

void
ConejoBalancer::setupVnodes(int localVnodes, int remoteVnodes)
{
  // The first time through we need to record the sizes, 
  // and build the balancer
  if ( m_balancer == 0 )
    {
      // Record the sizes.
      m_localVnodes = localVnodes;
      m_totalVnodes = localVnodes + remoteVnodes;

      // If we are processor zero, build the balancer object.
      if ( m_myProc == 0 )
	m_balancer = new MultiBalancer(m_procs,m_totalVnodes);
    }
  // If it isn't the first time, make sure this BareField
  // is consistent with the previous.
  else
    {
      PAssert(m_localVnodes == localVnodes);
      PAssert(m_totalVnodes == localVnodes + remoteVnodes);
      PAssert( (m_balancer!=0) == (m_myProc==0 ) );
    }
}

//////////////////////////////////////////////////////////////////////

//
// Tell the ConejoBalancer the number of vnodes on a given processor.
//
// INPUT:
//   count: an integer for the number of vnodes
//   proc:  an integer of the processor
//
// OUTPUT:
//   Updated state of the ConejoBalancer.
//   The first time it is called it records the count.
//   After that, it checks to make sure the counts agree.
//

void 
ConejoBalancer::recordVnodeCount(int count, int proc)
{
  // Make sure this processor number makes sense.
  PAssert( proc >= 0 );

  // Check to see if this is the first time it is being called.
  if ( m_vnodeCounts[proc] < 0 )
    {
      // First time, record the count.
      m_vnodeCounts[proc] = count;
    }
  // Not the first time.
  else
    {
      // Make sure this count agrees with record.
      PAssert( m_vnodeCounts[proc] == count );
    }
}

//////////////////////////////////////////////////////////////////////

//
// ConejoBalancer::broadcastVnodesToSend
//
// Scan through the data in the MultiBalancer and figure out
// which processors will be sending vnodes to where.
// Broadcast that information to the processors that will
// be sending.
//
// INPUT:
//   tag: an integer message tag.
//
// OUTPUT:
//   Messages sent.
//

void 
ConejoBalancer::broadcastVnodesToSend(int tag)
{
  // Get a pointer to the destination processors for each vnode.
  MultiBalancer::iterator vp = m_balancer->begin();

  // Loop over the source processors, sending data to each.
  for ( int sourceProc = 0; sourceProc < m_procs; ++sourceProc )
    {
      // Cache the number of vnodes on this processor.
      long int c = m_vnodeCounts[sourceProc];

      // Build the message.
      Message *msg = new Message;
      msg->put(c);
      putMessage(*msg, vp, vp + c);

      // Send the message.
      Ippl::Comm->send(msg,sourceProc,tag);

      // Increment the pointer to the destination processors.
      vp += c;
    }
}

//////////////////////////////////////////////////////////////////////

//
// ConejoBalancer::receiveVnodesToSend
//
// Every processor receives from processor zero the destinations
// for its vnodes.
//
// INPUT:
//   tag: an integer message tag.
//
// OUTPUT:
//   vnodeDestinations: A vector<int>& to hold the destinations for 
//                      all the vnodes on this proc.
//

void 
ConejoBalancer::receiveVnodesToSend(vector<int>& vnodeDestinations, int tag)
{
  // Receive the message from proc 0.
  int proc_zero = 0;
  Message *msg = Ippl::Comm->receive_block(proc_zero,tag);

  // Get the number of vnodes coming in.
  long int s;
  msg->get(s);

  // Extract the vnode destinations.
  vnodeDestinations.resize(s);
  getMessage_iter(*msg, vnodeDestinations.begin() );

  // Delete the message.
  delete msg;
}

/***************************************************************************
 * $RCSfile: ConejoBalancer_inst.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:27 $
 * IPPL_VERSION_ID: $Id: ConejoBalancer_inst.cpp,v 1.1.1.1 2003/01/23 07:40:27 adelmann Exp $ 
 ***************************************************************************/
