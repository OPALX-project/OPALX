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
#include "DataSource/ACLVISPtclBaseDataSource.h"
#include "DataSource/ACLVISOperations.h"
#include "DataSource/ACLVISDataConnect.h"
#include "DataSource/PtclAttribDataSource.h"
#include "Particle/IpplParticleBase.h"
#include "Message/Communicate.h"
#include "Utility/IpplInfo.h"
 



////////////////////////////////////////////////////////////////////////////
// constructor: name, connection method, transfer method, pbase
template<class PLayout>
ACLVISIpplParticleBaseDataSource<PLayout>::ACLVISIpplParticleBaseDataSource(const
								    char *nm,
         DataConnect *dc, int tm, IpplParticleBase<PLayout>& PB)
  : IpplParticleBaseDataSource(nm, dc, tm, &PB), MyIpplParticleBase(PB),
    IDMap(&IDMapA), NewIDMap(&IDMapB)
{

  // make sure we do not have a IpplParticleBase with positions greater than
  // 3 dimensions
  CTAssert(PLayout::Dimension < 4);

  // then check to make sure we're not trying to make a connection that we
  // don't support in this class
  std::string filestring = "aclvis";
  if (std::string(dc->DSID()) != filestring) {
    ERRORMSG("Illegal DataConnect object for ACLVIS Data Object." << endl);
    ACLVISConnection = 0;
    Connection = 0;
  } else if (tm != DataSource::OUTPUT) {
    ERRORMSG("ACLVIS data connections may only be of type OUTPUT." << endl);
    ACLVISConnection = 0;
    Connection = 0;
  } else {
    ACLVISConnection = dynamic_cast<ACLVISDataConnect *>(dc);
    checkin();
  }
}


////////////////////////////////////////////////////////////////////////////
// destructor
template<class PLayout>
ACLVISIpplParticleBaseDataSource<PLayout>::~ACLVISIpplParticleBaseDataSource() {

}


////////////////////////////////////////////////////////////////////////////
// make a connection using the given attribute.  Return success.
template<class PLayout>
bool ACLVISIpplParticleBaseDataSource<PLayout>::connect_attrib(
	   ParticleAttribDataSource *pa) {

  // on parent node, establish connection to ACLVIS for this base+attrib pair
  if (getConnection()->onConnectNode() && connected())
    // call ACLVIS API to register this data object
    ACLVISConnection->getConnection()->connect((void *)pa,
		       (char *)(pa->name()),
		       (ReadParticleTool *)(pa->getConnectStorage()),
		       ParticleDataType);

  return IpplParticleBaseDataSource::connect_attrib(pa);
}


////////////////////////////////////////////////////////////////////////////
// disconnect from the external agency the connection involving this
// particle base and the given attribute.  Return success.
template<class PLayout>
bool ACLVISIpplParticleBaseDataSource<PLayout>::disconnect_attrib(
	   ParticleAttribDataSource *pa) {

  // on parent node, remove connection to ACLVIS for this base+attrib pair
  if (getConnection()->onConnectNode() && connected())
#ifdef IPPL_LUX
    ACLVISConnection->getConnection()->disconnect((char *)(pa->name()));
#else
    ACLVISConnection->getConnection()->disconnect((void *)pa);
#endif

  return IpplParticleBaseDataSource::disconnect_attrib(pa);
}


////////////////////////////////////////////////////////////////////////////
// check to see if the given ParticleAttrib is in this IpplParticleBase's
// list of registered attributes.  Return true if this is so.
template<class PLayout>
bool ACLVISIpplParticleBaseDataSource<PLayout>::has_attrib(
	    ParticleAttribBase *pa) {

  // go through the list of registered attributes in our IpplParticleBase ...
  typename IpplParticleBase<PLayout>::attrib_iterator attr=MyIpplParticleBase.begin();
  typename IpplParticleBase<PLayout>::attrib_iterator endattr=MyIpplParticleBase.end();
  for ( ; attr != endattr; ++attr)
    if (pa == *attr)
      return true;

  // if here, not found
  return false;
}


////////////////////////////////////////////////////////////////////////////
// Indicate to the receiver that we're allowing them time to manipulate the
// data (e.g., for a viz program, to rotate it, change representation, etc.)
// This should only return when the manipulation is done.
template<class PLayout>
void ACLVISIpplParticleBaseDataSource<PLayout>::interact(const char *str) {
  
  

  // on parent node, hand off control to ACLVIS API.  If a command string
  // is given, instead call the viz API function to execute it as an
  // interactive command.
  if (getConnection()->onConnectNode() && connected()) {
    if (str != 0 && *str != 0) {
#ifndef IPPL_LUX
      ACLVISConnection->getConnection()->InterpretCommand(str);
#endif
    } else {
#ifdef IPPL_LUX
      ACLVISConnection->getConnection()->Interact();
#else
      ACLVISConnection->getConnection()->Interact(1);
#endif
    }
  }
}


////////////////////////////////////////////////////////////////////////////
// Update the object, that is, make sure the receiver of the data has a
// current and consistent snapshot of the current state.  Return success.
template<class PLayout>
bool ACLVISIpplParticleBaseDataSource<PLayout>::update() {
  
  

  //Inform dbgmsg("APB::update", INFORM_ALL_NODES);
  //dbgmsg << "Starting update on node " << Ippl::myNode() << endl;

  // only need to do update if we have connected attribs
  if (begin_attrib() == end_attrib())
    return true;

  // get the next tag for communication
  int n, tag = Ippl::Comm->next_tag(DS_PB_TAG, DS_CYCLE);
  unsigned N = Ippl::getNodes();
  unsigned myN = Ippl::myNode();
  unsigned locN = MyIpplParticleBase.getLocalNum();
  unsigned totN = MyIpplParticleBase.getTotalNum();
  unsigned currIndx = 0;

  typename AttribList_t::iterator attr;
  typename AttribList_t::iterator endattr;

  // clear out the id list
  NewIDMap->erase(NewIDMap->begin(), NewIDMap->end());

  // when running in parallel, we must send all data to parent nodes
  //dbgmsg << "About to send all the messages out ..." << endl;
  for (n=0; n < getConnection()->getNodes(); ++n) {
    if (myN != n) {
      // put number of local particles into message
      Message *msg = new Message;
      msg->put(locN);

      if (locN > 0) {
	// put position and ID data into a message
	MyIpplParticleBase.R.putMessage(*msg, locN, 0);
	MyIpplParticleBase.ID.putMessage(*msg, locN, 0);

	// put attribute data into a message
	for (attr=begin_attrib(), endattr=end_attrib();
	     attr != endattr;
	     ++attr)
	  (*attr)->putMessage(msg);
      }

      // send the data to parent
      //dbgmsg << "Sending " << locN << " particles to node " << n << endl;
      Ippl::Comm->send(msg, n, tag);
    }
  }

  // on the parent nodes, put in our own data and the data from other
  // nodes
  //dbgmsg << "About to receive in all the messages ..." << endl;
  if (getConnection()->onConnectNode()) {

    // on master nodes, prepare all the attributes for new data
    for (attr=begin_attrib(),endattr=end_attrib(); attr!=endattr; ++attr)
      (*attr)->prepare_data(true, totN);

    // process a message or local data from each node
    for (n=0; n < N; ++n) {
      if (myN == n) {
	// put in local data
	if (locN > 0) {
	  //dbgmsg << "Inserting local data ..." << endl;
	  for (attr=begin_attrib(),endattr=end_attrib();attr!=endattr;++attr){
	    (*attr)->insert_data(locN, currIndx, 0);
	    insert_pos(locN, currIndx, &(MyIpplParticleBase.R[0]),
		       &(MyIpplParticleBase.ID[0]), *attr);
	  }
	}

	// increment index, which indicates where the next particles go
	currIndx += locN;

      } else {
	// now receive message from current other node
	Message *msg = Ippl::Comm->receive_block(n, tag);
	unsigned remN;
	msg->get(remN);

	// extract attrib data into agency-specific data structure
	if (remN > 0) {
	  // get coordinate attribute from message
	  typename PLayout::SingleParticlePos_t* posdata =
	    (typename PLayout::SingleParticlePos_t *)(msg->remove());
	  typename PLayout::Index_t* iddata =
	    (typename PLayout::Index_t *)(msg->remove());

	  // get non-position attribute data from message, and put all the
	  // data into the ACLVIS data structures
	  //dbgmsg<<"Inserting remote data from node " << n << " ..." << endl;
	  for (attr=begin_attrib(),endattr=end_attrib();attr!=endattr;++attr){
	    (*attr)->insert_data(remN, currIndx, msg);
	    insert_pos(remN, currIndx, posdata, iddata, *attr);
	  }

	  // free the memory used to store position and ID values
	  free(static_cast<void *>(posdata));
	  free(static_cast<void *>(iddata));

	  // increment index, which indicates where the next particles go
	  currIndx += remN;
	}

	// we're done with this node
	delete msg;
      }
    }

    // tell ACLVIS API there is new data
    for (attr=begin_attrib(),endattr=end_attrib(); attr!=endattr; ++attr) {
      // on parent node, call ACLVIS API to update data for each attrib.
      (*attr)->prepare_data(false, totN);
#ifdef IPPL_LUX
      ACLVISConnection->getConnection()->update((char *)((*attr)->name()));
#else
      ACLVISConnection->getConnection()->update((void *)(*attr));
#endif
    }
  }

  // swap ID lists
  IDMap_t *tmpmap = IDMap;
  IDMap = NewIDMap;
  NewIDMap = tmpmap;

  // report success
  return true;
}


////////////////////////////////////////////////////////////////////////////
// copy data out of the given coordinates iterator into the data structure
// for the given ParticleAttribDataSource which represents the agency-
// specific storage.  Arguments = number of particles, starting index,
// pointer to value, pointer to ID's, and ParticleAttribDataSource
template<class PLayout>
void ACLVISIpplParticleBaseDataSource<PLayout>::insert_pos(
		      unsigned N, unsigned sIndx,
		      typename PLayout::SingleParticlePos_t* data,
		      typename PLayout::Index_t* iddata,
		      ParticleAttribDataSource* pa) { 

  //  Inform dbgmsg("APB::insert_pos", INFORM_ALL_NODES);

  // get ACLVIS data storage location for particle coordinates
  ReadParticleTool *vizdata =
    static_cast<ReadParticleTool *>(pa->getConnectStorage());

  // copy in the data from the storage
  unsigned maxIndx = sIndx + N;
  for ( ; sIndx < maxIndx; ++sIndx, ++data, ++iddata) {
    // set position for the Nth particle
    ACLVISTraits<ReadParticleTool, typename PLayout::SingleParticlePos_t>::
      setCoord(vizdata, sIndx, *data);

    // find where this particle was located during the last update
    typename IDMap_t::iterator idloc = IDMap->find(*iddata);
    int prevloc = (-1);
    if (idloc != IDMap->end())
      prevloc = (*idloc).second;

    // set ID information for the Nth particle
    ACLVISTraits<ReadParticleTool, typename PLayout::Index_t>::
      setID(vizdata, sIndx, prevloc);

    // save where this item is located in the list
    (*NewIDMap)[*iddata] = sIndx;

    //    dbgmsg << "  Inserted particle position " << sIndx << " = ";
    //    dbgmsg << *data << " with ID = " << prevloc << endl;
  }
}


/***************************************************************************
 * $RCSfile: ACLVISPtclBaseDataSource.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:24 $
 * IPPL_VERSION_ID: $Id: ACLVISPtclBaseDataSource.cpp,v 1.1.1.1 2003/01/23 07:40:24 adelmann Exp $ 
 ***************************************************************************/
